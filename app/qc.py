from __future__ import annotations

import os
import ast
import json
import re
import shutil
from pathlib import Path
from typing import List, Optional, TYPE_CHECKING

import pandas as pd

from .utils import log, log_err, run_cmd

if TYPE_CHECKING:
    from .entities import Config, BioProject, Sample, Dataset


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _safe_link_or_copy(src: Path, dst: Path, log_path: Path) -> None:
    """Hard-link if possible; fall back to symlink/copy. Always logs action."""
    try:
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        os.link(src, dst)
        log(f"[MultiQC] Hard-linked {src} -> {dst}", log_path)
        return
    except Exception as e:
        try:
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            os.symlink(src, dst)
            log(f"[MultiQC] Symlinked {src} -> {dst}", log_path)
            return
        except Exception as e2:
            shutil.copy2(src, dst)
            log(
                f"[MultiQC] Copied {src} -> {dst} "
                f"(hardlink/symlink failed: {e}; {e2})",
                log_path,
            )


# ─────────────────────────────────────────────────────────────────────────────
# MultiQC preparation (sanitized inputs so sample IDs are clean)
# ─────────────────────────────────────────────────────────────────────────────

def prepare_mqc_inputs_for_bp(
    bp: "BioProject",
    cfg: "Config",
) -> Optional[Path]:
    """
    Create <bp.path>/_mqc_inputs/<RUN> with:
      - run_info.json
      - abundance.tsv (and abundance.h5 if present)
      - <RUN>.fastp.json
      - kallisto log: per-sample log at Sample.log_path

    Returns the path, or None if nothing was prepared.

    Uses only the objects:
      • bp.path, bp.samples
      • s.outdir, s.id, s.log_path
    """
    bioproject_dir = bp.path.resolve()
    log_path = bp.log_path
    error_warnings = cfg.error_warnings

    tmp_root = bioproject_dir / "_mqc_inputs"

    if tmp_root.exists():
        shutil.rmtree(tmp_root)
    tmp_root.mkdir(parents=True, exist_ok=True)

    n_fastp = 0
    n_kallisto = 0

    # Each Sample in this BioProject is one SRR run
    for s in sorted(bp.samples, key=lambda x: x.id):
        run_id = s.id
        run_dir = s.outdir
        dst = tmp_root / run_id
        dst.mkdir(parents=True, exist_ok=True)

        # ------------------ fastp json ------------------
        # Prefer <RUN>.fastp.json; otherwise first *.fastp.json
        fjson = run_dir / "fastp.json"
        if not fjson.exists():
            hits = list(run_dir.glob("*fastp.json"))
            fjson = hits[0] if hits else None

        if fjson and fjson.exists():
            _safe_link_or_copy(fjson, dst / str(f"{run_id}.{fjson.name}"), log_path)
            n_fastp += 1
        else:
            log_err(
                error_warnings,
                log_path,
                f"[{run_id}] Missing fastp JSON; fastp totals may be absent in MultiQC",
            )

        # ------------------ kallisto artifacts ------------------
        have_k = False

        # run_info.json
        runinfo = run_dir / "run_info.json"
        if runinfo.exists():
            _safe_link_or_copy(runinfo, dst / "run_info.json", log_path)
            have_k = True

        # abundance.{tsv,h5}
        for name in ("abundance.tsv", "abundance.h5"):
            src = run_dir / name
            if src.exists():
                _safe_link_or_copy(src, dst / name, log_path)
                have_k = True

        # kallisto log: **only** the per-sample log
        klog_path = s.log_path  # <run_dir>/<run_id>_log.txt
        if klog_path.exists():
            _safe_link_or_copy(klog_path, dst / klog_path.name, log_path)
            have_k = True
        else:
            log_err(
                error_warnings,
                log_path,
                f"[{run_id}] Missing per-sample kallisto log: {klog_path.name}",
            )

        if have_k:
            n_kallisto += 1
        else:
            log_err(
                error_warnings,
                log_path,
                f"[{run_id}] Missing kallisto run_info.json and abundance file",
            )

    if (n_fastp + n_kallisto) == 0:
        log_err(
            error_warnings,
            log_path,
            f"[MultiQC sanitize] No inputs prepared in {bioproject_dir}",
        )
        shutil.rmtree(tmp_root, ignore_errors=True)
        return None

    log(
        f"[MultiQC sanitize] Prepared fastp={n_fastp}, kallisto={n_kallisto} "
        f"under {tmp_root}",
        log_path,
    )
    return tmp_root


# ─────────────────────────────────────────────────────────────────────────────
# MultiQC-derived read metrics (strict with fallback)
# ─────────────────────────────────────────────────────────────────────────────

def build_bp_metrics(
    bp: "BioProject",
    cfg: "Config",
    out_tsv: Optional[Path] = None,
) -> pd.DataFrame:
    """
    Build a per-SRR table for this BioProject:
        Sample, total_reads, high_quality_reads, pseudoaligned_reads

    Behaviour matches the original function; only the interface changed to
    use BioProject + Config instead of raw paths.
    """
    log_path = bp.log_path
    error_warnings = cfg.error_warnings

    def _log(msg: str) -> None:
        try:
            log(msg, log_path)
        except Exception:
            pass

    def _log_err(msg: str) -> None:
        try:
            log_err(error_warnings, log_path, msg)
        except Exception:
            pass

    def _srr(name: str) -> str:
        m = re.search(r"(SRR\d+)", str(name))
        return m.group(1) if m else str(name)

    def _to_int_or_na(x):
        try:
            if x is None or (isinstance(x, float) and pd.isna(x)):
                return pd.NA
            return int(round(float(x)))
        except Exception:
            return pd.NA

    bp_id = str(bp.id)
    bioproject_dir = bp.path
    mqc_dir = bioproject_dir / f"multiqc_{bp_id}_data"
    fastp_path = mqc_dir / "multiqc_fastp.txt"

    # ---------- FASTP ----------
    if not fastp_path.exists():
        _log_err(f"[{bp_id}] Missing MultiQC fastp file: {fastp_path}")
        return pd.DataFrame(
            columns=["Sample", "total_reads", "high_quality_reads", "pseudoaligned_reads"]
        )

    try:
        fastp = pd.read_csv(fastp_path, sep="\t", comment="#")
        # First column is the MultiQC sample id
        fastp = fastp.rename(columns={fastp.columns[0]: "Sample"})
        if "summary" not in fastp.columns:
            raise ValueError("'summary' column missing in MultiQC fastp table")

        # Keep only rows that actually have fastp payload (avoid stray *_1 lines)
        fastp = fastp[fastp["summary"].notna()].copy()

        # Normalize names to SRR########
        fastp["Sample"] = fastp["Sample"].astype(str).map(_srr)

        # Deduplicate after normalization
        fastp = fastp.drop_duplicates(subset=["Sample"], keep="first")

        # Extract totals from the 'summary' dict (stringified JSON or Python dict)
        def _parse_summary(cell):
            if isinstance(cell, dict):
                d = cell
            else:
                s = str(cell)
                try:
                    d = json.loads(s)
                except Exception:
                    d = ast.literal_eval(s)
            bf = d.get("before_filtering", {}).get("total_reads")
            af = d.get("after_filtering", {}).get("total_reads")
            return _to_int_or_na(bf), _to_int_or_na(af)

        totals = fastp["summary"].apply(_parse_summary)
        fastp["total_reads"] = totals.apply(lambda x: x[0]).astype("Int64")
        fastp["high_quality_reads"] = totals.apply(lambda x: x[1]).astype("Int64")

        f_sel = fastp[["Sample", "total_reads", "high_quality_reads"]].copy()
    except Exception as e:
        _log_err(f"[{bp_id}] Failed reading/parsing fastp table: {e}")
        return pd.DataFrame(
            columns=["Sample", "total_reads", "high_quality_reads", "pseudoaligned_reads"]
        )

    # ---------- KALLISTO from run_info.json ONLY ----------
    def _harvest_runinfo(root: Path) -> List[tuple[str, int, int]]:
        # returns (Sample, n_pseudoaligned, n_processed)
        rows: List[tuple[str, int, int]] = []
        if not root.exists():
            return rows
        for p in sorted(root.iterdir()):
            if not p.is_dir():
                continue
            run_id = _srr(p.name)
            ri = p / "run_info.json"
            if not ri.exists():
                continue
            try:
                with open(ri) as fh:
                    info = json.load(fh)
                npa = _to_int_or_na(info.get("n_pseudoaligned"))
                npr = _to_int_or_na(info.get("n_processed"))
                rows.append((run_id, npa, npr))
            except Exception as e:
                _log_err(f"[{bp_id}] Failed parsing {ri}: {e}")
        return rows

    rows: List[tuple[str, int, int]] = []
    rows += _harvest_runinfo(bioproject_dir / "_mqc_inputs")  # sanitized inputs
    rows += _harvest_runinfo(bioproject_dir)                  # raw SRR dirs (fallback)

    if rows:
        k_df = pd.DataFrame(rows, columns=["Sample", "pseudoaligned_reads", "_n_processed"])
        k_df["Sample"] = k_df["Sample"].astype(str).map(_srr)
        k_df = k_df.drop_duplicates(subset=["Sample"], keep="first")
    else:
        k_df = pd.DataFrame(columns=["Sample", "pseudoaligned_reads", "_n_processed"])

    # ---------- MERGE ----------
    try:
        out = (
            pd.merge(f_sel, k_df, on="Sample", how="left")
            .sort_values("Sample")
            .reset_index(drop=True)
        )
    except Exception as e:
        _log_err(f"[{bp_id}] Failed to assemble metrics table: {e}")
        return pd.DataFrame(
            columns=["Sample", "total_reads", "high_quality_reads", "pseudoaligned_reads"]
        )

    # ---------- ADJUST COUNTS TO FRAGMENTS WHEN DEFINITELY PAIRED ----------
    try:
        tr = out["total_reads"].astype("Int64")
        hq = out["high_quality_reads"].astype("Int64")
        np_ = out["_n_processed"].astype("Int64")

        # Vectorized paired/single detection – same logic as before
        paired_mask = (
            hq.notna()
            & np_.notna()
            & (abs(hq - (2 * np_)) <= (0.05 * (2 * np_.astype("float")))).fillna(False)
        )

        # Halve fastp counts only for confidently paired rows
        idx = paired_mask[paired_mask].index
        if len(idx) > 0:
            out.loc[idx, "total_reads"] = (
                out.loc[idx, "total_reads"].astype("float") / 2.0
            ).round().astype("Int64")
            out.loc[idx, "high_quality_reads"] = (
                out.loc[idx, "high_quality_reads"].astype("float") / 2.0
            ).round().astype("Int64")
    except Exception as e:
        _log_err(f"[{bp_id}] Failed while adjusting paired-end counts: {e}")

    # Drop helper (kallisto processed) before writing
    if "_n_processed" in out.columns:
        out = out.drop(columns=["_n_processed"])

    # Ensure dtypes are Int64
    for col in ("total_reads", "high_quality_reads", "pseudoaligned_reads"):
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce").astype("Int64")

    # ---------- WRITE ----------
    if out_tsv is None:
        out_tsv = bioproject_dir / "read_metrics.tsv"
    try:
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(out_tsv, sep="\t", index=False)
        _log(f"[{bp_id}] Read metrics written: {out_tsv}")
    except Exception as e:
        _log_err(f"[{bp_id}] Failed to write metrics TSV {out_tsv}: {e}")

    return out


# ─────────────────────────────────────────────────────────────────────────────
# Per-BioProject MultiQC (sanitized, using _mqc_inputs)
# ─────────────────────────────────────────────────────────────────────────────

def run_multiqc(
    bp: "BioProject",
    cfg: "Config",
    modules=("kallisto", "fastp"),
) -> Optional[Path]:
    """
    Per-BioProject MultiQC runner.

    This is the old `run_multiqc_sanitized_bp`, but now driven by objects:
      • `bp` (BioProject with samples)
      • `cfg` (Config with log + error_warnings)
    """
    log_path = bp.log_path
    error_warnings = cfg.error_warnings

    bioproject_dir = bp.path.resolve()
    inputs_root = prepare_mqc_inputs_for_bp(bp, cfg)
    if inputs_root is None:
        return None

    # Pick SRR dirs that contain run_info.json and abundance.tsv
    sra_dirs = []
    missing_msgs = []
    for d in sorted(p for p in inputs_root.iterdir() if p.is_dir()):
        has_runinfo = (d / "run_info.json").exists()
        has_tsv = (d / "abundance.tsv").exists()
        if has_runinfo and has_tsv:
            sra_dirs.append(d)
        else:
            reasons = []
            if not has_runinfo:
                reasons.append("run_info.json")
            if not has_tsv:
                reasons.append("abundance.tsv")
            missing_msgs.append(f"{d.name}: missing {', '.join(reasons)}")

    if not sra_dirs:
        log_err(
            error_warnings,
            log_path,
            "[MultiQC sanitize] No complete SRR dirs under "
            f"{inputs_root}.\n"
            + ("\n".join(missing_msgs) if missing_msgs else ""),
        )
        return None

    log(
        f"[MultiQC sanitize] Using {len(sra_dirs)} SRR dirs for parsing",
        log_path,
    )

    report_name = f"multiqc_{bioproject_dir.name}"
    cmd = [
        "multiqc",
        *map(str, sra_dirs),          # pass directories, not individual files
        "-o",
        str(bioproject_dir),
        "-n",
        report_name,
        "--force",
        "--ignore",
        "multiqc_*",
        "--ignore",
        "*/multiqc_*",
        "--ignore",
        "*.h5",                       # stop kallisto module from touching HDF5
        "-v",
    ]
    for m in modules:
        cmd += ["-m", m]

    # Run from _mqc_inputs (not strictly required, but keeps paths short)
    run_cmd(cmd, cwd=inputs_root, log_path=log_path)

    out = bioproject_dir / f"{report_name}_data"
    if not out.exists():
        log_err(error_warnings, log_path, f"[MultiQC] Expected data dir not found: {out}")
        return None
    return out


# ─────────────────────────────────────────────────────────────────────────────
# Global / generic MultiQC runner (for shared/global report)
# ─────────────────────────────────────────────────────────────────────────────

def run_multiqc_global(
    in_dir: Path,
    out_dir: Path,
    report_name: str,
    log_path: Path,
    modules=("kallisto", "fastp"),
) -> None:
    """
    Generic MultiQC runner for global/shared reports.

    This keeps the old API so the top-level pipeline can still call it
    path-based for the global MultiQC.
    """
    in_dir = Path(in_dir).resolve()
    out_dir = Path(out_dir)
    if not out_dir.is_absolute():
        out_dir = (in_dir / out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "multiqc",
        ".",
        "-o",
        str(out_dir),
        "-n",
        report_name,
        "--force",
        "--ignore",
        "multiqc_*",
        "--ignore",
        "*/multiqc_*",
    ]
    for m in modules:
        cmd += ["-m", m]

    run_cmd(cmd, cwd=in_dir, log_path=log_path)
