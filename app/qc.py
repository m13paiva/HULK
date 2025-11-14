import os
import ast
import json
import re
import pandas as pd
from pathlib import Path

from .utils import (
    log, log_err,
    run_cmd
)

# ─────────────────────────────────────────────────────────────────────────────
# MultiQC runner
# ─────────────────────────────────────────────────────────────────────────────

def run_multiqc(
    bioproject_dir: Path,
    log_path: Path,
    error_warnings: list[str],
    modules=('kallisto', 'fastp'),
) -> Path | None:
    """
    Run MultiQC directly on the SRR/sample directories inside a BioProject,
    without creating a _mqc_inputs/ staging area.

    We assume each sample directory may contain:
      - fastp.json or <RUN>.fastp.json (for fastp module)
      - run_info.json and/or kallisto logs (for kallisto module)
    """
    bioproject_dir = bioproject_dir.resolve()

    # Candidate sample dirs: immediate subdirs that are not MultiQC outputs
    sra_dirs = sorted(
        p for p in bioproject_dir.iterdir()
        if p.is_dir() and not p.name.startswith("multiqc_")
    )
    if not sra_dirs:
        log_err(
            error_warnings,
            log_path,
            f"[MultiQC] No sample directories found under {bioproject_dir}"
        )
        return None

    # Quick scan: how many dirs have fastp / kallisto inputs?
    n_fastp = 0
    n_kallisto = 0

    for d in sra_dirs:
        run_id = d.name

        # ---- fastp candidates ----
        # try <RUN>.fastp.json, fastp.json, then any *.fastp.json
        fjson_candidates = [
            d / f"{run_id}.fastp.json",
            d / "fastp.json",
        ]
        hits = [p for p in fjson_candidates if p.exists()]
        if not hits:
            hits = list(d.glob("*.fastp.json"))
        if hits:
            n_fastp += 1

        # ---- kallisto candidates ----
        has_kallisto = False

        # These are useful for metrics and some MultiQC behaviours
        if (d / "run_info.json").exists() or (d / "abundance.tsv").exists():
            has_kallisto = True

        # And any log that likely contains kallisto stdout
        if not has_kallisto:
            log_hits = list(d.glob("kallisto*.log")) + list(d.glob("*_log.txt"))
            if log_hits:
                has_kallisto = True

        if has_kallisto:
            n_kallisto += 1

    if n_fastp == 0 and n_kallisto == 0:
        log_err(
            error_warnings,
            log_path,
            f"[MultiQC] No fastp or kallisto inputs detected under {bioproject_dir}"
        )
        return None

    log(
        f"[MultiQC] Detected fastp={n_fastp}, kallisto={n_kallisto} candidate dirs under {bioproject_dir}",
        log_path,
    )

    report_name = f"multiqc_{bioproject_dir.name}"
    cmd = [
        "multiqc",
        *map(str, sra_dirs),
        "-o", str(bioproject_dir),
        "-n", report_name,
        "--force",
        "--ignore", "multiqc_*",
        "--ignore", "*/multiqc_*",
        "--ignore", "*.h5",  # keep kallisto module focused on text/logs
        "-v",
    ]
    for m in modules:
        cmd += ["-m", m]

    # Run from the BioProject directory
    run_cmd(cmd, cwd=bioproject_dir, log_path=log_path)

    out = bioproject_dir / f"{report_name}_data"
    if not out.exists():
        log_err(error_warnings, log_path, f"[MultiQC] Expected data dir not found: {out}")
        return None
    return out


def run_multiqc_global(
    in_dir: Path,
    out_dir: Path,
    report_name: str,
    log_path: Path,
    modules=('kallisto', 'fastp')
) -> None:
    """Generic MultiQC runner (used for the global/shared report)."""
    in_dir = in_dir.resolve()
    if not out_dir.is_absolute():
        out_dir = (in_dir / out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "multiqc", ".",
        "-o", str(out_dir),
        "-n", report_name,
        "--force",
        "--ignore", "multiqc_*",
        "--ignore", "*/multiqc_*",
        "--ignore", "*.h5",
    ]
    for m in modules:
        cmd += ["-m", m]
    run_cmd(cmd, cwd=in_dir, log_path=log_path)


# ─────────────────────────────────────────────────────────────────────────────
# MultiQC-derived read metrics (strict with fallback)
# ─────────────────────────────────────────────────────────────────────────────

def build_bp_metrics(
    bp_id: str,
    bioproject_dir: Path,
    log_path: Path,
    error_warnings: list[str],
    out_tsv: Path | None = None,
) -> pd.DataFrame:
    """
    Build a per-SRR table:
        Sample, total_reads, high_quality_reads, pseudoaligned_reads

    Rules:
    - fastp totals come from MultiQC fastp table (summary dict).
    - pseudoaligned_reads taken ONLY from run_info.json (kallisto).
    - Paired/single is inferred from kallisto's n_processed vs fastp high_quality:
        * if high_quality ~= 2 * n_processed  -> paired-end  -> divide fastp counts by 2
        * if high_quality ~= 1 * n_processed  -> single-end  -> leave as-is
      (5% tolerance; if unknown/ambiguous, leave counts as-is.)
    """

    def _log(msg: str):
        try:
            log(msg, log_path)
        except Exception:
            pass

    def _log_err(msg: str):
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

    bp_id = str(bp_id)
    bioproject_dir = Path(bioproject_dir)
    mqc_dir = bioproject_dir / f"multiqc_{bp_id}_data"
    fastp_path = mqc_dir / "multiqc_fastp.txt"

    # ---------- FASTP ----------
    if not fastp_path.exists():
        _log_err(f"[{bp_id}] Missing MultiQC fastp file: {fastp_path}")
        return pd.DataFrame(columns=["Sample","total_reads","high_quality_reads","pseudoaligned_reads"])

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
            af = d.get("after_filtering",  {}).get("total_reads")
            return _to_int_or_na(bf), _to_int_or_na(af)

        totals = fastp["summary"].apply(_parse_summary)
        fastp["total_reads"]        = totals.apply(lambda x: x[0]).astype("Int64")
        fastp["high_quality_reads"] = totals.apply(lambda x: x[1]).astype("Int64")

        f_sel = fastp[["Sample", "total_reads", "high_quality_reads"]].copy()
    except Exception as e:
        _log_err(f"[{bp_id}] Failed reading/parsing fastp table: {e}")
        return pd.DataFrame(columns=["Sample","total_reads","high_quality_reads","pseudoaligned_reads"])


    # ---------- KALLISTO from run_info.json ONLY ----------
    def _harvest_runinfo(root: Path) -> list[tuple[str, int, int]]:
        # returns (Sample, n_pseudoaligned, n_processed)
        rows: list[tuple[str,int,int]] = []
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

    # Now we only look in the BioProject SRR dirs (no _mqc_inputs anymore)
    rows = _harvest_runinfo(bioproject_dir)

    if rows:
        k_df = pd.DataFrame(rows, columns=["Sample","pseudoaligned_reads","_n_processed"])
        k_df["Sample"] = k_df["Sample"].astype(str).map(_srr)
        k_df = k_df.drop_duplicates(subset=["Sample"], keep="first")
    else:
        k_df = pd.DataFrame(columns=["Sample","pseudoaligned_reads","_n_processed"])

    # ---------- MERGE ----------
    try:
        out = (
            pd.merge(f_sel, k_df, on="Sample", how="left")
              .sort_values("Sample")
              .reset_index(drop=True)
        )
    except Exception as e:
        _log_err(f"[{bp_id}] Failed to assemble metrics table: {e}")
        return pd.DataFrame(columns=["Sample","total_reads","high_quality_reads","pseudoaligned_reads"])

    # ---------- ADJUST COUNTS TO FRAGMENTS WHEN DEFINITELY PAIRED ----------
    try:
        tr = out["total_reads"].astype("Int64")
        hq = out["high_quality_reads"].astype("Int64")
        np_ = out["_n_processed"].astype("Int64")

        # paired if hq ~= 2 * n_processed (5% tolerance)
        paired_mask  = hq.notna() & np_.notna()
        paired_mask  = paired_mask & (
            (abs(hq - (2 * np_)) <= (0.05 * (2 * np_.astype("float"))))
        ).fillna(False)

        # single if hq ~= 1 * n_processed (5% tolerance)
        single_mask  = hq.notna() & np_.notna()
        single_mask  = single_mask & (
            (abs(hq - np_) <= (0.05 * np_.astype("float")))
        ).fillna(False)

        # Halve fastp counts only for confidently paired rows
        idx = paired_mask[paired_mask].index
        if len(idx) > 0:
            out.loc[idx, "total_reads"]        = (out.loc[idx, "total_reads"].astype("float") / 2.0).round().astype("Int64")
            out.loc[idx, "high_quality_reads"] = (out.loc[idx, "high_quality_reads"].astype("float") / 2.0).round().astype("Int64")
    except Exception as e:
        _log_err(f"[{bp_id}] Failed while adjusting paired-end counts: {e}")

    # Drop helper column before writing
    if "_n_processed" in out.columns:
        out = out.drop(columns=["_n_processed"])

    # Ensure dtypes are Int64
    for col in ("total_reads","high_quality_reads","pseudoaligned_reads"):
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
