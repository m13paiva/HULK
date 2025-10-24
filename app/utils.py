# utils.py
import os
import sys
import math
import time
import shutil
import subprocess
import ast
import json
from pathlib import Path

import pandas as pd
from tqdm.auto import tqdm


# ─────────────────────────────────────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────────────────────────────────────

def log(msg: str, log_path: Path) -> None:
    """Append a single line to the shared log file."""
    with open(log_path, "a", buffering=1) as f:
        f.write(str(msg).rstrip() + "\n")

def log_err(error_warnings: list[str], log_path: Path, msg: str) -> None:
    """Log an error/warning and also append it to `error_warnings`."""
    log(msg, log_path)
    error_warnings.append(str(msg))


# ─────────────────────────────────────────────────────────────────────────────
# Subprocess helpers
# ─────────────────────────────────────────────────────────────────────────────

def run_cmd(cmd: list[str], cwd: Path | None, log_path: Path) -> None:
    """
    Run a command and tee stdout/stderr to the main log file.
    Raises CalledProcessError on non-zero exit.
    """
    with open(log_path, "a", buffering=1) as f:
        f.write(f"## cwd: {cwd}\n")
        f.write(">> " + " ".join(map(str, cmd)) + "\n")
        f.flush()
        subprocess.run(cmd, cwd=cwd, stdout=f, stderr=f, check=True)

def run_cmd_stream(cmd: list[str], cwd: Path | None, log_path: Path, side_log_path: Path | None = None) -> None:
    """
    Like `run_cmd`, but streams line-by-line and optionally duplicates output
    into a per-SRA side log (e.g., kallisto_<RUN>.log).
    """
    with open(log_path, "a", buffering=1) as flog:
        flog.write(f"## cwd: {cwd}\n")
        flog.write(">> " + " ".join(map(str, cmd)) + "\n")
        flog.flush()

        proc = subprocess.Popen(
            cmd, cwd=cwd,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1
        )
        side = open(side_log_path, "w", buffering=1) if side_log_path is not None else None
        try:
            assert proc.stdout is not None
            for line in proc.stdout:
                flog.write(line)
                if side is not None:
                    side.write(line)
        finally:
            if side is not None:
                side.close()

        ret = proc.wait()
        if ret != 0:
            raise subprocess.CalledProcessError(ret, cmd)


# ─────────────────────────────────────────────────────────────────────────────
# Filesystem helpers
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
            log(f"[MultiQC] Copied {src} -> {dst} (hardlink/symlink failed: {e}; {e2})", log_path)

def clean_fastq_files(run_dir: Path, recursive: bool = False) -> None:
    """Delete any files whose suffix includes 'fastq' (space saver)."""
    files = run_dir.rglob("*") if recursive else run_dir.glob("*")
    for f in files:
        if f.is_file() and "fastq" in "".join(f.suffixes).lower():
            try:
                f.unlink()
            except Exception:
                pass


# ─────────────────────────────────────────────────────────────────────────────
# Dataframe / planning helpers
# ─────────────────────────────────────────────────────────────────────────────

def df_to_dict(df: pd.DataFrame) -> dict[str, list[list[str]]]:
    """
    Convert input DataFrame into:
      { BioProject: [[Run, Model], ...] }
    Requires columns: 'Run', 'BioProject', 'Model'.
    """
    return {
        bp_id: df.loc[df['BioProject'] == bp_id, ['Run', 'Model']].to_numpy().tolist()
        for bp_id in df['BioProject'].unique()
    }

def get_available_threads() -> int:
    """Detect usable CPU count (Docker/CGroups aware)."""
    try:
        with open("/sys/fs/cgroup/cpu.max") as f:
            quota, period = f.read().split()
            if quota != "max":
                return max(1, math.ceil(int(quota) / int(period)))
    except Exception:
        pass
    return os.cpu_count() or 1

def plan_workers(total_threads: int, n_sras: int, min_threads: int = 4) -> tuple[int, int]:
    """
    Decide how many parallel jobs and threads per job to use.
    Returns (jobs, threads_each).
    """
    threads_each = max(1, min(min_threads, total_threads))
    jobs = max(1, min(n_sras, total_threads // threads_each))
    if jobs * threads_each < total_threads and jobs < n_sras:
        jobs += 1
    threads_each = max(1, total_threads // jobs)
    return jobs, threads_each


# ─────────────────────────────────────────────────────────────────────────────
# FASTQ / status helpers
# ─────────────────────────────────────────────────────────────────────────────

def transcriptome_suffixes(p: Path) -> str:
    """Return normalized transcriptome 'suffix' (handles .gz)."""
    suff = [s.lower() for s in p.suffixes]
    if not suff:
        return ""
    return "".join(suff[-2:]) if suff[-1] == ".gz" else suff[-1]

def is_sra_done(run_dir: Path) -> bool:
    """An SRA is considered done if kallisto's abundance.tsv exists."""
    #return False
    return (run_dir / "abundance.tsv").exists()

def detect_fastq_layout(run_id: str, outdir: Path):
    """
    Detect SINGLE vs PAIRED FASTQs produced by fasterq-dump.
    Checks both .fastq and .fastq.gz.
    """
    cands = [
        (outdir / f"{run_id}_1.fastq", outdir / f"{run_id}_2.fastq"),
        (outdir / f"{run_id}_1.fastq.gz", outdir / f"{run_id}_2.fastq.gz"),
    ]
    for r1, r2 in cands:
        if r1.exists() and r2.exists():
            return "PAIRED", r1, r2
    for se in (outdir / f"{run_id}.fastq", outdir / f"{run_id}.fastq.gz"):
        if se.exists():
            return "SINGLE", se, None
    return "MISSING", None, None

def pick_fastq(run_dir: Path, stem: str, log_path: Path, error_warnings: list[str]) -> Path | None:
    """
    Prefer trimmed files when present; otherwise fall back to raw and log once.
    """
    for suff in (".trim.fastq.gz", ".trim.fastq"):
        p = run_dir / f"{stem}{suff}"
        if p.exists():
            return p
    for suff in (".fastq.gz", ".fastq"):
        p = run_dir / f"{stem}{suff}"
        if p.exists():
            log_err(error_warnings, log_path,
                    f"[{stem}] No trimmed FASTQs found in {run_dir}; quantifying on untrimmed reads.")
            return p
    return None


# ─────────────────────────────────────────────────────────────────────────────
# MultiQC preparation (sanitized inputs so sample IDs are clean)
# ─────────────────────────────────────────────────────────────────────────────

def prepare_mqc_inputs_for_bp(
    bioproject_dir: Path,
    log_path: Path,
    error_warnings: list[str],
) -> Path | None:
    """
    Create <bioproject_dir>/_mqc_inputs/<RUN> with:
      - run_info.json
      - abundance.tsv (and abundance.h5 if present)
      - <RUN>.fastp.json
    Returns the path, or None if nothing was prepared.
    """
    tmp_root = bioproject_dir / "_mqc_inputs"
    if tmp_root.exists():
        shutil.rmtree(tmp_root)
    tmp_root.mkdir(parents=True, exist_ok=True)

    n_fastp = 0
    n_kallisto = 0

    for sra_dir in sorted(p for p in bioproject_dir.iterdir() if p.is_dir()):
        run_id = sra_dir.name
        dst = tmp_root / run_id
        dst.mkdir(parents=True, exist_ok=True)

        # fastp json (prefer exact; otherwise any *.fastp.json)
        fjson = sra_dir / f"{run_id}.fastp.json"
        if not fjson.exists():
            hits = list(sra_dir.glob("*.fastp.json"))
            fjson = hits[0] if hits else None
        if fjson and fjson.exists():
            _safe_link_or_copy(fjson, dst / fjson.name, log_path)
            n_fastp += 1
        else:
            log_err(error_warnings, log_path, f"[{run_id}] Missing fastp JSON; fastp totals may be absent in MultiQC")

        # kallisto: expose stdout log (what MultiQC actually parses) + run_info/abundance
        have_k = False

        runinfo = sra_dir / "run_info.json"
        if runinfo.exists():
            _safe_link_or_copy(runinfo, dst / "run_info.json", log_path)
            have_k = True

        for name in ("abundance.tsv", "abundance.h5"):
            src = sra_dir / name
            if src.exists():
                _safe_link_or_copy(src, dst / name, log_path)
                have_k = True

        # Also copy the kallisto stdout log so MultiQC can detect the module
        # Your pipeline writes: outdir / f"kallisto_{run_id}.log"
        klog = sra_dir / f"kallisto_{run_id}.log"
        if not klog.exists():
            # Fallback: grab any plausible kallisto log in the SRA dir
            hits = list(sra_dir.glob("kallisto*.log"))
            if hits:
                klog = hits[0]
        if klog.exists():
            _safe_link_or_copy(klog, dst / klog.name, log_path)
            have_k = True

        if have_k:
            n_kallisto += 1

        else:
            log_err(error_warnings, log_path, f"[{run_id}] Missing kallisto run_info.json and abundance file")

    if (n_fastp + n_kallisto) == 0:
        log_err(error_warnings, log_path, f"[MultiQC sanitize] No inputs prepared in {bioproject_dir}")
        shutil.rmtree(tmp_root, ignore_errors=True)
        return None

    log(f"[MultiQC sanitize] Prepared fastp={n_fastp}, kallisto={n_kallisto} under {tmp_root}", log_path)
    return tmp_root


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
    import ast, json, re
    import pandas as pd

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

    def _close(a, b, rtol=0.05):
        # relative closeness with tolerance
        try:
            return (pd.notna(a) and pd.notna(b) and abs(int(a) - int(b)) <= rtol * max(1, int(b)))
        except Exception:
            return False

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

    rows = []
    rows += _harvest_runinfo(bioproject_dir / "_mqc_inputs")  # sanitized inputs
    rows += _harvest_runinfo(bioproject_dir)                  # raw SRR dirs (fallback)

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
    # Decide per-row using kallisto ground truth (5% tolerance).
    # If ambiguous or n_processed missing, leave as-is (assume single-end / already fragment-scale).
    try:
        # Work on local copies to avoid pandas warnings
        tr = out["total_reads"].astype("Int64")
        hq = out["high_quality_reads"].astype("Int64")
        np_ = out["_n_processed"].astype("Int64")

        # Vectorized paired/single detection
        paired_mask  = hq.notna() & np_.notna() & hq.apply(lambda x: False if pd.isna(x) else True) & np_.apply(lambda x: False if pd.isna(x) else True)
        paired_mask  = paired_mask & (abs(hq - (2 * np_)) <= (0.05 * (2 * np_.astype("float")))).fillna(False)

        single_mask  = hq.notna() & np_.notna()
        single_mask  = single_mask & (abs(hq - np_) <= (0.05 * np_.astype("float"))).fillna(False)

        # Halve fastp counts only for confidently paired rows
        idx = paired_mask[paired_mask].index
        if len(idx) > 0:
            out.loc[idx, "total_reads"]        = (out.loc[idx, "total_reads"].astype("float") / 2.0).round().astype("Int64")
            out.loc[idx, "high_quality_reads"] = (out.loc[idx, "high_quality_reads"].astype("float") / 2.0).round().astype("Int64")
    except Exception as e:
        _log_err(f"[{bp_id}] Failed while adjusting paired-end counts: {e}")

    # Drop helper (kallisto processed) before writing
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



# ─────────────────────────────────────────────────────────────────────────────
# TPM merges & tximport helpers
# ─────────────────────────────────────────────────────────────────────────────

def merge_bioproject_tpm(bioproject_dir: Path, log_path: Path, error_warnings: list[str], output_file: Path | None = None) -> None:
    """Merge all kallisto TPMs for one BP into a single TSV."""
    abundance_files = list(bioproject_dir.glob("*/abundance.tsv"))
    if not abundance_files:
        log_err(error_warnings, log_path, f"No abundance.tsv found in {bioproject_dir}")
        return

    merged = None
    for ab_file in abundance_files:
        run_id = ab_file.parent.name
        df = pd.read_csv(ab_file, sep="\t", usecols=["target_id", "tpm"]).rename(columns={"tpm": run_id})
        merged = df if merged is None else merged.merge(df, on="target_id", how="outer")

    merged = merged.fillna(0.0)
    if output_file is None:
        output_file = bioproject_dir / f"{bioproject_dir.name}_TPM.tsv"
    merged.to_csv(output_file, sep="\t", index=False)

def smash():
    '''Makes Hulk Smash (easter egg)'''
    VIDEO_PATH = os.environ.get("HULK_SMASH_PATH", "/opt/hulk/hulk_smash.mp4")

    subprocess.run(
        [
            "mpv",
            "--vo=caca",
            "--really-quiet",
            "--contrast=100",
            "--brightness=-100",
            VIDEO_PATH,
        ],
        check=False,
        stderr=subprocess.DEVNULL,
    )


def merge_all_bioprojects_tpm(output_root: Path, log_path: Path, error_warnings: list[str], output_file: Path | None = None) -> None:
    """Merge all kallisto TPMs across all BioProjects into <output_root>/all_bioprojects_TPM.tsv."""
    abundance_files = list(output_root.glob("*/*/abundance.tsv"))
    if not abundance_files:
        log_err(error_warnings, log_path, f"No abundance.tsv files found in {output_root}")
        return

    merged = None
    for ab_file in abundance_files:
        run_id = ab_file.parent.name
        df = pd.read_csv(ab_file, sep="\t", usecols=["target_id", "tpm"]).rename(columns={"tpm": run_id})
        merged = df if merged is None else merged.merge(df, on="target_id", how="outer")

    merged = merged.fillna(0.0)
    if output_file is None:
        output_file = output_root / "all_bioprojects_TPM.tsv"
    merged.to_csv(output_file, sep="\t", index=False)


def _gene_counts(abundance_files: list[Path], sample_names: list[str], tx2gene_path: Path,
                 mode=None,ignore_tx_version=False) -> pd.DataFrame:
    """Compute gene-level counts from kallisto outputs with pytximport (quiet mode)."""
    prev_tqdm = os.environ.get("TQDM_DISABLE")
    os.environ["TQDM_DISABLE"] = "1"
    try:
        from pytximport import tximport
        with open(os.devnull, "w") as devnull:
            sys_stdout, sys_stderr = sys.stdout, sys.stderr
            try:
                sys.stdout = devnull
                sys.stderr = devnull
                ds = tximport(
                    [str(p) for p in abundance_files],
                    data_type="kallisto",
                    transcript_gene_map=str(tx2gene_path),
                    counts_from_abundance=mode,
                    output_type="xarray",
                    return_data=True,
                    existence_optional=False,
                    ignore_transcript_version=ignore_tx_version,
                )
            finally:
                sys.stdout, sys.stderr = sys_stdout, sys_stderr
    finally:
        if prev_tqdm is None:
            del os.environ["TQDM_DISABLE"]
        else:
            os.environ["TQDM_DISABLE"] = prev_tqdm

    df = ds["counts"].to_pandas()
    df.columns = sample_names
    return df


def bp_gene_counts(bioproject_dir: Path, sra_ids: list[str], tx2gene: Path, log_path: Path, error_warnings: list[str],tximport_opts=None) -> Path | None:
    """Write <bioproject_dir>/gene_counts.tsv (genes × samples). Logs missing SRAs but proceeds."""
    files, names = [], []
    for run_id in sra_ids:
        f = bioproject_dir / run_id / "abundance.tsv"
        if f.exists():
            files.append(f); names.append(run_id)
        else:
            log_err(error_warnings, log_path, f"[tximport] Missing file: {f}")

    if not files:
        log_err(error_warnings, log_path, f"[tximport] No abundance.tsv found in {bioproject_dir}, skipping gene table")
        return None
    _opts = {k: v for k, v in (tximport_opts or {}).items()
             if k in {"mode", "ignore_tx_version"} and v is not None}
    gdf = _gene_counts(files, names, tx2gene,**_opts)
    out = bioproject_dir / "gene_counts.tsv"
    gdf.to_csv(out, sep="\t")
    log(f"🧬 Gene counts (pytximport) written: {out}", log_path)
    return out

def global_gene_counts(output_root: Path, tx2gene: Path, log_path: Path, error_warnings: list[str], output_file: Path,tximport_opts=None) -> None:
    """
    Build a global gene count table across all BPs using pytximport directly.
    """
    files = list(output_root.glob("*/*/abundance.tsv"))
    if not files:
        log_err(error_warnings, log_path, f"[tximport] No abundance.tsv found under {output_root}")
        return
    names = [p.parent.name for p in files]
    _opts = {k: v for k, v in (tximport_opts or {}).items()
             if k in {"mode", "ignore_tx_version"} and v is not None}
    gdf = _gene_counts(files, names, tx2gene,**_opts)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_csv(output_file, sep="\t")
    log(f"🧬 Global gene counts (pytximport) written: {output_file}", log_path)
