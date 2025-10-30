# utils.py

import os
import math
import subprocess
from pathlib import Path

import pandas as pd


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
    Renames _1.fastq files to .fastq for single-end data.
    """
    cands = [
        (outdir / f"{run_id}_1.fastq", outdir / f"{run_id}_2.fastq"),
        (outdir / f"{run_id}_1.fastq.gz", outdir / f"{run_id}_2.fastq.gz"),
    ]

    # First check for paired-end
    for r1, r2 in cands:
        if r1.exists() and r2.exists():
            return "PAIRED", r1, r2

    # Then check for single-end files
    for se in (outdir / f"{run_id}.fastq", outdir / f"{run_id}.fastq.gz",
               outdir / f"{run_id}_1.fastq", outdir / f"{run_id}_1.fastq.gz"):
        if se.exists():
            # If file has _1 suffix, rename it to remove the suffix
            if "_1.fastq" in se.name:
                new_name = se.name.replace("_1.fastq", ".fastq")
                new_path = outdir / new_name
                se.rename(new_path)
                return "SINGLE", new_path, None
            else:
                return "SINGLE", se, None

    return "MISSING", None, None


# ─────────────────────────────────────────────────────────────────────────────
# TPM merges
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





