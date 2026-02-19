# utils.py

import os
import math
import subprocess
from pathlib import Path
from typing import List, Union, Iterable
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────────────────────────────────────
LogPathLike = Union[str, Path, Iterable[Union[str, Path]]]

def _normalize_log_paths(log_path: LogPathLike | None) -> list[Path]:
    if log_path is None:
        return []
    if isinstance(log_path, (str, Path)):
        return [Path(log_path)]
    # assume iterable of paths
    try:
        return [Path(p) for p in log_path]
    except TypeError:
        return [Path(log_path)]

def log(msg: str, log_path: LogPathLike | None) -> None:
    """Append a line to one or more log files."""
    paths = _normalize_log_paths(log_path)
    if not paths:
        return
    line = str(msg).rstrip() + "\n"
    for p in paths:
        try:
            p.parent.mkdir(parents=True, exist_ok=True)
            with open(p, "a", buffering=1) as f:
                f.write(line)
        except OSError:
            # disk full / permission errors: nothing else we can do safely
            continue

def log_err(error_warnings: list[str], log_path: LogPathLike | None, msg: str) -> None:
    """Record an error message in memory and logs."""
    error_warnings.append(str(msg))
    log(msg, log_path)


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

def scan_fastqs(directory: Path) -> List[Path]:
    exts = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    files = sorted(p for p in Path(directory).iterdir() if p.is_file() and p.name.lower().endswith(exts))
    return files


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

def pad_desc(name: str, width: int = 14) -> str:
    """Pad tqdm descriptions so bars line up neatly."""
    return name.ljust(width)


def generate_read_metrics_plot(dataset, out_dir: Path, log_path: Path) -> None:
    """
    Generates a bar plot of mean percentages for High Quality and Pseudoaligned reads
    per BioProject, aggregating the individual read metrics files.
    """
    log("Generating BioProject read metrics summary plot...", log_path)

    data_list = []

    # Iterate through the BioProjects in the Dataset
    for bp in dataset.bioprojects:
        # Construct path: BP_FOLDER / BP_ID_read_metrics.tsv
        # We assume the file naming convention matches the BP ID
        metric_file = bp.path / f"{bp.id}_read_metrics.tsv"

        if not metric_file.exists():
            # Try generic name if specific one fails
            metric_file = bp.path / "read_metrics.tsv"

        if not metric_file.exists():
            log(f"[Warn] Metrics file not found for {bp.id}, skipping.", log_path)
            continue

        try:
            df = pd.read_csv(metric_file, sep='\t')
        except Exception as e:
            log(f"[Error] Failed to read {metric_file}: {e}", log_path)
            continue

        # Validation
        required = {'total_reads', 'high_quality_reads', 'pseudoaligned_reads'}
        if not required.issubset(df.columns):
            log(f"[Warn] File {metric_file} missing required columns {required - set(df.columns)}", log_path)
            continue

        # Calculate Percentages
        # Avoid division by zero
        df = df[df['total_reads'] > 0].copy()
        df['High Quality'] = (df['high_quality_reads'] / df['total_reads']) * 100
        df['Pseudoaligned'] = (df['pseudoaligned_reads'] / df['total_reads']) * 100

        # Melt for Seaborn
        df_melted = df.melt(
            id_vars=['Sample'],
            value_vars=['High Quality', 'Pseudoaligned'],
            var_name='Metric',
            value_name='Percentage'
        )

        df_melted['Bioproject'] = bp.id
        data_list.append(df_melted)

    if not data_list:
        log("[Warn] No valid metrics data found to plot.", log_path)
        return

    final_df = pd.concat(data_list, ignore_index=True)

    # Save the source data for the plot
    plot_data_file = out_dir / "bioproject_mean_percentages.tsv"
    final_df.to_csv(plot_data_file, sep="\t", index=False)

    # Plotting
    try:
        plt.figure(figsize=(12, 8))  # Increased size slightly for readability

        sns.barplot(
            data=final_df,
            x='Bioproject',
            y='Percentage',
            hue='Metric',
            errorbar='sd',
            palette="viridis"  # Viridis is easier on the eyes than default
        )

        plt.title('Mean Percentage of High Quality and Pseudoaligned Reads per BioProject')
        plt.ylabel('Mean Percentage of Total Reads (%)')
        plt.xlabel('BioProject')

        # Rotation
        plt.xticks(rotation=90)
        plt.ylim(0, 105)
        plt.legend(title='Metric', loc='upper right')
        plt.tight_layout()

        output_file = out_dir / 'bioproject_mean_percentages.pdf'
        plt.savefig(output_file, format='pdf')
        plt.close()

        log(f"Saved metrics plot to {output_file}", log_path)

    except Exception as e:
        log(f"[Error] Failed during plotting: {e}", log_path)

def bp_seidr_batches(dataset: "Dataset", min_samples: int = 50) -> List[List["BioProject"]]:
    """
    Groups BioProjects into batches where each batch has at least min_samples.
    Maximizes the number of batches by greedy small-to-large pairing.
    """
    if dataset.mode != "SRR":
        # I'm assuming you aren't doing this for FASTQ mode since BioProjects don't exist there.
        # But knowing you, you might try.
        return []

    # 1. Separate the "big fish" from the "fry"
    big_bps = [bp for bp in dataset.bioprojects if bp.total() >= min_samples]
    small_bps = sorted(
        [bp for bp in dataset.bioprojects if bp.total() < min_samples],
        key=lambda x: x.total()
    )

    batches = []

    # 2. Every big BioProject gets its own batch. Easy.
    for bp in big_bps:
        batches.append([bp])

    # 3. Greedily pack the small ones to reach the threshold
    current_batch = []
    current_count = 0

    for bp in small_bps:
        current_batch.append(bp)
        current_count += bp.total()

        if current_count >= min_samples:
            batches.append(current_batch)
            current_batch = []
            current_count = 0

    # 4. Handle the leftovers.
    # If we have a final batch that didn't reach 50, we have to shove it
    # into the last successful batch to keep things valid.
    if current_batch:
        if batches:
            batches[-1].extend(current_batch)
        else:
            # If the entire dataset doesn't have 50 samples,
            # you get one sad, undersized batch.
            batches.append(current_batch)

    return batches

