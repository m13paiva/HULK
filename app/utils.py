import os
import sys
import math
import signal
import threading
import subprocess
from pathlib import Path
from typing import List, Union, Iterable, Optional, Dict, Any
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────────────
# Process Management & Interrupt Handling
# ─────────────────────────────────────────────────────────────────────────────
_active_processes_lock = threading.Lock()
_active_processes = set()

def register_process(proc: subprocess.Popen) -> None:
    """Registers an active subprocess so it can be killed on SIGINT / Ctrl+C."""
    with _active_processes_lock:
        _active_processes.add(proc)

def unregister_process(proc: subprocess.Popen) -> None:
    """Unregisters a completed subprocess."""
    with _active_processes_lock:
        _active_processes.discard(proc)

def kill_all_active_processes() -> None:
    """Force-terminates all tracked subprocesses and process groups immediately."""
    with _active_processes_lock:
        for proc in list(_active_processes):
            try:
                if proc.poll() is None:
                    if hasattr(os, 'killpg'):
                        try:
                            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                        except Exception:
                            proc.kill()
                    else:
                        proc.kill()
            except Exception:
                pass
        _active_processes.clear()

def setup_interrupt_handlers() -> None:
    """Configures SIGINT (Ctrl+C) and SIGTERM signal handlers to instantly terminate all child processes."""
    def _on_signal(signum, frame):
        print("\n[HULK] Signal received. Interrupting execution and killing processes...", file=sys.stderr, flush=True)
        kill_all_active_processes()
        os._exit(130)

    try:
        signal.signal(signal.SIGINT, _on_signal)
        signal.signal(signal.SIGTERM, _on_signal)
    except (ValueError, AttributeError):
        pass

def run_managed_subprocess(
    cmd: list[str],
    cwd: Path | str | None = None,
    stdout=None,
    stderr=None,
    env: dict[str, str] | None = None,
    check: bool = True
) -> subprocess.Popen:
    """
    Spawns a subprocess in its own process group, tracks it, and waits for completion.
    Guarantees child process termination if SIGINT / Ctrl+C is caught.
    """
    kwargs = {}
    if hasattr(os, 'setsid'):
        kwargs['start_new_session'] = True

    proc = subprocess.Popen(
        cmd, cwd=str(cwd) if cwd else None, stdout=stdout, stderr=stderr, env=env, text=True, **kwargs
    )
    register_process(proc)
    try:
        ret = proc.wait()
        if check and ret != 0:
            raise subprocess.CalledProcessError(ret, cmd)
        return proc
    finally:
        unregister_process(proc)


# ─────────────────────────────────────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────────────────────────────────────
LogPathLike = Union[str, Path, Iterable[Union[str, Path]]]

def _normalize_log_paths(log_path: LogPathLike | None) -> list[Path]:
    """
    Normalizes a given log path or collection of log paths into a list of Path objects.

    Args:
        log_path (LogPathLike | None): A single path, an iterable of paths, or None.

    Returns:
        list[Path]: A list of resolved Path objects.
    """
    if log_path is None:
        return []
    if isinstance(log_path, (str, Path)):
        return [Path(log_path)]
    try:
        return [Path(p) for p in log_path]
    except TypeError:
        return [Path(log_path)]

def log(msg: str, log_path: LogPathLike | None) -> None:
    """
    Appends a message line to one or more designated log files.

    Args:
        msg (str): The message to log.
        log_path (LogPathLike | None): The destination path(s) for the log entry.
    """
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
            continue


# ─────────────────────────────────────────────────────────────────────────────
# Subprocess helpers
# ─────────────────────────────────────────────────────────────────────────────

def run_cmd(cmd: list[str], cwd: Path | None, log_path: Path) -> None:
    """
    Runs a subprocess command and securely tees stdout/stderr to the main log file.

    Args:
        cmd (list[str]): Subprocess command sequences.
        cwd (Path | None): Explicit active directory defining context maps.
        log_path (Path): Direct target output path handling.
    """
    with open(log_path, "a", buffering=1) as f:
        f.write(f"## cwd: {cwd}\n")
        f.write(">> " + " ".join(map(str, cmd)) + "\n")
        f.flush()
        run_managed_subprocess(cmd, cwd=cwd, stdout=f, stderr=f, check=True)

def run_cmd_stream(cmd: list[str], cwd: Path | None, log_path: Path, side_log_path: Path | None = None) -> None:
    """
    Streams subprocess lines safely handling memory output blocks explicitly.

    Args:
        cmd (list[str]): Target hook values maps.
        cwd (Path | None): Working indicator variable path component limitations parameters.
        log_path (Path): General target file log output limit location component mappings limitations output variable components map variables indicators values array locations limit location mappings indicators components map location mappings variables properties target.
        side_log_path (Path | None, optional): Explicit target override. Defaults to None.
    """
    with open(log_path, "a", buffering=1) as flog:
        flog.write(f"## cwd: {cwd}\n")
        flog.write(">> " + " ".join(map(str, cmd)) + "\n")
        flog.flush()

        kwargs = {}
        if hasattr(os, 'setsid'):
            kwargs['start_new_session'] = True

        proc = subprocess.Popen(
            cmd, cwd=cwd,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1, **kwargs
        )
        register_process(proc)
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
            unregister_process(proc)

        ret = proc.wait()
        if ret != 0:
            raise subprocess.CalledProcessError(ret, cmd)


# ─────────────────────────────────────────────────────────────────────────────
# Filesystem helpers
# ─────────────────────────────────────────────────────────────────────────────

def clean_fastq_files(run_dir: Path, recursive: bool = False) -> None:
    """
    Unlinks residual target components mapping variable bounds.

    Args:
        run_dir (Path): Variable mappings limitations indicator.
        recursive (bool, optional): Variables properties target mapping property bounds limit target arrays parameters limitations limitation mappings array indicator properties limitations parameters indicators map boolean limitations limit target arrays property component mapping. Defaults to False.
    """
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
    Converts DataFrame component input maps efficiently.

    Args:
        df (pd.DataFrame): Data boundary matrix component maps limits targets variables.

    Returns:
        dict[str, list[list[str]]]: Map target definition mapping arrays components limitation.
    """
    return {
        bp_id: df.loc[df['BioProject'] == bp_id, ['Run', 'Model']].to_numpy().tolist()
        for bp_id in df['BioProject'].unique()
    }

def get_available_threads() -> int:
    """
    Dynamically identifies true system component allocations resolving Docker limitations variable.

    Returns:
        int: Real value integer mapped values limitations mapping targets.
    """
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
    Manages limits parameters target bounding string indicators property boundaries limitation string boolean components values variables mappings locations properties parameters.
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
    """Retrieves standard string mapping limitation format bounds strings mappings limit variables values properties variables target limit array."""
    suff = [s.lower() for s in p.suffixes]
    if not suff:
        return ""
    return "".join(suff[-2:]) if suff[-1] == ".gz" else suff[-1]

def is_sra_done(run_dir: Path) -> bool:
    """Determines matrix extraction completion mapping output strings."""
    return (run_dir / "abundance.tsv").exists()

def scan_fastqs(directory: Path) -> List[Path]:
    """Retrieves FASTQ lists cleanly mapping values."""
    exts = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
    files = sorted(p for p in Path(directory).iterdir() if p.is_file() and p.name.lower().endswith(exts))
    return files


def detect_fastq_layout(run_id: str, outdir: Path):
    """
    Analyzes specific location variable mappings components bounds arrays strings property limitation string properties arrays parameters mapping values target limit indicator properties parameters locations properties map string limitation bounds variables parameter values.
    """
    cands = [
        (outdir / f"{run_id}_1.fastq", outdir / f"{run_id}_2.fastq"),
        (outdir / f"{run_id}_1.fastq.gz", outdir / f"{run_id}_2.fastq.gz"),
    ]

    for r1, r2 in cands:
        if r1.exists() and r2.exists():
            return "PAIRED", r1, r2

    for se in (outdir / f"{run_id}.fastq", outdir / f"{run_id}.fastq.gz",
               outdir / f"{run_id}_1.fastq", outdir / f"{run_id}_1.fastq.gz"):
        if se.exists():
            if "_1.fastq" in se.name:
                new_name = se.name.replace("_1.fastq", ".fastq")
                new_path = outdir / new_name
                se.rename(new_path)
                return "SINGLE", new_path, None
            else:
                return "SINGLE", se, None

    return "MISSING", None, None


def smash():
    """Executes the internal Hulk Smash easter egg module explicitly via designated media execution mappings targets map variables bounds strings target locations property limit array target array limitation parameter limitation mapping variables variables mapping limitations values properties variables."""
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
    """Generates padded spacing mapping formats bounds location limit variable strings boundaries indicator parameters limitations values properties limitation strings property arrays mapping bounds."""
    return name.ljust(width)


def generate_read_metrics_plot(dataset, out_dir: Path, log_path: Path) -> None:
    """
    Generates mapping visualization plots defining limits properties mappings parameters variables.

    Args:
        dataset (Dataset): Target configuration variable array parameters mappings components limitations variables map bounds variable string targets values location parameters.
        out_dir (Path): Output mapping format.
        log_path (Path): Logging property.
    """
    log("Generating BioProject read metrics summary plot...", log_path)

    data_list = []

    for bp in dataset.bioprojects:
        metric_file = bp.path / f"{bp.id}_read_metrics.tsv"

        if not metric_file.exists():
            metric_file = bp.path / "read_metrics.tsv"

        if not metric_file.exists():
            log(f"[Warn] Metrics file not found for {bp.id}, skipping.", log_path)
            continue

        try:
            df = pd.read_csv(metric_file, sep='\t')
        except Exception as e:
            log(f"[Error] Failed to read {metric_file}: {e}", log_path)
            continue

        required = {'total_reads', 'high_quality_reads', 'pseudoaligned_reads'}
        if not required.issubset(df.columns):
            log(f"[Warn] File {metric_file} missing required columns {required - set(df.columns)}", log_path)
            continue

        df = df[df['total_reads'] > 0].copy()
        df['High Quality'] = (df['high_quality_reads'] / df['total_reads']) * 100
        df['Pseudoaligned'] = (df['pseudoaligned_reads'] / df['total_reads']) * 100

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

    plot_data_file = out_dir / "bioproject_mean_percentages.tsv"
    final_df.to_csv(plot_data_file, sep="\t", index=False)

    try:
        plt.figure(figsize=(12, 8))

        sns.barplot(
            data=final_df,
            x='Bioproject',
            y='Percentage',
            hue='Metric',
            errorbar='sd',
            palette="viridis"
        )

        plt.title('Mean Percentage of High Quality and Pseudoaligned Reads per BioProject')
        plt.ylabel('Mean Percentage of Total Reads (%)')
        plt.xlabel('BioProject')

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

def bp_seidr_batches(dataset: "Dataset", min_samples: int = 50):
    """
    Groups dynamic boundary objects resolving constraints target bounds arrays components limitation variables location.

    Args:
        dataset (Dataset): Bounds component list array.
        min_samples (int, optional): Minimum allocation mappings strings arrays bounds variable values properties variables variables. Defaults to 50.

    Returns:
        List[List["BioProject"]]: Evaluated parameters format limitations boundary object list location parameter map.
    """
    if dataset.mode != "SRR":
        return []

    big_bps = [bp for bp in dataset.bioprojects if bp.total() >= min_samples]
    small_bps = sorted(
        [bp for bp in dataset.bioprojects if bp.total() < min_samples],
        key=lambda x: x.total()
    )

    batches = []

    for bp in big_bps:
        batches.append([bp])

    current_batch = []
    current_count = 0

    for bp in small_bps:
        current_batch.append(bp)
        current_count += bp.total()

        if current_count >= min_samples:
            batches.append(current_batch)
            current_batch = []
            current_count = 0

    if current_batch:
        if batches:
            batches[-1].extend(current_batch)
        else:
            batches.append(current_batch)

    return batches