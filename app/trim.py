# trim.py

from pathlib import Path

# ─────────────────────────────────────────────────────────────────────────────
# Fastp commands
# ─────────────────────────────────────────────────────────────────────────────

def fastp_single_cmd(run_dir: Path, run_id: str, r: Path, threads: int,window_size=4,mean_quality=20) -> list[str]:
    """fastp command (SINGLE-END)."""
    return [
        "fastp",
        "-i", r.name,
        "-o", f"{run_id}.trim.fastq.gz",
        "-w", str(threads),
        "--cut_front", "--cut_front_window_size", str(window_size), "--cut_front_mean_quality", str(mean_quality),
        "--cut_tail",  "--cut_tail_window_size",  str(window_size), "--cut_tail_mean_quality",  str(mean_quality),
        "-l", "50",
        "-j", f"{run_id}.fastp.json",
        "-h", f"{run_id}.fastp.html",
    ]

def fastp_paired_cmd(run_dir: Path, run_id: str, r1: Path, r2: Path, threads: int,window_size=4,mean_quality=20) -> list[str]:
    """fastp command (PAIRED-END)."""
    return [
        "fastp",
        "-i", r1.name, "-I", r2.name,
        "-o", f"{run_id}_1.trim.fastq.gz", "-O", f"{run_id}_2.trim.fastq.gz",
        "--detect_adapter_for_pe",
        "-w", str(threads),
        "--cut_front", "--cut_front_window_size", str(window_size), "--cut_front_mean_quality", str(mean_quality),
        "--cut_tail",  "--cut_tail_window_size",  str(window_size), "--cut_tail_mean_quality",  str(mean_quality),
        "-l", "50",
        "-j", f"{run_id}.fastp.json",
        "-h", f"{run_id}.fastp.html",
    ]