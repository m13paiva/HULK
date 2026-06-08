from __future__ import annotations

import re
from dataclasses import dataclass
from types import MappingProxyType
from pathlib import Path

from .utils import (
    run_cmd,
)


# ─────────────────────────────────────────────────────────────────────────────
# Kallisto fragment length presets by sequencer model
# ─────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class FragParams:
    """
    Data class representing fragment length parameters for single-end alignment.

    Attributes:
        mean (int): The mean fragment length.
        sd (int): The standard deviation of the fragment length.
    """
    mean: int
    sd: int


_MODEL_PARAMS = MappingProxyType({
    # NovaSeq
    "ILLUMINA NOVASEQ 6000": FragParams(200, 20),
    "ILLUMINA NOVASEQ X": FragParams(200, 20),
    "ILLUMINA NOVASEQ X PLUS": FragParams(200, 20),
    # HiSeq family
    "HISEQ X TEN": FragParams(200, 20),
    "ILLUMINA HISEQ X": FragParams(200, 20),
    "ILLUMINA HISEQ X TEN": FragParams(200, 20),
    "ILLUMINA HISEQ 4000": FragParams(200, 20),
    "ILLUMINA HISEQ 3000": FragParams(200, 20),
    "ILLUMINA HISEQ 2500": FragParams(200, 20),
    "ILLUMINA HISEQ 2000": FragParams(200, 20),
    "ILLUMINA HISEQ 1500": FragParams(200, 20),
    "ILLUMINA HISEQ 1000": FragParams(200, 20),
    # NextSeq
    "NEXTSEQ 1000": FragParams(200, 20),
    "NEXTSEQ 500": FragParams(200, 20),
    "NEXTSEQ 550": FragParams(200, 20),
    # BGI/MGI
    "DNBSEQ-G400": FragParams(170, 15),
    "DNBSEQ-T7": FragParams(170, 15),
    "BGISEQ-500": FragParams(170, 15),
    "MGISEQ-2000RS": FragParams(170, 15),
    # Older Illumina GA
    "ILLUMINA GENOME ANALYZER II": FragParams(160, 20),
    "ILLUMINA GENOME ANALYZER IIX": FragParams(160, 20),
})

# Regular expression patterns for fallback fragment length matching based on platform keywords.
_FALLBACK_RULES = [
    (re.compile(r"\bNOVASEQ\b", re.I), FragParams(200, 20)),
    (re.compile(r"\bHISEQ\b", re.I), FragParams(200, 20)),
    (re.compile(r"\bNEXTSEQ\b", re.I), FragParams(200, 20)),
    (re.compile(r"\b(DNBSEQ|BGISEQ|MGISEQ)\b", re.I), FragParams(170, 15)),
    (re.compile(r"\b(GENOME ANALYZER|GAII)\b", re.I), FragParams(160, 20)),
]
_DEFAULT_PARAMS = FragParams(200, 20)


def get_frag_params(platform: str | None) -> FragParams:
    """
    Determines suggested fragment length parameters for Kallisto single-end quantification
    based on the provided sequencer model or platform string.

    Args:
        platform (str | None): The sequencing platform string (e.g., 'ILLUMINA NOVASEQ 6000').

    Returns:
        FragParams: The inferred mean and standard deviation for the fragment length.
    """
    if not platform:
        return _DEFAULT_PARAMS

    key = platform.strip().upper()

    # Exact match check
    if key in _MODEL_PARAMS:
        return _MODEL_PARAMS[key]

    # Fallback heuristic checks
    for rx, params in _FALLBACK_RULES:
        if rx.search(platform):
            return params

    return _DEFAULT_PARAMS


def list_known_seq_techs() -> list[str]:
    """
    Retrieves a sorted list of all exact sequencing technologies recognized by the pipeline.

    Returns:
        list[str]: Alphabetically sorted list of known sequencer platform identifiers.
    """
    return sorted(_MODEL_PARAMS.keys())


# ─────────────────────────────────────────────────────────────────────────────
# Kallisto Index
# ─────────────────────────────────────────────────────────────────────────────

def build_transcriptome_index(transcriptome: Path, shared: Path, log_path: Path) -> Path:
    """
    Builds or rebuilds a Kallisto index from a provided transcriptome FASTA file.

    This function should typically be called once at the top level of the pipeline,
    with the resulting index path passed to worker processes.

    Args:
        transcriptome (Path): Path to the reference transcriptome FASTA file.
        shared (Path): Directory where the generated index should be stored.
        log_path (Path): Path to the execution log file.

    Returns:
        Path: The absolute path to the generated Kallisto index file (`transcripts.idx`).
    """
    idx = shared / "transcripts.idx"
    run_cmd(["kallisto", "index", "-i", str(idx), str(transcriptome)], shared, log_path)
    return idx


# ─────────────────────────────────────────────────────────────────────────────
# Kallisto Quant
# ─────────────────────────────────────────────────────────────────────────────

def kallisto_single_cmd(
        run_dir: Path,
        run_id: str,
        index_path: Path,
        threads: int,
        platform: str | None,
        log_path: Path,
        error_warnings: list[str],
        *,
        bootstraps: int = 100,
) -> list[str]:
    """
    Constructs the `kallisto quant` command array for single-end FASTQ data.

    Args:
        run_dir (Path): The working directory containing the FASTQ files.
        run_id (str): The unique identifier for the run/sample.
        index_path (Path): Path to the Kallisto index.
        threads (int): Number of threads to allocate for the process.
        platform (str | None): Sequencing platform, used to infer fragment length parameters.
        log_path (Path): Path to the log file.
        error_warnings (list[str]): List reference to append warnings during file detection.
        bootstraps (int, optional): Number of bootstrap samples. Defaults to 100.

    Returns:
        list[str]: The complete Kallisto command array ready for subprocess execution.
    """
    # Attempt to locate trimmed FASTQs, falling back to raw if necessary.
    fq = pick_fastq(run_dir, run_id, log_path, error_warnings)
    fl = get_frag_params(platform)

    return [
        "kallisto", "quant",
        "--plaintext",
        "-i", str(Path(index_path).resolve()),
        "-o", str(run_dir.resolve()),
        "-t", str(threads),
        "-b", str(int(bootstraps)),
        "--single", "-l", str(fl.mean), "-s", str(fl.sd),
        str(fq.resolve()),
    ]


def kallisto_paired_cmd(
        run_dir: Path,
        run_id: str,
        index_path: Path,
        threads: int,
        log_path: Path,
        error_warnings: list[str],
        *,
        bootstraps: int = 100,
) -> list[str]:
    """
    Constructs the `kallisto quant` command array for paired-end FASTQ data.

    Args:
        run_dir (Path): The working directory containing the FASTQ files.
        run_id (str): The unique identifier for the run/sample.
        index_path (Path): Path to the Kallisto index.
        threads (int): Number of threads to allocate for the process.
        log_path (Path): Path to the log file.
        error_warnings (list[str]): List reference to append warnings during file detection.
        bootstraps (int, optional): Number of bootstrap samples. Defaults to 100.

    Returns:
        list[str]: The complete Kallisto command array ready for subprocess execution.
    """
    # Locate read 1 and read 2 FASTQ files.
    r1 = pick_fastq(run_dir, f"{run_id}_1", log_path, error_warnings)
    r2 = pick_fastq(run_dir, f"{run_id}_2", log_path, error_warnings)

    return [
        "kallisto", "quant",
        "--plaintext",
        "-i", str(Path(index_path).resolve()),
        "-o", str(run_dir.resolve()),
        "-t", str(threads),
        "-b", str(int(bootstraps)),
        str(r1.resolve()), str(r2.resolve()),
    ]


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def pick_fastq(run_dir: Path, stem: str, log_path: Path, error_warnings: list[str]) -> Path:
    """
    Selects the most appropriate FASTQ file for a given sample stem.
    It prioritizes trimmed files but will fall back to raw files if trimmed variants are unavailable.

    Args:
        run_dir (Path): The directory containing the target FASTQ files.
        stem (str): The base name stem of the FASTQ file (e.g., sample ID or sample_1).
        log_path (Path): Path to the log file.
        error_warnings (list[str]): List reference to append warnings.

    Raises:
        FileNotFoundError: If neither trimmed nor raw FASTQ files can be found matching the stem.

    Returns:
        Path: The path to the resolved FASTQ file.
    """
    # 1. Search for trimmed variants first.
    for suff in (".trim.fastq.gz", ".trim.fastq"):
        p = run_dir / f"{stem}{suff}"
        if p.exists():
            return p

    # 2. Search for untrimmed variants as a fallback mechanism.
    for suff in (".fastq.gz", ".fastq"):
        p = run_dir / f"{stem}{suff}"
        if p.exists():
            print(f"[{stem}] No trimmed FASTQs found in {run_dir}; quantifying on untrimmed reads.")
            return p

    raise FileNotFoundError(f"No FASTQ found for stem '{stem}' in {run_dir}")