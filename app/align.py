# align.py

import re
from dataclasses import dataclass
from types import MappingProxyType
from pathlib import Path

from .utils import (
    log_err,
    run_cmd
)


# ─────────────────────────────────────────────────────────────────────────────
# Kallisto fragment length presets by sequencer model
# ─────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class FragParams:
    mean: int
    sd:   int

_MODEL_PARAMS = MappingProxyType({
    # NovaSeq
    "ILLUMINA NOVASEQ 6000": FragParams(200, 20),
    "ILLUMINA NOVASEQ X":    FragParams(200, 20),
    "ILLUMINA NOVASEQ X PLUS": FragParams(200, 20),
    # HiSeq family
    "HISEQ X TEN":           FragParams(200, 20),
    "ILLUMINA HISEQ X":      FragParams(200, 20),
    "ILLUMINA HISEQ X TEN":  FragParams(200, 20),
    "ILLUMINA HISEQ 4000":   FragParams(200, 20),
    "ILLUMINA HISEQ 3000":   FragParams(200, 20),
    "ILLUMINA HISEQ 2500":   FragParams(200, 20),
    "ILLUMINA HISEQ 2000":   FragParams(200, 20),
    "ILLUMINA HISEQ 1500":   FragParams(200, 20),
    "ILLUMINA HISEQ 1000":   FragParams(200, 20),
    # NextSeq
    "NEXTSEQ 1000":          FragParams(200, 20),
    "NEXTSEQ 500":           FragParams(200, 20),
    "NEXTSEQ 550":           FragParams(200, 20),
    # BGI/MGI
    "DNBSEQ-G400":           FragParams(170, 15),
    "DNBSEQ-T7":             FragParams(170, 15),
    "BGISEQ-500":            FragParams(170, 15),
    "MGISEQ-2000RS":         FragParams(170, 15),
    # Older Illumina GA
    "ILLUMINA GENOME ANALYZER II":  FragParams(160, 20),
    "ILLUMINA GENOME ANALYZER IIX": FragParams(160, 20),
})

_FALLBACK_RULES = [
    (re.compile(r"\bNOVASEQ\b", re.I),                 FragParams(200, 20)),
    (re.compile(r"\bHISEQ\b", re.I),                   FragParams(200, 20)),
    (re.compile(r"\bNEXTSEQ\b", re.I),                 FragParams(200, 20)),
    (re.compile(r"\b(DNBSEQ|BGISEQ|MGISEQ)\b", re.I),  FragParams(170, 15)),
    (re.compile(r"\b(GENOME ANALYZER|GAII)\b", re.I),  FragParams(160, 20)),
]
_DEFAULT_PARAMS = FragParams(200, 20)

def get_frag_params(platform: str | None) -> FragParams:
    """Return suggested fragment length params for kallisto --single based on sequencer model."""
    if not platform:
        return _DEFAULT_PARAMS
    key = platform.strip().upper()
    if key in _MODEL_PARAMS:
        return _MODEL_PARAMS[key]
    for rx, params in _FALLBACK_RULES:
        if rx.search(platform):
            return params
    return _DEFAULT_PARAMS


# ─────────────────────────────────────────────────────────────────────────────
# Kallisto Index
# ─────────────────────────────────────────────────────────────────────────────

def build_transcriptome_index(transcriptome: Path, shared: Path, log_path: Path) -> Path:
    """Build (or rebuild) a kallisto index at <shared>/transcripts.idx."""
    idx = shared / 'transcripts.idx'
    run_cmd(["kallisto", "index", "-i", str(idx), str(transcriptome)], shared, log_path)
    return idx

# ─────────────────────────────────────────────────────────────────────────────
# Kallisto Quant
# ─────────────────────────────────────────────────────────────────────────────

def kallisto_single_cmd(run_dir: Path, run_id: str, index_path: Path, threads: int,
                        platform: str | None, log_path, error_warnings) -> list[str]:
    fq = pick_fastq(run_dir, run_id, log_path, error_warnings)
    fl = get_frag_params(platform)
    return [
        "kallisto", "quant",
        "--plaintext",                 # ← ensure abundance.tsv is written
        "-i", str(Path(index_path).resolve()),
        "-o", str(run_dir.resolve()),
        "-t", str(threads),
        "--single", "-l", str(fl.mean), "-s", str(fl.sd),
        str(fq.resolve()),
    ]

def kallisto_paired_cmd(run_dir: Path, run_id: str, index_path: Path, threads: int,
                        log_path, error_warnings) -> list[str]:
    r1 = pick_fastq(run_dir, f"{run_id}_1", log_path, error_warnings)
    r2 = pick_fastq(run_dir, f"{run_id}_2", log_path, error_warnings)
    return [
        "kallisto", "quant",
        "--plaintext",                 # ← ensure abundance.tsv is written
        "-i", str(Path(index_path).resolve()),
        "-o", str(run_dir.resolve()),
        "-t", str(threads),
        str(r1.resolve()), str(r2.resolve()),
    ]

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

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