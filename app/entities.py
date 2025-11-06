from __future__ import annotations
import os
import json
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, Optional, List

# ------------------------------- Config -------------------------------
class Config:
    """
    Central configuration object

    Combines:
      • Runtime options passed from CLI (input/output paths, threads, etc.)
      • Persistent defaults loaded from a JSON config (~/.hulk.json or HULK_CONFIG env)
      • Tool-specific options (fastp, kallisto, tximport, etc.)

    This class provides unified access to all configuration settings
    and automatically creates directories when needed.
    """

    # Defaults from CLI-level constants
    CFG_ENV = "HULK_CONFIG"
    CFG_DEFAULT_NAME = ".hulk.json"

    def __init__(
        self,
        *,
        input_path: Optional[Path] = None,
        reference_path: Optional[Path] = None,
        output_dir: Path,
        min_threads: int = 4,
        max_threads: int = 10,
        verbose: bool = False,
        force: bool = False,
        aggregate: bool = False,
        dry_run: bool = False,
        tx2gene_path: Optional[Path] = None,
        keep_fastq: bool = False,
    ):
        # --------------------------- Runtime (CLI) ---------------------------
        self.input_path = input_path
        self.reference_path = reference_path
        self.output_dir = Path(output_dir).expanduser().resolve()
        self.min_threads = min_threads
        self.max_threads = max_threads
        self.verbose = verbose
        self.force = force
        self.aggregate = aggregate
        self.dry_run = dry_run
        self.tx2gene_path = Path(tx2gene_path).expanduser().resolve() if tx2gene_path else None
        self.keep_fastq = keep_fastq

        # --------------------------- Persistent JSON ---------------------------
        self.cfg_path = self._resolve_cfg_path()
        self.persisted_cfg = self._load_json_cfg(self.cfg_path)

        # --------------------------- Tool options ---------------------------
        self.trim_opts = self.persisted_cfg.get("trim", {})
        self.tximport_opts = self.persisted_cfg.get("tximport", {})

        # --------------------------- Directory setup ---------------------------
        # <output_dir>/
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # <output_dir>/shared/
        self.shared_dir = self.output_dir / "shared"
        self.shared_dir.mkdir(parents=True, exist_ok=True)

        # <output_dir>/shared/log.txt
        self.log_path = self.shared_dir / "log.txt"
        if not self.log_path.exists():
            with open(self.log_path, "w", encoding="utf-8") as f:
                f.write(f"\n===== HULK start {datetime.now().isoformat()} =====\n")

    # ======================================================================
    # Internal helpers
    # ======================================================================

    def _resolve_cfg_path(self) -> Path:
        """Resolve which configuration file to use."""
        env_path = os.environ.get(self.CFG_ENV)
        if env_path:
            return Path(env_path).expanduser().resolve()
        return Path.cwd() / self.CFG_DEFAULT_NAME

    def _load_json_cfg(self, path: Path) -> Dict[str, Any]:
        """Load persistent JSON config (safe)."""
        if not path.exists():
            return {}
        try:
            with open(path, "r", encoding="utf-8") as fh:
                return json.load(fh)
        except Exception:
            return {}

    # ======================================================================
    # Convenience accessors
    # ======================================================================

    @property
    def threads(self) -> int:
        """Return max threads (alias for convenience)."""
        return self.max_threads

    def get_tool_opts(self, tool: str) -> Dict[str, Any]:
        """Generic accessor for tool-specific options."""
        return self.persisted_cfg.get(tool, {})

    def summary(self) -> str:
        """Return a formatted string summarizing this configuration."""
        lines = [
            f"Output dir:   {self.output_dir}",
            f"Threads:      min={self.min_threads}, max={self.max_threads}",
            f"Force:        {self.force}",
            f"Verbose:      {self.verbose}",
            f"Aggregate:    {self.aggregate}",
            f"Dry run:      {self.dry_run}",
            f"Keep FASTQ:   {self.keep_fastq}",
        ]
        if self.tx2gene_path:
            lines.append(f"tx2gene:      {self.tx2gene_path}")
        if self.trim_opts:
            lines.append(f"Trim opts:    {self.trim_opts}")
        if self.tximport_opts:
            lines.append(f"Tximport:     {self.tximport_opts}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"<Config output={self.output_dir} threads={self.min_threads}-{self.max_threads}>"

# ------------------------------- Sample -------------------------------

class Sample:
    """
    Represents a single RNA-seq sample (generic):
    - type: "SRR" (SRA accession) or "FASTQ" (local fastqs)
    - outdir is namespaced per-sample under its parent dataset's outdir
    """

    def __init__(
        self,
        sample_id: str,
        sample_type: str,                   # "SRR" or "FASTQ"
        dataset_outdir: Path,
        config: Config,
        fastq_paths: Optional[List[Path]] = None,
        metadata: Optional[Dict[str, Any]] = None,
        status: str = "pending",            # "pending" | "ready" | "done" | "failed"
    ):
        self.id = str(sample_id)
        self.type = sample_type.upper()
        if self.type not in {"SRR", "FASTQ"}:
            raise ValueError(f"Unsupported sample_type: {sample_type}")
        self.outdir = Path(dataset_outdir) / self.id
        self.config = config
        self.fastq_paths: List[Path] = list(fastq_paths or [])
        self.metadata: Dict[str, Any] = dict(metadata or {})
        self.status = status

        self.outdir.mkdir(parents=True, exist_ok=True)

    def is_srr(self) -> bool:
        return self.type == "SRR"

    def is_fastq(self) -> bool:
        return self.type == "FASTQ"

    def __repr__(self) -> str:
        return f"<Sample id={self.id} type={self.type} status={self.status}>"


# ------------------------------- Dataset ------------------------------

class Dataset:
    """
    General collection of samples (SRR and/or FASTQ) under one logical run.
    This is the primary container the pipeline will operate on.
    """

    def __init__(self, dataset_id: str, config: Config):
        self.id = str(dataset_id)
        self.config = config
        self.outdir = self.config.outdir / self.id
        self.samples: List[Sample] = []

        self.outdir.mkdir(parents=True, exist_ok=True)

    # ---- sample management

    def add_sample(
        self,
        sample_id: str,
        sample_type: str,                    # "SRR" or "FASTQ"
        fastq_paths: Optional[List[Path]] = None,
        metadata: Optional[Dict[str, Any]] = None,
        status: str = "pending",
    ) -> Sample:
        s = Sample(
            sample_id=sample_id,
            sample_type=sample_type,
            dataset_outdir=self.outdir,
            config=self.config,
            fastq_paths=fastq_paths,
            metadata=metadata,
            status=status,
        )
        self.samples.append(s)
        return s

    def add_srr(self, srr_id: str, metadata: Optional[Dict[str, Any]] = None) -> Sample:
        return self.add_sample(sample_id=srr_id, sample_type="SRR", metadata=metadata)

    def add_fastq(
        self,
        sample_id: str,
        fastq_paths: List[Path],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> Sample:
        return self.add_sample(
            sample_id=sample_id, sample_type="FASTQ", fastq_paths=fastq_paths, metadata=metadata
        )

    # ---- queries

    def get_sample(self, sample_id: str) -> Optional[Sample]:
        for s in self.samples:
            if s.id == sample_id:
                return s
        return None

    def iter_srr(self):
        for s in self.samples:
            if s.is_srr():
                yield s

    def iter_fastq(self):
        for s in self.samples:
            if s.is_fastq():
                yield s

    def __repr__(self) -> str:
        return f"<Dataset id={self.id} n_samples={len(self.samples)}>"


# ------------------------------ BioProject ----------------------------

class BioProject(Dataset):
    """
    Specialized Dataset for SRA-only collections grouped by a BioProject (e.g., PRJNAxxxxxx).
    Enforces that all added samples are SRR accessions and provides SRA-specific helpers.
    """

    def __init__(self, bioproject_id: str, config: Config):
        super().__init__(dataset_id=bioproject_id, config=config)
        self.bioproject_id = bioproject_id  # alias for clarity

    # Enforce SRR-only semantics
    def add_sample(
        self,
        sample_id: str,
        sample_type: str,
        fastq_paths: Optional[List[Path]] = None,
        metadata: Optional[Dict[str, Any]] = None,
        status: str = "pending",
    ) -> Sample:
        if sample_type.upper() != "SRR":
            raise ValueError("BioProject only accepts SRR samples.")
        return super().add_sample(sample_id, "SRR", fastq_paths=None, metadata=metadata, status=status)

    def add_srr(self, srr_id: str, metadata: Optional[Dict[str, Any]] = None) -> Sample:
        return super().add_srr(srr_id, metadata=metadata)

    # Optional helpers you can flesh out later:
    # - ingest_from_table(tsv_path, srr_column="Run", extra_cols=[...])
    # - resolve_runinfo()
    # - prefetch_all(), fasterq_all(), etc., if you prefer orchestration here.

    def __repr__(self) -> str:
        return f"<BioProject id={self.bioproject_id} n_samples={len(self.samples)}>"
