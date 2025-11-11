from __future__ import annotations
import os
import json
import shutil
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, Optional, List, Tuple

FASTQ_EXTS = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}

def _is_fastq(p: Path) -> bool:
    name = p.name.lower()
    return any(name.endswith(ext) for ext in FASTQ_EXTS)

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
        outdir: Path,
        min_threads: int = 4,
        max_threads: int = 10,
        verbose: bool = False,
        force: bool = False,
        aggregate: bool = False,
        dry_run: bool = False,
        tx2gene: Optional[Path] = None,
        keep_fastq: bool = False,
    ):
        # --------------------------- Runtime (CLI) ---------------------------
        self.input_path = input_path
        self.reference_path = reference_path
        self.outdir = Path(outdir).expanduser().resolve()
        self.min_threads = min_threads
        self.max_threads = max_threads
        self.verbose = verbose
        self.force = force
        self.aggregate = aggregate
        self.dry_run = dry_run
        self.tx2gene = Path(tx2gene).expanduser().resolve() if tx2gene else None
        self.keep_fastq = keep_fastq
        self.error_warnings = []

        # --------------------------- Persistent JSON ---------------------------
        self.cfg_path = self._resolve_cfg_path()
        self.persisted_cfg = self._load_json_cfg(self.cfg_path)

        # --------------------------- Tool options ---------------------------
        self.trim_opts = self.persisted_cfg.get("trim", {})
        self.tximport_opts = self.persisted_cfg.get("tximport", {})

        # --------------------------- Directory setup ---------------------------
        # <outdir>/
        self.outdir = Path(self.outdir).expanduser().resolve()
        self.outdir.mkdir(parents=True, exist_ok=True)

        # <outdir>/shared/
        self.shared = (self.outdir / "shared").resolve()
        if self.shared.exists():
            shutil.rmtree(self.shared)
        self.shared.mkdir(parents=True, exist_ok=True)

        # <outdir>/shared/cache/
        self.cache = (self.shared / "cache").resolve()
        if self.cache.exists():
            shutil.rmtree(self.cache)
        self.cache.mkdir(parents=True, exist_ok=True)

        # <outdir>/shared/log.txt
        self.log = (self.shared / "log.txt").resolve()

        with open(self.log, "w", encoding="utf-8") as f:
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
            f"Output dir:   {self.outdir}",
            f"Threads:      min={self.min_threads}, max={self.max_threads}",
            f"Force:        {self.force}",
            f"Verbose:      {self.verbose}",
            f"Aggregate:    {self.aggregate}",
            f"Dry run:      {self.dry_run}",
            f"Keep FASTQ:   {self.keep_fastq}",
        ]
        if self.tx2gene:
            lines.append(f"tx2gene:      {self.tx2gene}")
        if self.trim_opts:
            lines.append(f"Trim opts:    {self.trim_opts}")
        if self.tximport_opts:
            lines.append(f"Tximport:     {self.tximport_opts}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"<Config output={self.outdir} threads={self.min_threads}-{self.max_threads}>"

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
        fastq_paths: Optional[List[Path]] = None,
        metadata: Optional[Dict[str, Any]] = None,
        bioproject: Optional[BioProject] = None,
        outdir: Optional[Path] = None,
        status: str = "pending",            # "pending" | "ready" | "done" | "failed"
    ):
        self.id = str(sample_id)
        self.type = sample_type.upper()
        if self.type not in {"SRR", "FASTQ"}:
            raise ValueError(f"Unsupported sample_type: {sample_type}")
        if self.type == "SRR":
            self.outdir = bioproject.path / self.id
        else:
            self.outdir = outdir / self.id
        self.fastq_paths: List[Path] = list(fastq_paths or [])
        self.metadata: Dict[str, Any] = dict(metadata or {})
        self.status = status

    def is_srr(self) -> bool:
        return self.type == "SRR"

    def is_fastq(self) -> bool:
        return self.type == "FASTQ"

    def is_done(self) -> bool:
        """
        Return True if this sample has completed processing.
        For SRR samples, that means 'abundance.tsv' exists in the sample's outdir.
        For FASTQ samples, uses the same condition (abundance.tsv).
        """
        if (self.outdir / "abundance.tsv").exists():
            self.status="done"
            return True

    def __repr__(self) -> str:
        return f"<Sample id={self.id} type={self.type} status={self.status}>"

# ------------------------------ BioProject ----------------------------

class BioProject:
    """
    SRR-only container under <output>/<BioProject>.
    Provides a .path used by Sample(type="SRR") to place per-SRR outdirs.
    """

    def __init__(self, bioproject_id: str, base_outdir: Optional[Path] = None):
        self.id = str(bioproject_id)
        self.path = base_outdir / self.id
        self.path.mkdir(parents=True, exist_ok=True)
        self.status: str = "pending"
        self.samples: List["Sample"] = []

    def add_srr(
        self,
        srr_id: str,
        *,
        metadata: Optional[Dict[str, Any]] = None,
        status: str = "pending",
    ) -> "Sample":
        """Create an SRR Sample under this BioProject."""
        s = Sample(
            sample_id=srr_id,
            sample_type="SRR",
            fastq_paths=None,            # will be produced later by fasterq-dump
            metadata=metadata,
            bioproject=self,             # drives <BioProject>/<SRR> placement
            outdir=None,                 # ignored for SRR (bioproject.path is used)
            status=status,
        )
        self.samples.append(s)
        return s

    def total(self) -> int:
        """Total number of samples in this bioproject."""
        return len(self.samples)

    def done(self) -> int:
        """Number of completed samples."""
        return sum(1 for s in self.samples if s.status == "done")

    def remaining(self) -> int:
        """Number of unfinished samples."""
        return sum(1 for s in self.samples if s.status not in {"done", "failed"})

    def update_status(self) -> None:
        """Refresh project status based on contained samples."""
        rem = self.remaining()
        if rem == 0:
            self.status = "done"
        elif self.done() > 0:
            self.status = "active"
        else:
            self.status = "pending"

    def __repr__(self):
        return f"<BioProject id={self.id} status={self.status} samples={len(self.samples)}>"


# ------------------------------- Dataset ------------------------------


class Dataset:
    """
    Single-mode collection of samples for a run.

    mode:
      - "SRR": expects SRR accessions grouped by BioProject
      - "FASTQ": expects local FASTQ samples (each file becomes a sample)
    """

    def __init__(self, dataset_id: str, config: "Config", mode: str):
        self.id = str(dataset_id)
        self.config = config
        self.mode = mode.upper()
        if self.mode not in {"SRR", "FASTQ"}:
            raise ValueError("Dataset mode must be 'SRR' or 'FASTQ'.")

        self.path = config.outdir
        self.path.mkdir(parents=True, exist_ok=True)

        self.fastq_root = self.path / "fastq_samples"
        if self.mode == "FASTQ":
            self.fastq_root.mkdir(parents=True, exist_ok=True)

        self.samples: List[Sample] = []
        self.bioprojects: List[BioProject] = []   # <-- list instead of dict

    def __len__(self):
        return len(self.samples)

    # ------------------- SRR helpers -------------------
    def _find_bioproject(self, bioproject_id: str) -> Optional[BioProject]:
        for bp in self.bioprojects:
            if getattr(bp, "id", None) == bioproject_id:
                return bp
        return None

    def get_or_create_bioproject(self, bioproject_id: str) -> BioProject:
        bp = self._find_bioproject(bioproject_id)
        if bp is None:
            bp = BioProject(bioproject_id, base_outdir=self.path)
            self.bioprojects.append(bp)
        return bp

    def add_srr(
        self,
        bioproject_id: str,
        srr_id: str,
        *,
        metadata: Optional[Dict[str, Any]] = None,
        status: str = "pending",
    ) -> Sample:

        bp = self.get_or_create_bioproject(bioproject_id)
        s = bp.add_srr(srr_id, metadata=metadata, status=status)
        self.samples.append(s)
        return s

    # ------------------- FASTQ helpers -------------------
    def add_fastq(
        self,
        sample_id: str,
        fastq_paths: List[Path],
        *,
        outdir: Optional[Path] = None,
        metadata: Optional[Dict[str, Any]] = None,
        status: str = "pending",
    ) -> Sample:

        target_outdir = Path(outdir) if outdir else self.fastq_root
        target_outdir.mkdir(parents=True, exist_ok=True)
        s = Sample(
            sample_id=sample_id,
            sample_type="FASTQ",
            fastq_paths=fastq_paths,
            metadata=metadata,
            bioproject=None,
            outdir=target_outdir,
            status=status,
        )
        self.samples.append(s)
        return s

    def update_status(self):
        if self.mode == "FASTQ":
            for s in self.samples:
                if s.is_done():
                    s.status = "done"
        else:  # SRR
            for bp in self.bioprojects:
                for s in bp.samples:
                    if s.is_done():
                        s.status = "done"
                bp.update_status()

    def to_do_pairs(self) -> List[Tuple["BioProject", "Sample"]]:
        """
        SRR mode: priority-ordered (bp, sample) needing work.
        FASTQ mode: not used (call .to_do() instead).
        """
        if self.mode != "SRR":
            return []
        # keep BP statuses fresh
        for bp in self.bioprojects:
            bp.update_status()
        # fewest remaining first (closes BPs sooner)
        ordered_bps = sorted(
            (bp for bp in self.bioprojects if bp.remaining() > 0),
            key=lambda b: (b.remaining(), b.total(), b.id)
        )
        pairs: List[Tuple["BioProject", "Sample"]] = []
        for bp in ordered_bps:
            for s in bp.samples:
                if s.status not in {"done", "failed"}:
                    pairs.append((bp, s))
        return pairs

    def to_do(self) -> List["Sample"]:
        """
        Unified pending samples list (no tuples).
        - SRR: returns samples in BP-prioritized order
        - FASTQ: returns pending FASTQ samples
        """
        if self.mode == "SRR":
            return [s for (_bp, s) in self.to_do_pairs()]
        # FASTQ
        return [s for s in self.samples if s.status not in {"done", "failed"}]

    def bp_done(self):
        return [bp for bp in self.bioprojects if bp.status == "done"]

    def done(self):
        return [s for s in self.samples if s.status == "done"]

    def to_do(self):
        if self.mode=="SRR":
            return [
                sample
                for bp in self.bioprojects
                if bp.status != "done"
                for sample in bp.samples
                if sample.status != "done"
                ]
        if self.mode=="FASTQ":
            return [
                sample
                for sample in self.samples
                if sample.status != "done"
                ]


    # ------------------- Constructors -------------------
    @classmethod
    def from_dataframe(cls, df, cfg: "Config") -> "Dataset":
        """
        Build an SRR-mode Dataset from a table with columns: Run, BioProject, Model.
        """
        ds = cls(dataset_id="run", config=cfg, mode="SRR")
        for _, row in df.iterrows():
            srr = str(row["Run"])
            bp_id = str(row["BioProject"])
            meta = {"Model": row.get("Model", None)} if "Model" in row else {}
            ds.add_srr(bp_id, srr_id=srr, metadata=meta)
        return ds

    @classmethod
    def from_fastq_dir(cls, directory: Path, cfg: "Config") -> "Dataset":
        """
        Build a FASTQ-mode Dataset from a directory of FASTQ files.
        Each FASTQ file is treated as a separate sample.
        Sample ID = file name (must be unique in the directory).
        """
        directory = Path(directory).expanduser().resolve()
        if not directory.is_dir():
            raise ValueError(f"FASTQ input must be a directory: {directory}")

        files = sorted(p for p in directory.iterdir() if p.is_file() and _is_fastq(p))
        if not files:
            raise ValueError(f"No FASTQ files found in: {directory}")

        names = [p.name for p in files]
        if len(set(names)) != len(names):
            raise ValueError("Duplicate FASTQ filenames detected; filenames must be unique for sample IDs.")

        ds = cls(dataset_id="run", config=cfg, mode="FASTQ")
        for f in files:
            ds.add_fastq(sample_id=f.name, fastq_paths=[f])
        return ds



    def __repr__(self) -> str:
        return f"<Dataset id={self.id} mode={self.mode} n_samples={len(self.samples)} n_bioprojects={len(self.bioprojects)}>"
