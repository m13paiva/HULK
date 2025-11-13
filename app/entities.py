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
    Central configuration object.
    """

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
        # NEW: feature / tool config that should live *in* Config
        trim_window_size: int = 4,
        trim_mean_quality: int = 20,
        align_method: str = "kallisto",
        kallisto_bootstrap: int = 100,
        tximport_mode: str = "raw_counts",
        tximport_ignore_tx_version: bool = False,
        cache_gb: Optional[int] = None,    # -c/--cache (GiB)
        no_cache: bool = False,            # --no-cache
    ):
        # --------------------------- Runtime (CLI) ---------------------------
        self.input_path = Path(input_path).expanduser().resolve() if input_path is not None else None
        self.reference_path = Path(reference_path).expanduser().resolve() if reference_path is not None else None
        self.outdir = Path(outdir).expanduser().resolve()
        self.min_threads = min_threads
        self.max_threads = max_threads
        self.verbose = verbose
        self.force = force
        self.aggregate = aggregate
        self.dry_run = dry_run
        self.tx2gene = Path(tx2gene).expanduser().resolve() if tx2gene else None
        self.keep_fastq = keep_fastq
        self.error_warnings: List[str] = []

        # --------------------------- Persistent JSON ---------------------------
        # (still available if others want to inspect raw config)
        self.cfg_path = self._resolve_cfg_path()
        self.persisted_cfg = self._load_json_cfg(self.cfg_path)

        # --------------------------- Tool options (effective) ------------------
        # trim
        self.trim_window_size = int(trim_window_size)
        self.trim_mean_quality = int(trim_mean_quality)
        self.trim_opts: Dict[str, Any] = {
            "window_size": self.trim_window_size,
            "mean_quality": self.trim_mean_quality,
        }

        # align
        self.align_method = (align_method or "kallisto").lower()
        self.kallisto_bootstrap = int(kallisto_bootstrap)

        # tximport
        self.tximport_mode = tximport_mode or "raw_counts"
        self.tximport_ignore_tx_version = bool(tximport_ignore_tx_version)
        self.tximport_opts: Dict[str, Any] = {
            "mode": self.tximport_mode,
            "ignore_tx_version": self.tximport_ignore_tx_version,
        }

        # cache
        self.cache_gb: Optional[int] = cache_gb  # may be None → auto
        self.no_cache: bool = bool(no_cache)

        # --------------------------- Directory setup ---------------------------
        self.outdir.mkdir(parents=True, exist_ok=True)

        self.shared = (self.outdir / "shared").resolve()
        if self.shared.exists() and self.force:
            shutil.rmtree(self.shared)
        self.shared.mkdir(parents=True, exist_ok=True)

        self.cache = (self.shared / "cache").resolve()
        if self.cache.exists() and self.force:
            shutil.rmtree(self.cache)
        self.cache.mkdir(parents=True, exist_ok=True)

        self.log = (self.shared / "log.txt").resolve()
        with open(self.log, "w", encoding="utf-8") as f:
            f.write(f"\n===== HULK start {datetime.now().isoformat()} =====\n")

    # ======================================================================
    # Internal helpers
    # ======================================================================

    def _resolve_cfg_path(self) -> Path:
        env_path = os.environ.get(self.CFG_ENV)
        if env_path:
            return Path(env_path).expanduser().resolve()
        return Path.cwd() / self.CFG_DEFAULT_NAME

    def _load_json_cfg(self, path: Path) -> Dict[str, Any]:
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
        return self.max_threads

    def get_tool_opts(self, tool: str) -> Dict[str, Any]:
        """
        Effective tool options, preferring the resolved attributes.
        """
        if tool == "trim":
            return dict(self.trim_opts)
        if tool == "tximport":
            return dict(self.tximport_opts)
        if tool == "align":
            return {
                "method": self.align_method,
                "bootstrap": self.kallisto_bootstrap,
            }
        return self.persisted_cfg.get(tool, {})

    @property
    def cache_high_gb(self) -> int:
        """
        High watermark for cache in GiB.
        If --no-cache, returns 0.
        If -c/--cache was given, use that. Otherwise default = 300.
        """
        if self.no_cache:
            return 0
        if self.cache_gb is not None:
            return int(self.cache_gb)
        return 300

    @property
    def cache_low_gb(self) -> int:
        """
        Low watermark in GiB: 80% of high, or 0 if cache_high_gb <= 0.
        """
        high = self.cache_high_gb
        if high <= 0:
            return 0
        return max(1, int(high * 0.8))

    def summary(self) -> str:
        lines = [
            f"Output dir:   {self.outdir}",
            f"Threads:      min={self.min_threads}, max={self.max_threads}",
            f"Force:        {self.force}",
            f"Verbose:      {self.verbose}",
            f"Aggregate:    {self.aggregate}",
            f"Dry run:      {self.dry_run}",
            f"Keep FASTQ:   {self.keep_fastq}",
            f"Trim:         ws={self.trim_window_size}, mq={self.trim_mean_quality}",
            f"Align:        method={self.align_method}, boot={self.kallisto_bootstrap}",
            f"tximport:     mode={self.tximport_mode}, ignore_ver={self.tximport_ignore_tx_version}",
            f"Cache:        no_cache={self.no_cache}, high={self.cache_high_gb} GiB, low={self.cache_low_gb} GiB",
        ]
        if self.tx2gene:
            lines.append(f"tx2gene:      {self.tx2gene}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"<Config output={self.outdir} threads={self.min_threads}-{self.max_threads}>"

# ------------------------------- Sample -------------------------------

class Sample:
    def __init__(
        self,
        sample_id: str,
        sample_type: str,                   # "SRR" or "FASTQ"
        fastq_paths: Optional[List[Path]] = None,
        metadata: Optional[Dict[str, Any]] = None,
        bioproject: Optional["BioProject"] = None,
        outdir: Optional[Path] = None,
        status: str = "pending",            # "pending" | "ready" | "processing" | "done" | "failed"
    ):
        self.id = str(sample_id)
        self.type = sample_type.upper()
        if self.type not in {"SRR", "FASTQ"}:
            raise ValueError(f"Unsupported sample_type: {sample_type}")
        if self.type == "SRR":
            if bioproject is None:
                raise ValueError("SRR sample requires a BioProject.")
            self.outdir = bioproject.path / self.id
        else:
            if outdir is None:
                raise ValueError("FASTQ sample requires an outdir.")
            self.outdir = outdir / self.id
        self.fastq_paths: List[Path] = list(fastq_paths or [])
        self.metadata: Dict[str, Any] = dict(metadata or {})
        self.status = status
        # dynamically set during prefetch if applicable
        self.sra_path: Optional[Path] = None

    def is_srr(self) -> bool:
        return self.type == "SRR"

    def is_fastq(self) -> bool:
        return self.type == "FASTQ"

    def is_done(self) -> bool:
        ok = (self.outdir / "abundance.tsv").exists()
        if ok:
            self.status = "done"
            return True
        return False

    def __repr__(self) -> str:
        return f"<Sample id={self.id} type={self.type} status={self.status}>"

# ------------------------------ BioProject ----------------------------

class BioProject:
    def __init__(self, bioproject_id: str, base_outdir: Optional[Path] = None):
        self.id = str(bioproject_id)
        self.path = base_outdir / self.id
        self.path.mkdir(parents=True, exist_ok=True)
        self.status: str = "pending"
        self.samples: List["Sample"] = []

    def add_srr(self, srr_id: str, *, metadata: Optional[Dict[str, Any]] = None, status: str = "pending") -> "Sample":
        s = Sample(
            sample_id=srr_id,
            sample_type="SRR",
            fastq_paths=None,
            metadata=metadata,
            bioproject=self,
            outdir=None,
            status=status,
        )
        self.samples.append(s)
        return s

    def total(self) -> int:
        return len(self.samples)

    def done(self) -> int:
        return sum(1 for s in self.samples if s.status == "done")

    def remaining(self) -> int:
        return sum(1 for s in self.samples if s.status not in {"done", "failed"})

    def update_status(self) -> None:
        rem = self.remaining()
        if rem == 0:
            self.status = "done"
        elif self.done() > 0:
            self.status = "active"
        else:
            self.status = "pending"

    def run_postprocessing(self, cfg):
        """
        Run per-BioProject finalization steps:
          1. MultiQC
          2. TPM merge
          3. tximport gene counts (only if user provided -g)
          4. Read metrics table
        """
        from .utils import log, log_err
        from .qc import run_multiqc, build_bp_metrics
        from .tx2gene import merge_bioproject_tpm, bp_gene_counts

        bp_id = self.id
        bp_dir = cfg.outdir / bp_id
        log_path = cfg.log

        errors: list[str] = []

        log(f"[{bp_id}] === Starting BioProject post-processing ===", log_path)

        # ---------------------------------------------------------
        # 1) MULTIQC
        # ---------------------------------------------------------
        try:
            mqc = run_multiqc(bp_dir, log_path, errors)
            if mqc:
                log(f"[{bp_id}] MultiQC complete: {mqc}", log_path)
            else:
                log_err(errors, log_path, f"[{bp_id}] MultiQC returned no output")
        except Exception as e:
            log_err(errors, log_path, f"[{bp_id}] MultiQC failed: {e}")

        # ---------------------------------------------------------
        # 2) TPM MERGE (always run)
        # ---------------------------------------------------------
        try:
            merge_bioproject_tpm(bp_dir, log_path, errors)
            log(f"[{bp_id}] TPM merge complete.", log_path)
        except Exception as e:
            log_err(errors, log_path, f"[{bp_id}] TPM merge failed: {e}")

        # ---------------------------------------------------------
        # 3) TXIMPORT GENE COUNTS  (--gene-counts required)
        # ---------------------------------------------------------
        if cfg.tx2gene:
            try:
                run_ids = [s.id for s in self.samples]
                bp_gene_counts(
                    bp_dir,
                    run_ids,
                    cfg.tx2gene,
                    log_path,
                    errors,
                    tximport_opts=cfg.tximport_opts,
                )
                log(f"[{bp_id}] tximport gene counts complete.", log_path)
            except Exception as e:
                log_err(errors, log_path, f"[{bp_id}] tximport gene counts failed: {e}")
        else:
            log(f"[{bp_id}] No tx2gene provided → skipping gene counts.", log_path)

        # ---------------------------------------------------------
        # 4) READ METRICS TABLE (multiqc-derived)
        # ---------------------------------------------------------
        try:
            build_bp_metrics(bp_id, bp_dir, log_path, errors)
            log(f"[{bp_id}] Read metrics table written.", log_path)
        except Exception as e:
            log_err(errors, log_path, f"[{bp_id}] read metrics failed: {e}")

        log(f"[{bp_id}] === BioProject post-processing complete ===", log_path)

    def __repr__(self):
        return f"<BioProject id={self.id} status={self.status} samples={len(self.samples)}>"

# ------------------------------- Dataset ------------------------------

class Dataset:
    """
    Single-mode collection of samples for a run.
    mode: "SRR" (grouped by BioProject) or "FASTQ" (local fastqs)
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
        self.bioprojects: List[BioProject] = []

    def __len__(self):
        return len(self.samples)

    # SRR helpers
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

    def add_srr(self, bioproject_id: str, srr_id: str, *, metadata: Optional[Dict[str, Any]] = None, status: str = "pending") -> Sample:
        bp = self.get_or_create_bioproject(bioproject_id)
        s = bp.add_srr(srr_id, metadata=metadata, status=status)
        self.samples.append(s)
        return s

    # FASTQ helpers
    def add_fastq(self, sample_id: str, fastq_paths: List[Path], *, outdir: Optional[Path] = None,
                  metadata: Optional[Dict[str, Any]] = None, status: str = "pending") -> Sample:
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
        else:
            for bp in self.bioprojects:
                for s in bp.samples:
                    if s.is_done():
                        s.status = "done"
                bp.update_status()

    def to_do_pairs(self) -> List[Tuple["BioProject", "Sample"]]:
        if self.mode != "SRR":
            return []
        for bp in self.bioprojects:
            bp.update_status()
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
        if self.mode == "SRR":
            return [s for (_bp, s) in self.to_do_pairs()]
        return [s for s in self.samples if s.status not in {"done", "failed"}]

    def bp_done(self):
        return [bp for bp in self.bioprojects if bp.status == "done"]

    def done(self):
        return [s for s in self.samples if s.status == "done"]

    @classmethod
    def from_dataframe(cls, df, cfg: "Config") -> "Dataset":
        ds = cls(dataset_id="run", config=cfg, mode="SRR")
        for _, row in df.iterrows():
            srr = str(row["Run"])
            bp_id = str(row["BioProject"])
            meta = {"Model": row.get("Model", None)} if "Model" in row else {}
            ds.add_srr(bp_id, srr_id=srr, metadata=meta)
        return ds

    @classmethod
    def from_fastq_dir(cls, directory: Path, cfg: "Config") -> "Dataset":
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
