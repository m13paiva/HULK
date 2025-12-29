from __future__ import annotations
import os
import json
import shutil
from pathlib import Path
from datetime import datetime
from typing import Any, Dict, Optional, List, Tuple

_FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".trim.fastq")


def _is_fastq_path(p: Path) -> bool:
    """
    Return True if p looks like a FASTQ/FASTQ.GZ file.
    """
    if not p.is_file():
        return False
    name = p.name.lower()
    return any(name.endswith(suf) for suf in _FASTQ_SUFFIXES)


# ------------------------------- Config -------------------------------

import os
import json
from pathlib import Path
from typing import Optional, List, Dict, Any


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
            dry_run: bool = False,
            tx2gene: Optional[Path] = None,
            keep_fastq: bool = False,
            seq_tech: Optional[str] = None,

            # Tool config
            trim_window_size: int = 4,
            trim_mean_quality: int = 20,
            align_method: str = "kallisto",
            kallisto_bootstrap: int = 100,
            tximport_mode: str = "raw_counts",
            tximport_ignore_tx_version: bool = False,
            cache_gb: Optional[int] = None,
            no_cache: bool = False,

            # Analysis flow
            rem_missing_bps: bool = False,
            target_genes_files: Optional[List[Path]] = None,

            # Plotting options
            plot_pca: bool = True,
            plot_heatmap: bool = True,
            plot_var_heatmap: bool = True,
            plot_sample_cor: bool = True,
            plot_dispersion: bool = True,
            top_n_vars: int = 500,

            plots_only_mode: bool = False,
            tximport_only_mode: bool = False,
            bioproject_filter: Optional[str] = None,

            # Post-processing flags
            no_bp_postprocessing: bool = False,
            no_global_postprocessing: bool = False,
    ):
        # --------------------------- Runtime (CLI) ---------------------------
        self.input_path = Path(input_path).expanduser().resolve() if input_path is not None else None
        self.reference_path = Path(reference_path).expanduser().resolve() if reference_path is not None else None
        self.outdir = Path(outdir).expanduser().resolve()
        self.min_threads = min_threads
        self.max_threads = max_threads
        self.verbose = verbose
        self.force = force
        self.dry_run = dry_run
        self.tx2gene = Path(tx2gene).expanduser().resolve() if tx2gene else None
        self.keep_fastq = keep_fastq
        self.error_warnings: List[str] = []
        self.rem_missing_bps = rem_missing_bps
        self.seq_tech: Optional[str] = seq_tech

        # --------------------------- Persistent JSON ---------------------------
        self.cfg_path = self._resolve_cfg_path()
        self.persisted_cfg = self._load_json_cfg(self.cfg_path)

        # --------------------------- Tool options ---------------------------
        # Trim
        self.trim_window_size = int(trim_window_size)
        self.trim_mean_quality = int(trim_mean_quality)
        self.trim_opts: Dict[str, Any] = {
            "window_size": self.trim_window_size,
            "mean_quality": self.trim_mean_quality,
        }

        # Align
        self.align_method = (align_method or "kallisto").lower()
        self.kallisto_bootstrap = int(kallisto_bootstrap)

        # Tximport
        self.tximport_mode = tximport_mode or "raw_counts"
        self.tximport_ignore_tx_version = bool(tximport_ignore_tx_version)
        self.tximport_opts: Dict[str, Any] = {
            "mode": self.tximport_mode,
            "ignore_tx_version": self.tximport_ignore_tx_version,
        }

        # --------------------------- DESeq2 / Expression ---------------------
        deseq2_cfg = self.persisted_cfg.get("deseq2", {})
        self.deseq2_vst_enabled: bool = bool(deseq2_cfg.get("enabled", True))
        self.deseq2_var_threshold: float = float(deseq2_cfg.get("var_threshold", 0.1))
        self.expr_use_matrix: str = "vst"
        self.drop_nonvarying_genes: bool = True

        # Plot Flags
        self.plot_pca: bool = bool(plot_pca)
        self.plot_heatmap: bool = bool(plot_heatmap)
        self.plot_var_heatmap: bool = bool(plot_var_heatmap)
        self.plot_sample_cor: bool = bool(plot_sample_cor)
        self.plot_dispersion: bool = bool(plot_dispersion)
        self.top_n_vars: int = int(top_n_vars)

        self.plots_only_mode: bool = bool(plots_only_mode)
        self.tximport_only_mode: bool = bool(tximport_only_mode)

        # Post-processing control
        self.no_bp_postprocessing = no_bp_postprocessing
        self.no_global_postprocessing = no_global_postprocessing

        # Handle multiple target files (CLI overrides Persisted)
        self.target_genes_files: List[Path] = []
        if target_genes_files:
            for tf in target_genes_files:
                self.target_genes_files.append(Path(tf).expanduser().resolve())

        self.bioproject_filter = bioproject_filter

        # --------------------------- Seidr / Network ---------------------------
        seidr_cfg = self.persisted_cfg.get("seidr", {})

        self.seidr_enabled = bool(seidr_cfg.get("enabled", True))
        self.seidr_preset = seidr_cfg.get("preset", "BALANCED").upper()
        self.seidr_algos = seidr_cfg.get("algorithms", [])
        self.seidr_backbone = float(seidr_cfg.get("backbone", 1.28))
        self.seidr_workers = int(seidr_cfg.get("workers", 2))
        self.seidr_target_mode = seidr_cfg.get("target_mode", "targeted_only")

        # <--- LOOK HERE: Load persisted force preference
        self.seidr_force = bool(seidr_cfg.get("force", False))

        # Persisted defaults for targets (if any)
        self.seidr_persisted_targets = [Path(t) for t in seidr_cfg.get("targets", [])]

        # Effective Seidr Targets: CLI targets take precedence if present
        self.seidr_targets = self.target_genes_files if self.target_genes_files else self.seidr_persisted_targets

        # Cache
        self.cache_gb: Optional[int] = cache_gb
        self.no_cache: bool = bool(no_cache)

        # --------------------------- Directory setup ---------------------------
        self.shared = (self.outdir / "shared").resolve()
        self.cache = (self.shared / "cache").resolve()
        self.log = (self.shared / "log.txt").resolve()

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
                "seq_tech": self.seq_tech,
            }
        if tool == "seidr":
            # Paths are standard in the pipeline flow:
            return {
                "enabled": self.seidr_enabled,
                "preset": self.seidr_preset,
                "algorithms": self.seidr_algos,
                "backbone": self.seidr_backbone,
                "workers": self.seidr_workers,
                "target_mode": self.seidr_target_mode,

                # <--- LOOK HERE: This line is required to pass the setting to Seidr!
                "force": self.seidr_force,

                "targets": [str(t) for t in self.seidr_targets],
                # Standard input locations generated by post_processing.R
                "genes_file": self.shared / "deseq2" / "genes.txt",
                "expression_file": self.shared / "deseq2" / "expression.tsv",
                "outdir": self.shared / "seidr",
                "aggregate": "irp",
                "no_full": False
            }
        return self.persisted_cfg.get(tool, {})

    @property
    def cache_high_gb(self) -> int:
        if self.no_cache:
            return 0
        if self.cache_gb is not None:
            return int(self.cache_gb)
        return 300

    @property
    def cache_low_gb(self) -> int:
        high = self.cache_high_gb
        if high <= 0:
            return 0
        return max(1, int(high * 0.8))

    def __repr__(self) -> str:
        return f"<Config output={self.outdir} threads={self.max_threads}>"

# ------------------------------- Sample -------------------------------

class Sample:
    def __init__(
            self,
            sample_id: str,
            sample_type: str,  # "SRR" or "FASTQ"
            fastq_paths: Optional[List[Path]] = None,
            metadata: Optional[Dict[str, Any]] = None,
            bioproject: Optional["BioProject"] = None,
            outdir: Optional[Path] = None,
            status: str = "pending",  # "pending" | "ready" | "processing" | "done" | "failed" | "prefetched"
            cache_dir: Optional[Path] = None,
            no_cache: bool = False,
    ):
        self.id = str(sample_id)
        self.type = sample_type.upper()
        if self.type not in {"SRR", "FASTQ"}:
            raise ValueError(f"Unsupported sample_type: {sample_type}")

        self.bioproject = bioproject

        # -------------------------
        # Decide per-sample outdir
        # -------------------------
        if self.type == "SRR":
            if bioproject is None:
                raise ValueError("SRR sample requires a BioProject.")
            self.outdir = bioproject.path / self.id
        else:
            if outdir is None:
                raise ValueError("FASTQ sample requires an outdir.")
            self.outdir = outdir / self.id

        self.outdir.mkdir(parents=True, exist_ok=True)

        # per-sample log
        self.log_path = self.outdir / f"{self.id}_log.txt"

        # Always ensure directory exists
        self.log_path.parent.mkdir(parents=True, exist_ok=True)

        # Only create header if file doesn't exist
        if not self.log_path.exists():
            with open(self.log_path, "w") as f:
                f.write(f"# Log for Sample {self.id}\n")

        self.fastq_paths: List[Path] = list(fastq_paths or [])
        self.metadata: Dict[str, Any] = dict(metadata or {})
        self.status = status

        # dynamically set during prefetch / detection
        self.sra_path: Optional[Path] = None

        # ------------------------------------------------------
        # Auto-detect state from disk (done / prefetched)
        # ------------------------------------------------------
        if self.type == "SRR":
            # 1) Fully done? (abundance.tsv present)
            if self.is_done():
                # is_done() already sets status="done"
                return

            # 2) Already-prefetched SRA?
            candidates: List[Path] = []

            # cache mode: shared cache/SRR.sra
            if not no_cache and cache_dir is not None:
                candidates.append(cache_dir / f"{self.id}.sra")

            # no-cache or legacy: <BP>/<SRR>/<SRR>.sra
            candidates.append(self.outdir / f"{self.id}.sra")

            for sra_file in candidates:
                try:
                    if sra_file.exists() and sra_file.stat().st_size > 0:
                        self.sra_path = sra_file
                        if self.status not in {"done", "failed"}:
                            self.status = "prefetched"
                        break
                except OSError:
                    # If stat() fails, just ignore and try others
                    continue

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
    def __init__(
            self,
            bioproject_id: str,
            base_outdir: Optional[Path] = None,
            cache_dir: Optional[Path] = None,
            no_cache: bool = False,
    ):
        if base_outdir is None:
            raise ValueError("BioProject requires a base_outdir.")
        self.id = str(bioproject_id)
        self.path = base_outdir / self.id
        self.path.mkdir(parents=True, exist_ok=True)

        self.status: str = "pending"
        self.samples: List["Sample"] = []

        # cache settings (for Sample autodetect)
        self.cache_dir = cache_dir
        self.no_cache = no_cache

        # per-BioProject log
        self.log_path = self.path / f"{self.id}_log.txt"
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.log_path, "w") as f:
            f.write(f"# Log for BioProject {self.id}\n")

    def add_srr(
            self,
            srr_id: str,
            *,
            metadata: Optional[Dict[str, Any]] = None,
            status: str = "pending",
    ) -> "Sample":
        s = Sample(
            sample_id=srr_id,
            sample_type="SRR",
            fastq_paths=None,
            metadata=metadata,
            bioproject=self,
            outdir=None,
            status=status,
            cache_dir=None if self.no_cache else self.cache_dir,
            no_cache=self.no_cache,
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

    def get_sample_ids(self):
        return [sample.id for sample in self.samples]


# ------------------------------- Dataset ------------------------------

class Dataset:
    """
    Single-mode collection of samples for a run.
    mode: "SRR" (grouped by BioProject) or "FASTQ" (local fastqs)
    """

    def __init__(self, config: "Config", mode: str):
        self.config = config
        self.mode = mode.upper()
        if self.mode not in {"SRR", "FASTQ"}:
            raise ValueError("Dataset mode must be 'SRR' or 'FASTQ'.")
        self.path = config.outdir
        self.path.mkdir(parents=True, exist_ok=True)

        # Root for per-sample outputs in FASTQ mode
        self.fastq_root = self.path / "samples"
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
        cache_dir = getattr(self.config, "cache", None)
        no_cache = bool(getattr(self.config, "no_cache", False))
        if bp is None:
            bp = BioProject(
                bioproject_id,
                base_outdir=self.path,
                cache_dir=None if no_cache else cache_dir,
                no_cache=no_cache,
            )
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

    # FASTQ helpers
    def add_fastq(
            self,
            sample_id: str,
            fastq_paths: List[Path],
            *,
            outdir: Optional[Path] = None,
            metadata: Optional[Dict[str, Any]] = None,
            status: str = "pending",
    ) -> Sample:
        """
        FASTQ helpers.

        In FASTQ mode we let Sample create its own directory as:
            <parent_outdir>/<sample_id>/

        So here we just pass the *parent* (typically <outdir>/samples).
        """
        # Parent directory under which samples live
        if outdir is None:
            parent_outdir = self.fastq_root
        else:
            parent_outdir = Path(outdir).resolve()

        s = Sample(
            sample_id=sample_id,
            sample_type="FASTQ",
            fastq_paths=fastq_paths,
            metadata=metadata,
            bioproject=None,
            outdir=parent_outdir,
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
                    s.status = "pending"
        else:
            for bp in self.bioprojects:
                for s in bp.samples:
                    if s.is_done():
                        s.status = "done"
                    else:
                        s.status = "pending"
                bp.update_status()

    def to_do_pairs(self) -> List[Tuple["BioProject", "Sample"]]:
        if self.mode != "SRR":
            return []
        for bp in self.bioprojects:
            bp.update_status()
        ordered_bps = sorted(
            (bp for bp in self.bioprojects if bp.remaining() > 0),
            key=lambda b: (b.remaining(), b.total(), b.id),
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
        """
        Build a SRR-mode dataset from a metadata dataframe with
        columns 'Run', 'BioProject', and optionally 'Model'.
        """
        ds = cls(config=cfg, mode="SRR")
        for _, row in df.iterrows():
            srr = str(row["Run"])
            bp_id = str(row["BioProject"])
            meta = {"Model": row.get("Model", None)} if "Model" in row else {}
            ds.add_srr(bp_id, srr_id=srr, metadata=meta)
        return ds

    @classmethod
    def from_fastq_dir(cls, directory: Path, cfg: "Config") -> "Dataset":
        """
        Build a FASTQ-mode dataset from a directory structured as:

            <root>/
              sampleA/
                sampleA_R1.fastq.gz
                sampleA_R2.fastq.gz
              sampleB/
                sampleB.trim.fastq    # single-end also OK
              emptyFolder/
                # no FASTQs → skipped with a [WARNING]

        Each immediate subdirectory of <root> is treated as one potential sample.
        Subdirectories without FASTQs are skipped with a warning.
        """
        directory = Path(directory).expanduser().resolve()
        if not directory.is_dir():
            raise ValueError(f"FASTQ input must be a directory: {directory}")

        ds = cls(config=cfg, mode="FASTQ")

        sample_dirs = sorted(p for p in directory.iterdir() if p.is_dir())
        if not sample_dirs:
            raise ValueError(
                f"No sample subdirectories found in FASTQ root: {directory}\n"
                "Expected layout: <root>/<sample_id>/(R1/R2 or R.fastq[.gz])"
            )

        for sample_dir in sample_dirs:
            fastqs = sorted(
                p for p in sample_dir.iterdir()
                if _is_fastq_path(p)
            )
            if not fastqs:
                # Just warn and continue; this BioProject folder may contain other stuff
                print(f"[WARNING] No FASTQ files found in sample directory: {sample_dir}")
                continue

            sample_id = sample_dir.name
            ds.add_fastq(sample_id=sample_id, fastq_paths=fastqs)

        if not ds.samples:
            # None of the subdirs had FASTQs → hard error
            raise ValueError(
                f"No FASTQ files found in any subdirectory of: {directory}"
            )

        return ds

    @classmethod
    def reconstruct_from_output(cls, cfg: Config) -> 'Dataset':
        """
        Scans the output directory defined in Config to reconstruct a Dataset object
        based on finished processing results (presence of abundance.tsv).
        """
        if not cfg.outdir.exists():
            raise FileNotFoundError(f"Output directory not found: {cfg.outdir}")

        all_files = cfg.outdir.rglob("abundance.tsv")

        # Exclude unwanted directories
        blacklist_markers = {"multiqc", "mqc", "shared", "cache", "logs", "slurm_logs"}
        found_files = []
        for f in all_files:
            if not any(marker in part for part in f.parts for marker in blacklist_markers):
                found_files.append(f)

        if not found_files:
            raise FileNotFoundError(f"No valid processed samples (abundance.tsv) found in {cfg.outdir}")

        # Determine mode
        detected_mode = "SRR"
        for f in found_files:
            # Check structure for FASTQ mode: .../fastq_samples/<sample_id>/abundance.tsv
            # f.parent = sample_id, f.parent.parent = "fastq_samples"
            if f.parent.parent.name == "fastq_samples":
                detected_mode = "FASTQ"
                break

        ds = cls(config=cfg, mode=detected_mode)

        for f in found_files:
            sample_dir = f.parent
            bp_dir = sample_dir.parent
            s_id = sample_dir.name
            bp_id = bp_dir.name

            if detected_mode == "FASTQ":
                # For FASTQ mode, we rely on add_fastq
                # Pass 'fastq_paths=[]' because we are post-processing,
                # we don't need the original fastqs to run DESeq2.
                ds.add_fastq(sample_id=s_id, fastq_paths=[], outdir=bp_dir, status="done")
            else:
                # For SRR mode, use add_srr to auto-wire the BioProject
                ds.add_srr(bioproject_id=bp_id, srr_id=s_id, status="done")

        return ds

    def __repr__(self) -> str:
        return (
            f"<Dataset mode={self.mode} "
            f"n_samples={len(self.samples)} "
            f"n_bioprojects={len(self.bioprojects)}>"
        )