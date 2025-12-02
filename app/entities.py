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
        # NEW: sequencing technology (FASTQ mode)
        seq_tech: Optional[str] = None,
        # Existing feature / tool config
        trim_window_size: int = 4,
        trim_mean_quality: int = 20,
        align_method: str = "kallisto",
        kallisto_bootstrap: int = 100,
        tximport_mode: str = "raw_counts",
        tximport_ignore_tx_version: bool = False,
        cache_gb: Optional[int] = None,    # -c/--cache (GiB)
        no_cache: bool = False,            # --no-cache
        # NEW: DESeq2 / expression matrix / plotting options (used by Seidr subcmd)
        deseq2_vst_enabled: bool = True,   # main toggle for DESeq2+VST stage
        expr_use_matrix: str = "vst",      # "vst" or "normalized"
        drop_nonvarying_genes: bool = True,
        plot_pca: bool = False,
        plot_heatmap: bool = False,
        plot_var_heatmap: bool = False,
        plots_only_mode: bool = False,
        tximport_only_mode: bool = False,
        target_genes_file: Optional[Path] = None,
        bioproject_filter: Optional[str] = None,
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

        # Sequencing tech (mostly relevant in FASTQ mode; passed from CLI)
        self.seq_tech: Optional[str] = seq_tech

        # --------------------------- Persistent JSON ---------------------------
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

        # DESeq2 / expression exports / plotting (used by R script & Seidr subcommand)
        self.deseq2_vst_enabled: bool = bool(deseq2_vst_enabled)

        self.expr_use_matrix: str = (expr_use_matrix or "vst").lower()
        if self.expr_use_matrix not in {"vst", "normalized"}:
            self.expr_use_matrix = "vst"

        self.drop_nonvarying_genes: bool = bool(drop_nonvarying_genes)
        self.plot_pca: bool = bool(plot_pca)
        self.plot_heatmap: bool = bool(plot_heatmap)
        self.plot_var_heatmap: bool = bool(plot_var_heatmap)
        self.plots_only_mode: bool = bool(plots_only_mode)
        self.tximport_only_mode: bool = bool(tximport_only_mode)
        self.target_genes_file = Path(target_genes_file).expanduser().resolve() if target_genes_file else None
        self.bioproject_filter = bioproject_filter

        # Bundle for the Seidr subcommand / expression-network step
        self.expr_network_opts: Dict[str, Any] = {
            "deseq2_vst_enabled": self.deseq2_vst_enabled,
            "use_matrix": self.expr_use_matrix,
            "drop_nonvarying": self.drop_nonvarying_genes,
            "pca": self.plot_pca,
            "heatmap": self.plot_heatmap,
            "var_heatmap": self.plot_var_heatmap,
            "plots_only": self.plots_only_mode,
            "tximport_only": self.tximport_only_mode,
            "target_genes": str(self.target_genes_file) if self.target_genes_file else None,
            "bioproject": self.bioproject_filter,
        }

        # cache
        self.cache_gb: Optional[int] = cache_gb  # may be None → auto
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
            # seq_tech is optional; align/quant code can branch on it if needed
            return {
                "method": self.align_method,
                "bootstrap": self.kallisto_bootstrap,
                "seq_tech": self.seq_tech,
            }
        if tool == "seidr":
            # Options used by the Seidr subcommand / R pipeline
            return dict(self.expr_network_opts)
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
            f"Seq. tech:    {self.seq_tech or '-'}",
            f"Trim:         ws={self.trim_window_size}, mq={self.trim_mean_quality}",
            f"Align:        method={self.align_method}, boot={self.kallisto_bootstrap}",
            f"tximport:     mode={self.tximport_mode}, ignore_ver={self.tximport_ignore_tx_version}",
            f"DESeq2/VST:   enabled={self.deseq2_vst_enabled}, "
            f"matrix={self.expr_use_matrix}, drop_nonzvar={self.drop_nonvarying_genes}",
            f"Plots:        PCA={self.plot_pca}, HM={self.plot_heatmap}, "
            f"VarHM={self.plot_var_heatmap}, plots_only={self.plots_only_mode}, "
            f"txi_only={self.tximport_only_mode}",
        ]
        if self.bioproject_filter or self.target_genes_file:
            lines.append(
                f"Expr scope:   bioproject={self.bioproject_filter or '-'}, "
                f"target_genes={self.target_genes_file or '-'}"
            )
        lines.append(
            f"Cache:        no_cache={self.no_cache}, high={self.cache_high_gb} GiB, low={self.cache_low_gb} GiB",
        )
        if self.tx2gene:
            lines.append(f"tx2gene:      {self.tx2gene}")
        return "\n".join(lines)

    def prepare_directories(self) -> None:
        if self.outdir.exists() and self.force:
            shutil.rmtree(self.outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)

        if self.shared.exists() and self.force:
            shutil.rmtree(self.shared)
        self.shared.mkdir(parents=True, exist_ok=True)

        if self.cache.exists() and self.force:
            shutil.rmtree(self.cache)
        self.cache.mkdir(parents=True, exist_ok=True)

        with open(self.log, "w", encoding="utf-8") as f:
            f.write(f"\n===== HULK start {datetime.now().isoformat()} =====\n")

    def __repr__(self) -> str:
        return f"<Config output={self.outdir} threads={self.min_threads}-{self.max_threads}>"


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
        if not self.log_path.exists():
            self.log_path.parent.mkdir(parents=True, exist_ok=True)
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

    def run_postprocessing(
            self,
            cfg: "Config"
    ) -> None:
        """
        Per-BioProject post-processing:
            • MultiQC (fastp + kallisto)
            • R-based expression post-processing (tximport / DESeq2 / plots)
            • Read metrics (fastp + run_info.json)
        """
        bp_dir = self.path
        log_path = cfg.log
        errors = cfg.error_warnings
        if not bp_dir.exists():
            from .utils import log
            log(f"[{self.id}] No output directory found; skipping post-processing", log_path)
            return

        # Internal imports — prevent circular import loops
        from .utils import log, log_err
        from .qc import run_multiqc, build_bp_metrics
        from .post_processing import run_postprocessing_bp

        run_ids = [s.id for s in self.samples]

        # ─────────────────────────────────────────
        # 1) MultiQC
        # ─────────────────────────────────────────
        try:
            mqc_data = run_multiqc(
                self,
                cfg,
                modules=("kallisto", "fastp"),
            )
            if mqc_data is None:
                log(f"[{self.id}] MultiQC returned no output", log_path)
        except Exception as e:
            log_err(errors, log_path, f"[{self.id}] MultiQC failed: {e}")
            mqc_data = None

        # ─────────────────────────────────────────
        # 2) R-based BP-level expression post-processing
        #     - tximport-only if DESeq2 is off / tximport-only mode
        #     - full DESeq2+VST+plots ("full treatment") if DESeq2 enabled
        # ─────────────────────────────────────────
        try:
            if getattr(cfg, "tx2gene", None) is not None:
                counts_path = run_postprocessing_bp(self, cfg)
                if counts_path is not None:
                    log(f"[{self.id}] Gene counts complete: {counts_path}", log_path)
                else:
                    # Full DESeq2/VST mode or no counts file; still a success if R didn’t error.
                    log(f"[{self.id}] R post-processing complete (no dedicated gene_counts.tsv).", log_path)
            else:
                log(f"[{self.id}] No tx2gene provided; skipping R-based post-processing.", log_path)
        except Exception as e:
            log_err(errors, log_path, f"[{self.id}] R-based BP post-processing failed: {e}")

        # ─────────────────────────────────────────
        # 3) Read metrics table
        # ─────────────────────────────────────────
        try:
            build_bp_metrics(
                self,
                cfg,
                out_tsv=bp_dir / "read_metrics.tsv",
            )
            log(f"[{self.id}] Read metrics table written.", log_path)
        except Exception as e:
            log_err(errors, log_path, f"[{self.id}] Failed to build read metrics: {e}")

        log(f"[{self.id}] === BioProject post-processing complete ===", log_path)


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
        Register a FASTQ-based sample.

        fastq_paths should contain 1 file (single-end) or 2 files (paired-end);
        they point to the original FASTQ files (under the user-provided FASTQ root).
        """
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
                sampleB_R.fastq.gz    # single-end

        Each immediate subdirectory of <root> is treated as one sample.
        All FASTQ files directly inside that subdirectory are attached
        to that sample's fastq_paths.
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
                if p.is_file() and _is_fastq(p)
            )
            if not fastqs:
                # You can choose to skip empty dirs instead of raising:
                raise ValueError(
                    f"No FASTQ files found in sample directory: {sample_dir}"
                )

            # sample_id is the folder name
            sample_id = sample_dir.name
            ds.add_fastq(sample_id=sample_id, fastq_paths=fastqs)

        if not ds.samples:
            raise ValueError(
                f"No FASTQ files found in any subdirectory of: {directory}"
            )

        return ds

    def __repr__(self) -> str:
        return (
            f"<Dataset mode={self.mode} "
            f"n_samples={len(self.samples)} "
            f"n_bioprojects={len(self.bioprojects)}>"
        )
