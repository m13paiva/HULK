from __future__ import annotations


import shutil
import click
from pathlib import Path
from datetime import datetime
from .utils import log
from .align import build_transcriptome_index
from .qc import run_multiqc_global
from .entities import Config, Dataset
from .orchestrator import run_download_and_process
from .post_processing import run_postprocessing as run__postprocessing


def prepare_runtime_environment(cfg: Config, dataset: Dataset) -> None:
    """
    Sets up the output directory structure and log files.

    If 'rem-missing-bps' is enabled, this function performs a targeted
    purge of BioProject folders in the output directory that do not match
    the current input dataset.
    """

    # ------------------------------------------------------------------
    # 1. The Purge (rem-missing-bps)
    # ------------------------------------------------------------------
    # We do this BEFORE creating new directories to ensure a clean slate,
    # but only if the output dir actually exists.
    if cfg.rem_missing_bps and cfg.outdir.exists():

        # SRA Mode check: FASTQ mode usually doesn't output BioProject folders like this.
        if getattr(dataset, "mode", "SRA") != "SRA":
            log("[WARNING] --rem-missing-bps ignored: Only applicable in SRA mode.", None)
        else:
            expected_bps = set(getattr(dataset, "bioprojects", []))

            # SAFETY CHECK: If the input table is empty, do NOT wipe the folder.
            if not expected_bps:
                click.secho(
                    "\n[SAFETY ABORT] Input table contains no BioProjects. "
                    "Skipping --rem-missing-bps to prevent total deletion of output directory.",
                    fg="red", bold=True
                )
            else:
                click.secho(f"\n[CLEANUP] Scanning {cfg.outdir} for extraneous BioProjects...", fg="yellow")

                # Folders we NEVER delete automatically
                blacklist = {"shared", "fastq_samples", "logs", "slurm_logs", "multiqc_data"}

                # List existing directories
                existing_items = [p for p in cfg.outdir.iterdir() if p.is_dir()]

                for folder in existing_items:
                    # Skip blacklisted or MultiQC folders
                    if folder.name in blacklist or folder.name.endswith("_mqc"):
                        continue

                    if folder.name not in expected_bps:
                        msg = f"[DANGER] Deleting extraneous BioProject folder: {folder.name}"
                        click.secho(msg, fg="red", bold=True)

                        if not cfg.dry_run:
                            try:
                                shutil.rmtree(folder)
                            except Exception as ex:
                                click.secho(f"Failed to remove {folder}: {ex}", fg="red")

    # ------------------------------------------------------------------
    # 2. Standard Directory Setup (Moved from Config)
    # ------------------------------------------------------------------

    # Create Root Output
    cfg.outdir.mkdir(parents=True, exist_ok=True)

    # Clean/Create Shared Directory
    if cfg.shared.exists() and cfg.force:
        # If forcing, we might want to clear shared logs/cache,
        # though usually we just want to overwrite.
        # Kept strict cleaning as per your previous code:
        try:
            shutil.rmtree(cfg.shared)
        except Exception as e:
            click.secho(f"Warning: Could not clear shared dir: {e}", fg="yellow")

    cfg.shared.mkdir(parents=True, exist_ok=True)

    # Clean/Create Cache Directory
    if cfg.cache.exists() and cfg.force:
        try:
            shutil.rmtree(cfg.cache)
        except Exception:
            pass
    cfg.cache.mkdir(parents=True, exist_ok=True)

    # Initialize Log
    # We do this last to ensure the 'shared' dir exists
    if not cfg.dry_run:
        with open(cfg.log, "w", encoding="utf-8") as f:
            f.write(f"\n===== HULK start {datetime.now().isoformat()} =====\n")
            if cfg.rem_missing_bps:
                f.write("!! WARNING: Run started with --rem-missing-bps (Destructive Mode) !!\n")

def pipeline(data: "Dataset", cfg: "Config") -> None:

    prepare_runtime_environment(cfg,data)
    data.update_status()
    outdir: Path = cfg.outdir
    shared: Path = cfg.shared
    cache_dir: Path = cfg.cache
    log_path: Path = cfg.log
    reference: Path = cfg.reference_path
    tximport_opts = getattr(cfg, "tximport_opts", {})

    temp_dir: Path = getattr(cfg, "temp_dir", shared / "tmp")
    temp_dir.mkdir(parents=True, exist_ok=True)

    # Basic header
    log(f"Output directory: {outdir}", log_path)
    log(f"Mode: {getattr(data, 'mode', '-')}", log_path)

    total_samples = len(data)
    total_done = 0 if getattr(cfg, "force", False) else len(getattr(data, "done", lambda: [])())
    log(f"Total samples: {total_samples} | done: {total_done}", log_path)

    if getattr(data, "mode", None) == "SRR":
        bp_total = len(getattr(data, "bioprojects", []))
        bp_done = len(getattr(data, "bp_done", lambda: [])())
        bp_ids = [bp.id for bp in getattr(data, "bioprojects", [])]
        log(f"BioProjects ({bp_total}, done {bp_done}): {', '.join(sorted(bp_ids))}", log_path)

    # Build / locate kallisto index
    if reference and reference.suffix.lower() != ".idx":
        transcripts_index = build_transcriptome_index(reference, shared, log_path)
    else:
        transcripts_index = reference

    transcripts_index = Path(transcripts_index).resolve()
    log(f"Using index: {transcripts_index}", log_path)

    # Inject index path into sample metadata
    if getattr(data, "mode", None) == "FASTQ":
        samples_iter = getattr(data, "samples", [])
    else:
        samples_iter = [s for bp in getattr(data, "bioprojects", []) for s in bp.samples]

    for s in samples_iter:
        s.metadata.setdefault("kallisto_index", str(transcripts_index))

    # Log per-project plan
    if getattr(data, "mode", None) == "SRR":
        for bp in getattr(data, "bioprojects", []):
            log(f"[PLAN] BioProject {bp.id}: {len(bp.samples)} SRR(s) scheduled.", log_path)
    else:
        for s in getattr(data, "samples", []):
            log(f"[PLAN] FASTQ sample {s.id} -> {s.outdir}", log_path)

    # Run orchestrator (prefetch + processing) — uses cfg internally (incl. bootstraps)
    run_download_and_process(
        dataset=data,
        cfg=cfg,
        cache_dir=cache_dir,
        work_root=outdir,
        temp_dir=temp_dir,
        log_path=log_path,
    )

    # Global MultiQC
    run_multiqc_global(outdir, shared, "multiqc_shared", log_path, modules=("kallisto", "fastp"))

    # Post-processing (R-based: tximport + DESeq2/VST + plots/exports)
    if getattr(cfg, "tx2gene", None) is not None:
        # This will respect cfg.deseq2_vst_enabled and the plot flags.
        # If DESeq2 is disabled, it will automatically fall back to tximport-only.
        run__postprocessing(data, cfg, skip_bp=True)

    # Warnings, if any
    if getattr(cfg, "error_warnings", None):
        log("\n\n=================== WARNINGS ===================\n", log_path)
        for m in cfg.error_warnings:
            log(m, log_path)
