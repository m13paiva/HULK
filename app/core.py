from __future__ import annotations


import shutil
import click
from pathlib import Path
from datetime import datetime
from .utils import log, log_err
from .align import build_transcriptome_index
from .qc import run_multiqc_global
from .entities import Config, Dataset
from .orchestrator import run_download_and_process
from .post_processing import run_postprocessing as run__postprocessing
from .seidr import run_seidr


def prepare_runtime_environment(cfg: Config, dataset: Dataset) -> None:
    """
    Sets up the output directory structure and log files.
    """
    import shutil # ensure import is available

    # 1. The Purge (rem-missing-bps) - UNCHANGED logic
    if cfg.rem_missing_bps and cfg.outdir.exists():
        if getattr(dataset, "mode", "SRA") != "SRA":
            print("[WARNING] --rem-missing-bps ignored: Only applicable in SRA mode.") # Changed log call to print/click
        else:
            expected_bps = set(getattr(dataset, "bioprojects", []))
            if not expected_bps:
                print("\n[SAFETY ABORT] Input table empty. Skipping cleanup.")
            else:
                print(f"\n[CLEANUP] Scanning {cfg.outdir}...")
                blacklist = {"shared", "fastq_samples", "logs", "slurm_logs", "multiqc_data"}
                existing_items = [p for p in cfg.outdir.iterdir() if p.is_dir()]

                for folder in existing_items:
                    if folder.name in blacklist or folder.name.endswith("_mqc"):
                        continue
                    if folder.name not in expected_bps:
                        print(f"[DANGER] Deleting extraneous BioProject folder: {folder.name}")
                        if not cfg.dry_run:
                            try:
                                shutil.rmtree(folder)
                            except Exception as ex:
                                print(f"Failed to remove {folder}: {ex}")

    # 2. Standard Directory Setup
    cfg.outdir.mkdir(parents=True, exist_ok=True)

    # --- CHANGED BLOCK START ---
    # We DO NOT delete 'shared' on force anymore.
    # This preserves the cache and the main log history.
    cfg.shared.mkdir(parents=True, exist_ok=True)
    cfg.cache.mkdir(parents=True, exist_ok=True)
    # --- CHANGED BLOCK END ---

    # Initialize Log
    if not cfg.dry_run:
        # Open in Append mode ('a') so we don't wipe history even on force
        mode = 'a'
        with open(cfg.log, mode, encoding="utf-8") as f:
            f.write(f"\n\n{'='*60}\n")
            f.write(f"===== HULK start {datetime.now().isoformat()} =====\n")
            if cfg.force:
                 f.write("!! Run mode: FORCE (Overwriting sample data) !!\n")
            if cfg.rem_missing_bps:
                f.write("!! Run mode: REM_MISSING_BPS (Destructive Cleanup) !!\n")
            f.write(f"{'='*60}\n")

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
    if not cfg.no_global_postprocessing:
        try:
            log("Generating Global MultiQC report...", log_path)
            # Pass 'dataset' to trigger global aggregation mode
            run_multiqc_global(outdir, shared, "multiqc_shared", log_path, modules=("kallisto", "fastp"))
        except Exception as e:
            log_err(cfg.error_warnings, log_path, f"Global MultiQC failed: {e}")
    else:
        log("Skipping Global MultiQC (--no-global-postprocessing).", log_path)


    # Post-processing (R-based: tximport + DESeq2/VST + plots/exports)
    if getattr(cfg, "tx2gene", None) is not None:
        # This will respect cfg.deseq2_vst_enabled and the plot flags.
        # If DESeq2 is disabled, it will automatically fall back to tximport-only.
        run__postprocessing(data, cfg, skip_bp=True)

    try:
        run_seidr(cfg)
    except Exception as e:
        log_err(cfg.error_warnings, log_path, f"[Seidr] Pipeline step failed: {e}")

        # End of pipeline
    log("Pipeline finished.", log_path)

    # Warnings, if any
    if getattr(cfg, "error_warnings", None):
        log("\n\n=================== WARNINGS ===================\n", log_path)
        for m in cfg.error_warnings:
            log(m, log_path)
