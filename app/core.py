# core.py
from __future__ import annotations

import shutil
from pathlib import Path
from datetime import datetime
from .utils import log, generate_read_metrics_plot
from .align import build_transcriptome_index
from .qc import run_multiqc_global
from .entities import Config, Dataset
from .orchestrator import run_download_and_process
from .post_processing import run_postprocessing as run__postprocessing
from .seidr import run_seidr_single

def prepare_runtime_environment(cfg: Config, dataset: Dataset) -> None:
    """
    Sets up the output directory structure and log files.
    """

    # 1. The Purge (rem-missing-bps)
    if cfg.rem_missing_bps and cfg.outdir.exists():
        if getattr(dataset, "mode", "SRR") != "SRR":
            print("[WARNING] --rem-missing-bps ignored: Only applicable in SRA mode.")
        else:
            raw_bps = getattr(dataset, "bioprojects", [])
            expected_bps = {str(bp.id) for bp in raw_bps if hasattr(bp, "id")}

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
    cfg.shared.mkdir(parents=True, exist_ok=True)
    cfg.cache.mkdir(parents=True, exist_ok=True)

    # Initialize Log
    if not cfg.dry_run:
        mode = 'a'
        with open(cfg.log, mode, encoding="utf-8") as f:
            f.write(f"\n\n{'=' * 60}\n")
            f.write(f"===== HULK start {datetime.now().isoformat()} =====\n")
            if cfg.force:
                f.write("!! Run mode: FORCE (Overwriting sample data) !!\n")
            if cfg.rem_missing_bps:
                f.write("!! Run mode: REM_MISSING_BPS (Destructive Cleanup) !!\n")
            f.write(f"{'=' * 60}\n")

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

    if reference and reference.suffix.lower() != ".idx":
        transcripts_index = build_transcriptome_index(reference, shared, log_path)
    else:
        transcripts_index = reference

    transcripts_index = Path(transcripts_index).resolve()
    log(f"Using index: {transcripts_index}", log_path)

    if getattr(data, "mode", None) == "FASTQ":
        samples_iter = getattr(data, "samples", [])
    else:
        samples_iter = [s for bp in getattr(data, "bioprojects", []) for s in bp.samples]

    for s in samples_iter:
        s.metadata.setdefault("kallisto_index", str(transcripts_index))

    if getattr(data, "mode", None) == "SRR":
        for bp in getattr(data, "bioprojects", []):
            log(f"[PLAN] BioProject {bp.id}: {len(bp.samples)} SRR(s) scheduled.", log_path)
    else:
        for s in getattr(data, "samples", []):
            log(f"[PLAN] FASTQ sample {s.id} -> {s.outdir}", log_path)

    # Run orchestrator
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
            run_multiqc_global(outdir, shared, "multiqc_shared", log_path, modules=("kallisto", "fastp"))
            generate_read_metrics_plot(data,cfg.shared / "plots", cfg.log)
        except Exception as e:
            print(f"Global MultiQC failed: {e}")
    else:
        log("Skipping Global MultiQC (--no-global-postprocessing).", log_path)


    # Post-processing
    if getattr(cfg, "tx2gene", None) is not None:
        run__postprocessing(data, cfg, skip_bp=True)

    try:
        run_seidr_single(cfg)
    except Exception as e:
        print(f"[Seidr] Single Network pipeline failed: {e}")

    log("Pipeline finished.", log_path)

    if getattr(cfg, "error_warnings", None):
        log("\n\n=================== WARNINGS ===================\n", log_path)
        for m in cfg.error_warnings:
            log(m, log_path)