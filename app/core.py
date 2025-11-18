from __future__ import annotations

from pathlib import Path

from .utils import log
from .align import build_transcriptome_index
from .qc import run_multiqc_global
from .tx2gene import global_gene_counts
from .entities import Config, Dataset
from .orchestrator import run_download_and_process


def pipeline(data: "Dataset", cfg: "Config") -> None:
    """
    Top-level runner:
      • (re)build index if needed & inject into samples
      • orchestrate prefetch + processing
      • run global MultiQC
      • optional aggregation (TPM / gene counts)
    """
    cfg.prepare_directories()
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

    # Aggregation
    if getattr(cfg, "aggregate", False):
        if getattr(cfg, "tx2gene", None) is not None:
            global_gene_counts(
                outdir,
                cfg.tx2gene,
                log_path,
                cfg.error_warnings,
                shared / "global_gene_counts.tsv",
                tximport_opts,
            )
        else:
            from .utils import merge_all_bioprojects_tpm
            merge_all_bioprojects_tpm(
                outdir,
                log_path,
                cfg.error_warnings,
                shared / "overall_TPM.tsv",
            )

    # Warnings, if any
    if getattr(cfg, "error_warnings", None):
        log("\n\n=================== WARNINGS ===================\n", log_path)
        for m in cfg.error_warnings:
            log(m, log_path)
