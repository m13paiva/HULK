# core.py
import sys
import time
import shutil
import subprocess
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from tqdm.auto import tqdm

from .utils import (
    log, log_err,
    run_cmd, run_cmd_stream,
    df_to_dict, get_available_threads, plan_workers,
    is_sra_done, detect_fastq_layout, clean_fastq_files,
    merge_bioproject_tpm, merge_all_bioprojects_tpm,
)
from .trim import (
    fastp_single_cmd,
    fastp_paired_cmd
)

from .align import (
    build_transcriptome_index,
    kallisto_single_cmd,
    kallisto_paired_cmd
)

from .qc import (
    run_multiqc,
    build_bp_metrics
)

from .tx2gene import (
    bp_gene_counts,
    global_gene_counts
)


# ─────────────────────────────────────────────────────────────────────────────
# Per-SRR worker
# ─────────────────────────────────────────────────────────────────────────────

def pipeline_one_srr(run, transcripts_index: Path, threads: int, outdir: Path,
                     log_path: Path, overwrite: bool, keep_fastq: bool,error_warnings: list[str],
                     trim_opts=None,tximport_opts=None) -> None:
    """
    Process one SRR:
      prefetch → fasterq-dump → fastp → kallisto → clean FASTQs
    """
    outdir = outdir.expanduser().resolve()
    run_id, run_model = run

    if outdir.exists():
        if overwrite:
            shutil.rmtree(outdir)
        else:
            if is_sra_done(outdir):
                log(f"SKIP [{run_id}]: {outdir} exists (use --overwrite to redo)", log_path)
                return
            else:
                shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # prefetch: prefer local .sra (raise max cap), fall back to accession streaming
    sra_file = outdir / f"{run_id}.sra"
    prefetch_try1 = ["prefetch", "--force", "all", "-X", "200G", "--type", "sra", "--output-file", str(sra_file), run_id]
    prefetch_try2 = ["prefetch", "--force", "all", "-X", "200G", "--output-directory", ".", run_id]

    def _normalize_sra():
        alt_dir, alt_file = outdir / run_id, outdir / run_id / f"{run_id}.sra"
        if (sra_file is not None and not sra_file.exists()) and alt_file.exists():
            alt_file.replace(sra_file)
            shutil.rmtree(alt_dir, ignore_errors=True)

    for i, cmd in enumerate((prefetch_try1, prefetch_try2), start=1):
        try:
            run_cmd(cmd, outdir, log_path)
            _normalize_sra()
            break
        except subprocess.CalledProcessError as e:
            log(f"[{run_id}] prefetch attempt {i} failed (exit {e.returncode}).", log_path)
            if i == 2:
                log(f"[{run_id}] proceeding without local .sra (will stream via fasterq-dump)", log_path)
                sra_file = None
    _normalize_sra()

    target = str(sra_file) if (sra_file is not None and sra_file.exists()) else run_id
    if target == run_id:
        log(f"WARN: .sra not found; using accession {run_id} (remote/cache mode)", log_path)

    fqd_threads = min(threads, 16)  # v3 is touchy with huge thread counts
    run_cmd([
        "fasterq-dump", target,
        "--skip-technical",
        "--split-files",
        "--threads", str(fqd_threads),
        "--outdir", ".",
        "--temp", ".",
    ], outdir, log_path)

    try:
        if sra_file is not None and sra_file.exists():
            sra_file.unlink()
    except Exception:
        pass

    log(f"✅ Finished downloading {run_id} in {outdir}", log_path)

    layout = detect_fastq_layout(run_id, outdir)
    _opts = {k: v for k, v in (trim_opts or {}).items()
             if k in {"window_size", "mean_quality"} and v is not None}
    if layout[0] == 'SINGLE':
        r = layout[1]
        run_cmd(fastp_single_cmd(outdir, run_id, r, threads,**_opts), outdir, log_path)
        run_cmd_stream(
            kallisto_single_cmd(outdir, run_id, transcripts_index, threads, run_model, log_path, error_warnings),
            outdir, log_path, outdir / f"kallisto_{run_id}.log"
        )

    elif layout[0] == 'PAIRED':
        r1, r2 = layout[1], layout[2]
        run_cmd(fastp_paired_cmd(outdir, run_id, r1, r2, threads,**_opts), outdir, log_path)
        run_cmd_stream(
            kallisto_paired_cmd(outdir, run_id, transcripts_index, threads, log_path, error_warnings),
            outdir, log_path, outdir / f"kallisto_{run_id}.log"
        )

    else:
        log_err(error_warnings, log_path, f"[{run_id}] No FASTQ files detected in {outdir}")
        return

    if not keep_fastq:
        clean_fastq_files(outdir)
    log(f"✅ Finished processing {run_id}", log_path)


# ─────────────────────────────────────────────────────────────────────────────
# Per-BioProject driver
# ─────────────────────────────────────────────────────────────────────────────

def pipeline_one_bioproject(
    bioproject: str,
    sras: list[list[str]],
    outdir: Path,
    log_path: Path,
    min_threads: int,
    max_threads: int,
    transcripts_index: Path,
    overwrite: bool,
    keep_fastq: bool,
    tx2gene: Path | None,
    error_warnings: list[str],
    all_bar, all_durations, disable_pb: bool,
    trim_opts=None,
    tximport_opts=None
) -> None:
    """
    Run all SRRs in one BioProject with auto-scaling threads, then:
      • per-BP MultiQC (sanitized inputs)
      • per-BP read metrics
      • per-BP TPM merge or gene counts (if tx2gene)
    """
    bioproject_dir = outdir / bioproject
    bioproject_dir.mkdir(parents=True, exist_ok=True)

    if overwrite:
        done_in_bp, todo_in_bp = [], sras
    else:
        done_in_bp = [r for r in sras if is_sra_done((outdir / bioproject) / r[0])]
        todo_in_bp = [r for r in sras if r not in done_in_bp]

    # progress for this BP
    bp_bar_inner = tqdm(total=len(sras), desc=f"{bioproject} SRRs", position=2, leave=False,
                        disable=disable_pb, unit="SRR", initial=len(done_in_bp))

    # If nothing to run, still regenerate per-BP artifacts
    if not todo_in_bp:
        mqc_data = run_multiqc(bioproject_dir, log_path, error_warnings)
        if mqc_data is not None:
            _ = build_bp_metrics(bioproject, bioproject_dir, log_path, error_warnings, out_tsv=bioproject_dir / "read_metrics.tsv")
        if tx2gene is not None:
            sra_ids = [r[0] for r in sras]
            bp_gene_counts(bioproject_dir, sra_ids, tx2gene, log_path, error_warnings)
        else:
            merge_bioproject_tpm(bioproject_dir, log_path, error_warnings)
        bp_bar_inner.close()
        return

    # thread planning
    total_threads = max_threads if max_threads else get_available_threads()
    jobs, threads_each = plan_workers(total_threads, len(todo_in_bp), min_threads=min_threads)
    log(f'BIOPROJECT: {bioproject}', log_path)
    log(f"🖥️ Available threads: {total_threads}", log_path)
    log(f"➡️ Running {jobs} jobs in parallel, {threads_each} threads per job", log_path)

    start_times: dict[str, float] = {}
    bp_durations: list[float] = []

    with ThreadPoolExecutor(max_workers=jobs) as ex:
        fut2run = {}
        for run in todo_in_bp:
            run_id = run[0]
            run_dir = bioproject_dir / run_id
            start_times[run_id] = time.time()
            fut = ex.submit(pipeline_one_srr,
                            run,
                            transcripts_index,
                            threads_each,
                            run_dir,
                            log_path,
                            overwrite,
                            keep_fastq,
                            error_warnings,
                            trim_opts,
                            tximport_opts)
            fut2run[fut] = run_id

        for fut in as_completed(fut2run):
            run_id = fut2run[fut]
            dur = time.time() - start_times[run_id]
            bp_durations.append(dur)
            all_durations.append(dur)

            # per-BP ETA
            bp_mean = sum(bp_durations) / max(1, len(bp_durations))
            bp_remaining = len(sras) - (len(done_in_bp) + len(bp_durations))
            bp_bar_inner.set_postfix_str(f"avg {int(bp_mean)}s | ETA {int(bp_mean * bp_remaining)}s")
            bp_bar_inner.update(1)

            # global ETA (simple)
            g_mean = sum(all_durations) / max(1, len(all_durations))
            all_bar.set_postfix_str(f"{all_bar.n + 1}/{all_bar.total} | avg {int(g_mean)}s")
            all_bar.update(1)

            try:
                fut.result()  # re-raise if worker failed
            except Exception as e:
                log_err(error_warnings, log_path, f"[{run_id}] FAILED: {e}")
                continue

    # per-BP MultiQC (sanitized) + metrics
    mqc_data = run_multiqc(bioproject_dir, log_path, error_warnings)
    if mqc_data is not None:
        _ = build_bp_metrics(bioproject, bioproject_dir, log_path, error_warnings, out_tsv=bioproject_dir / "read_metrics.tsv")

    # per-BP gene table or TPM merge
    if tx2gene is not None:
        sra_ids = [run[0] for run in sras]
        bp_gene_counts(bioproject_dir, sra_ids, tx2gene, log_path, error_warnings,tximport_opts)
    else:
        merge_bioproject_tpm(bioproject_dir, log_path, error_warnings)

    bp_bar_inner.close()


# ─────────────────────────────────────────────────────────────────────────────
# Whole-run orchestrator
# ─────────────────────────────────────────────────────────────────────────────

def pipeline(
    data_df: pd.DataFrame,
    outdir: Path,
    transcriptome: Path,
    min_threads: int,
    max_threads: int,
    verbose: bool,
    overwrite: bool,
    keep_fastq:bool,
    overall_table: bool,
    tx2gene: Path | None,
    trim_opts=None,
    tximport_opts=None,
) -> None:
    """
    Top-level runner:
      • build index if needed
      • iterate BioProjects
      • per-BP pipeline + per-BP MultiQC + metrics
      • global MultiQC + optional overall tables
    """
    outdir = outdir.expanduser().resolve()
    error_warnings: list[str] = []

    data = df_to_dict(data_df)

    shared = (outdir / 'shared').expanduser().resolve()
    shared.mkdir(parents=True, exist_ok=True)
    log_path = (shared / "log.txt").expanduser().resolve()

    with open(log_path, "w", buffering=1) as f:
        f.write(f"\n===== HULK start {datetime.now().isoformat()} =====\n")

    log(f"Loaded {len(data_df)} entries", log_path)
    log(f"First few rows:\n{data_df.head()}", log_path)
    log(f"Output directory: {str(outdir)}", log_path)

    # global progress bar
    total_sras = sum(len(v) for v in data.values())
    already_done = 0 if overwrite else sum(
        1 for bp, runs in data.items() for run in runs if is_sra_done(outdir / bp / run[0])
    )
    disable_pb = (not verbose) or (not sys.stdout.isatty())
    all_bar = tqdm(total=total_sras, desc="All SRRs", position=1, leave=True,
                   disable=disable_pb, unit="SRR", initial=already_done)
    all_durations: list[float] = []

    # index
    if transcriptome.suffix.lower() != '.idx':
        transcripts_index = build_transcriptome_index(transcriptome, shared, log_path)
    else:
        transcripts_index = transcriptome

    # per-BP
    for bioproject, sras in data.items():
        pipeline_one_bioproject(
            bioproject, sras, outdir, log_path,
            min_threads, max_threads, transcripts_index,
            overwrite, keep_fastq,tx2gene, error_warnings,
            all_bar, all_durations, disable_pb,trim_opts,tximport_opts
        )

    # global MultiQC (fastp + kallisto)
    run_multiqc(outdir, shared, "multiqc_shared", log_path, modules=('kallisto','fastp'))

    # overall tables
    if overall_table:
        if tx2gene is not None:
            global_gene_counts(outdir, tx2gene, log_path, error_warnings, shared / "global_gene_counts.tsv",tximport_opts)
        else:
            merge_all_bioprojects_tpm(outdir, log_path, error_warnings, shared / "overall_TPM.tsv")

    if error_warnings:
        log('\n\n===================WARNINGS===================\n', log_path)
        for m in error_warnings:
            log(m, log_path)

    all_bar.close()
