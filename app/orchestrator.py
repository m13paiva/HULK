from __future__ import annotations

import math
import os
import sys
import time
from dataclasses import dataclass
from threading import Thread
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple

from tqdm.auto import tqdm

from .utils import log, log_err, pad_desc, bp_seidr_batches
from .prefetcher import prefetch
from .processor import process
from .cache_manager import CacheGate
from .qc import run_multiqc, build_bp_metrics
from .post_processing import run_postprocessing_bp, run_postprocessing_batch
from .seidr import run_seidr_aggregation


def force_test_seidr_aggregation(cfg: "Config") -> None:
    """
    A sledgehammer to bypass the orchestrator and test both R aggregation and Seidr
    on a specific multi-BP batch.
    """
    import sys
    from pathlib import Path
    from types import SimpleNamespace
    from .seidr import run_seidr_aggregation
    from .post_processing import run_postprocessing_batch

    print(f"\n{'=' * 50}", flush=True)
    print(f"[TEST MODE] Forcing Post-processing & Seidr on [PRJNA1174730, PRJNA1198503]", flush=True)
    print(f"{'=' * 50}\n", flush=True)

    # 1. Mock the BioProject objects because the function expects objects, not strings
    bp1 = SimpleNamespace(id="PRJNA1174730", path=cfg.outdir / "PRJNA1174730")
    bp2 = SimpleNamespace(id="PRJNA1198503", path=cfg.outdir / "PRJNA1198503")
    fake_batch = [bp1, bp2]

    # 2. Run the post-processing step
    print("\n[TEST MODE] Step 1: Running run_postprocessing_batch...", flush=True)
    batch_dir = run_postprocessing_batch(fake_batch, cfg)

    if not batch_dir or not batch_dir.exists():
        print("\n[TEST MODE ERROR] Post-processing failed or returned nothing. Look at the R output above.", flush=True)
        sys.exit(1)

    print(f"\n[TEST MODE] Post-processing succeeded. Batch dir: {batch_dir}", flush=True)

    # 3. Run Seidr
    print("\n[TEST MODE] Step 2: Running run_seidr_aggregation...", flush=True)
    try:
        run_seidr_aggregation(batch_dir, cfg)
        print("\n[TEST MODE] Seidr finished successfully.", flush=True)
    except Exception as e:
        print(f"\n[TEST MODE] Seidr crashed: {e}", flush=True)
    finally:
        print("[TEST MODE] Exiting pipeline so the real orchestrator doesn't run.", flush=True)
        sys.exit(0)


def _finalize_bioproject(bp: "BioProject", cfg: "Config") -> None:
    """
    Per-BioProject post-processing logic.
    """
    log_path = cfg.log

    # --- EARLY EXIT IF DISABLED ---
    if cfg.no_bp_postprocessing:
        log(f"[{bp.id}] Skipping BioProject post-processing (MultiQC, R, Metrics) due to flag.", log_path)
        return

    errors = cfg.error_warnings

    # 1. MultiQC
    try:
        mqc_data = run_multiqc(bp, cfg, modules=("kallisto", "fastp"))
        if mqc_data is None:
            log(f"[{bp.id}] MultiQC returned no output", log_path)
    except Exception as e:
        log_err(errors, log_path, f"[{bp.id}] MultiQC failed: {e}")

    # 2. R-based Post-processing
    try:
        if getattr(cfg, "tx2gene", None) is not None:
            counts_path = run_postprocessing_bp(bp, cfg)
            if counts_path is not None:
                log(f"[{bp.id}] Gene counts complete: {counts_path}", log_path)
            else:
                log(f"[{bp.id}] R post-processing complete.", log_path)
        else:
            log(f"[{bp.id}] No tx2gene provided; skipping R steps.", log_path)
    except Exception as e:
        log_err(errors, log_path, f"[{bp.id}] R-based BP post-processing failed: {e}")

    # 3. Read Metrics
    try:
        build_bp_metrics(bp, cfg, out_tsv=bp.path / "read_metrics.tsv")
        log(f"[{bp.id}] Read metrics table written.", log_path)
    except Exception as e:
        log_err(errors, log_path, f"[{bp.id}] Failed to build read metrics: {e}")

    # --- MARKER CREATION ---
    try:
        marker = bp.path / ".postprocessing_done"
        marker.touch()
    except Exception as e:
        log_err(errors, log_path, f"[{bp.id}] Failed to create completion marker: {e}")

    log(f"[{bp.id}] === BioProject post-processing complete ===", log_path)


def _start_bp_progress(dataset, cfg, *, start_position: int = 2, poll_secs: float = 0.5):
    """
    Create one tqdm bar per BioProject and return a monitor thread.
    Takes the full dataset object to allow batch organization utilities to inspect dataset attributes.
    """
    bioprojects = dataset.bioprojects
    bp_bars = {}
    last_done = {}
    pos = start_position

    terminal = {"done", "failed", "skipped"}

    # Track which bars are finished and closed to stop updating them
    closed_bars = set()

    # --- BATCH LOGIC INITIALIZATION ---
    # Only organize batches if we are in 'aggregated' or 'both' mode.
    # Defaults to 'single' if missing.
    seidr_mode = getattr(cfg, "seidr_construction_mode", "single")

    batches = []
    if seidr_mode in ("aggregated", "both"):
        batches = bp_seidr_batches(dataset, min_samples=50)
        if batches:
            log(f"[monitor] Organized {len(batches)} batches for aggregated Seidr inference.", cfg.log)

    processed_batches = set()

    # -------------------------------
    # Create bars with initial offset
    # -------------------------------
    for bp in bioprojects:
        total = len(bp.samples)

        # Handle empty BioProjects gracefully
        if total == 0:
            closed_bars.add(bp.id)
            continue

        bar = tqdm(
            total=total,
            desc=pad_desc(bp.id),
            unit="Sample",
            position=pos,
            leave=True,  # Bar remains visible after closing (completed state)
            mininterval=0,
        )

        # initial offset = already completed samples
        init = sum(
            1 for s in bp.samples
            if getattr(s, "status", None) in terminal
        )
        if init:
            bar.update(init)
            bar.refresh()

        bp_bars[bp.id] = bar
        last_done[bp.id] = init
        pos += 1

    # -------------------------------
    # Monitor thread to update bars
    # -------------------------------
    def _monitor():
        # Track which BPs have already had postprocessing run in THIS pipeline run
        already_postprocessed = set()
        force_run = getattr(cfg, "force", False)

        while True:
            # If all bars are closed AND (if batching is active) all batches are done
            if len(closed_bars) == len(bp_bars) and len(processed_batches) == len(batches):
                break

            for bp in bioprojects:
                # Skip bars that are already finished and closed
                if bp.id in closed_bars:
                    continue

                # Safety check if bp wasn't initialized (e.g. empty)
                if bp.id not in bp_bars:
                    continue

                done_now = sum(
                    1 for s in bp.samples
                    if getattr(s, "status", None) in terminal
                )
                delta = done_now - last_done[bp.id]

                if delta > 0:
                    bp_bars[bp.id].update(delta)
                    bp_bars[bp.id].refresh()
                    last_done[bp.id] = done_now

                # All samples terminal -> run postprocessing ONCE per run
                if done_now >= len(bp.samples):
                    # 1. Run Post-processing (if not already done)
                    if bp.id not in already_postprocessed:
                        marker = bp.path / ".postprocessing_done"

                        # CHECK MARKER
                        if marker.exists() and not force_run:
                            # Skip execution, but mark as done so batch logic sees it
                            log(f"[{bp.id}] Found .postprocessing_done marker. Skipping BP post-processing.", cfg.log)
                            bp.status = "done"
                            already_postprocessed.add(bp.id)
                        else:
                            # Execute actual processing
                            bp.status = "done"
                            try:
                                _finalize_bioproject(bp, cfg)
                            except Exception as e:
                                # Log error but don't crash the monitor thread
                                if hasattr(bp, 'log_path'):
                                    log(f"[{bp.id}] Postprocessing failed: {e}", bp.log_path)
                            finally:
                                already_postprocessed.add(bp.id)

                    # 2. Close the bar
                    bp_bars[bp.id].close()
                    closed_bars.add(bp.id)

            # -------------------------------
            # Batch Completion Check
            # -------------------------------
            for idx, batch in enumerate(batches):
                if idx in processed_batches:
                    continue

                # Check if all BioProjects in this specific batch are 'done'
                if all(getattr(bp, "status", None) == "done" for bp in batch):
                    log(f"[batch-monitor] Batch {idx} ready. Aggregating expression data...", cfg.log)

                    # Call the batch post-processing function (R aggregation)
                    batch_dir = run_postprocessing_batch(batch, cfg)

                    if batch_dir:
                        log(f"[batch-monitor] Batch {idx} data prepared at {batch_dir}", cfg.log)

                        # --- TRIGGER SEIDR INFERENCE ---
                        try:
                            run_seidr_aggregation(batch_dir, cfg)
                        except Exception as e:
                            log_err(cfg.error_warnings, cfg.log, f"[batch-monitor] Seidr failed for batch {idx}: {e}")

                    else:
                        log(f"[batch-monitor] Batch {idx} failed data aggregation.", cfg.log)

                    processed_batches.add(idx)

            time.sleep(poll_secs)

        # Ensure all bars are explicitly closed at exit (redundant safety)
        for bid, bar in bp_bars.items():
            if bid not in closed_bars:
                bar.close()

    t = Thread(target=_monitor, daemon=True)
    t.start()

    return bp_bars, t


def _cfg(cfg, name, default=None):
    return getattr(cfg, name, default)


def _clamp(v, lo, hi):
    return max(lo, min(hi, v))


@dataclass
class ThreadPlan:
    bundle_concurrency: int
    bundle_threads: int
    dump_cap: Optional[int]
    fastp_cap: Optional[int]
    kallisto_cap: Optional[int]
    prefetch_workers: int
    logical_cpus: int
    user_max_threads: Optional[int]
    reserved_for_os: int
    usable_threads: int


def _plan_threads(cfg) -> ThreadPlan:
    logical = os.cpu_count() or 8
    reserve_default = _clamp(math.floor(logical * 0.10), 1, min(4, max(1, logical - 1)))
    user_max = _cfg(cfg, "max_threads", _cfg(cfg, "threads", None))
    if user_max is not None:
        user_max = _clamp(int(user_max), 1, logical)
        usable = logical - reserve_default if user_max > logical - reserve_default else user_max
    else:
        usable = logical - reserve_default
    usable = max(1, usable)

    default_bundles = (1 if usable < 8 else 2 if usable < 24 else 3 if usable < 48 else 4 if usable < 96 else 6)
    bundles = _clamp(int(_cfg(cfg, "max_bundles", default_bundles)), 1, usable)

    forced_bt = _cfg(cfg, "bundle_threads", None)
    bt = max(1, int(forced_bt)) if forced_bt is not None else max(1, usable // bundles)

    dump_cap = _cfg(cfg, "dump_cap", 1)
    fastp_cap = _cfg(cfg, "fastp_cap", 8)
    kallisto_cap = _cfg(cfg, "kallisto_cap", 32)

    default_pf = _clamp(max(4, logical // 8), 2, 24)
    pf_workers = int(_cfg(cfg, "prefetch_workers", default_pf))

    return ThreadPlan(
        bundle_concurrency=bundles,
        bundle_threads=bt,
        dump_cap=None if dump_cap is None else int(dump_cap),
        fastp_cap=None if fastp_cap is None else int(fastp_cap),
        kallisto_cap=None if kallisto_cap is None else int(kallisto_cap),
        prefetch_workers=pf_workers,
        logical_cpus=logical,
        user_max_threads=(int(_cfg(cfg, "max_threads")) if _cfg(cfg, "max_threads", None) is not None else None),
        reserved_for_os=reserve_default,
        usable_threads=usable,
    )


def run_download_and_process(
        dataset: "Dataset",
        *,
        cfg: "Config",
        cache_dir: Path,
        work_root: Path,
        temp_dir: Path,
        log_path: Path
):
    """
    Orchestrate the two operational modes:

      CACHE MODE:
        • Prefetch thread runs in parallel
        • Processing daemon consumes cached .sra files

      NO-CACHE MODE (--no-cache):
        • No prefetch thread
        • Processing daemon calls prefetch_one() internally per-SRR
        • Still parallel across bundles
    """
    # ----------------------------------------
    # HIJACK PIPELINE FOR TESTING
    # ----------------------------------------
    #force_test_seidr_aggregation(cfg)

    gate = CacheGate(
        cache_dir,
        cfg.cache_high_gb * (1024 ** 3),
        cfg.cache_low_gb * (1024 ** 3),
    )

    cache_enabled = not gate.disabled
    mode = "cache" if cache_enabled else "local"

    log(
        f"[cache] enabled={cache_enabled}, "
        f"high={getattr(gate, 'high', 0)}, "
        f"low={getattr(gate, 'low', 0)}",
        log_path,
    )

    plan = _plan_threads(cfg)
    msg_user = "auto" if plan.user_max_threads is None else f"{plan.user_max_threads} (user)"

    log(
        (
            f"[orchestrator] CPU: logical={plan.logical_cpus} | reserve(OS)={plan.reserved_for_os} | "
            f"usable={plan.usable_threads} (max={msg_user}) | "
            f"bundles={plan.bundle_concurrency} | bundle_threads={plan.bundle_threads} | "
            f"caps(dump={plan.dump_cap}, fastp={plan.fastp_cap}, kallisto={plan.kallisto_cap}) | "
            f"prefetch_workers={plan.prefetch_workers} | "
            f"cache_mode={'ENABLED' if cache_enabled else 'DISABLED'}"
        ),
        log_path,
    )

    # ----------------------------------------
    # BioProject progress bars
    # ----------------------------------------
    if getattr(dataset, "bioprojects", None):
        _bp_bars, t_bpmon = _start_bp_progress(dataset, cfg, start_position=2, poll_secs=0.5)
    else:
        _bp_bars, t_bpmon = {}, None

    # ----------------------------------------
    # Start PROCESSING daemon
    # ----------------------------------------
    t_proc = Thread(
        target=process,
        kwargs=dict(
            dataset=dataset,
            cfg=cfg,
            cache_dir=cache_dir,
            work_root=work_root,
            temp_dir=temp_dir,
            log_path=log_path,
            max_bundles=plan.bundle_concurrency,
            dump_threads=plan.dump_cap,
            fastp_threads=min(plan.fastp_cap, plan.bundle_threads),
            kallisto_threads=min(plan.kallisto_cap, plan.bundle_threads),
        ),
        daemon=True,
    )
    t_proc.start()

    # ----------------------------------------
    # Start PREFETCH daemon (only in CACHE MODE)
    # ----------------------------------------
    if cache_enabled and dataset.mode == "SRR":
        t_pref = Thread(
            target=prefetch,
            kwargs=dict(
                dataset=dataset,
                cache_dir=cache_dir,
                log_path=log_path,
                max_workers=plan.prefetch_workers,
                cache_high_gb=cfg.cache_high_gb,
                cache_low_gb=cfg.cache_low_gb,
                mode="cache",
            ),
            daemon=True,
        )
        t_pref.start()
    else:
        t_pref = None

    # ----------------------------------------
    # WAIT
    # ----------------------------------------
    if t_pref:
        t_pref.join()

    t_proc.join()

    if t_bpmon:
        t_bpmon.join()

    log("[orchestrator] all done", log_path)