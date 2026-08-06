# orchestrator.py

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

from .utils import log, pad_desc
from .prefetcher import prefetch
from .processor import process
from .cache_manager import CacheGate
from .qc import run_multiqc, build_bp_metrics
from .post_processing import run_postprocessing_bp


def _finalize_bioproject(bp: "BioProject", cfg: "Config") -> None:
    """
    Executes specific post-processing validation logic isolating metrics calculations locally
    per independent subset. Resolves dependency blocks.

    Args:
        bp (BioProject): State representation tracking a subset container.
        cfg (Config): Running profile logic.
    """
    log_path = cfg.log

    if cfg.no_bp_postprocessing:
        log(f"[{bp.id}] Skipping BioProject post-processing (MultiQC, R, Metrics) due to flag.", log_path)
        return

    # 1. MultiQC Initialization
    try:
        mqc_data = run_multiqc(bp, cfg, modules=("kallisto", "fastp"))
        if mqc_data is None:
            log(f"[{bp.id}] MultiQC returned no output", log_path)
    except Exception as e:
        print(f"[{bp.id}] MultiQC failed: {e}")

    # 2. R-Script Matrix Calculation and Resolution Processing
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
        print(f"[{bp.id}] R-based BP post-processing failed: {e}")

    # 3. Read Depth/Metric Calculation Module
    try:
        build_bp_metrics(bp, cfg, out_tsv=bp.path / "read_metrics.tsv")
        log(f"[{bp.id}] Read metrics table written.", log_path)
    except Exception as e:
        print(f"[{bp.id}] Failed to build read metrics: {e}")

    # Explicit marker creation guaranteeing persistence across interrupted runs
    try:
        marker = bp.path / ".postprocessing_done"
        marker.touch()
    except Exception as e:
        print(f"[{bp.id}] Failed to create completion marker: {e}")

    log(f"[{bp.id}] === BioProject post-processing complete ===", log_path)


def _start_bp_progress(dataset, cfg, *, start_position: int = 2, poll_secs: float = 0.5):
    """
    Allocates thread logic controlling CLI terminal bar outputs syncing standard events.

    Args:
        dataset (Dataset): Target group logic array.
        cfg (Config): Execution context maps.
        start_position (int): Tqdm offset rendering positioning relative to screen origins.
        poll_secs (float): Time threshold validating check delays against CPU burn loops.

    Returns:
        tuple: (bp_bars, thread reference).
    """
    bioprojects = dataset.bioprojects
    bp_bars = {}
    last_done = {}
    pos = start_position

    terminal = {"done", "failed", "skipped"}
    closed_bars = set()

    # Create bars with initial offset mapping resolving previously completed jobs
    for bp in bioprojects:
        total = len(bp.samples)
        if total == 0:
            closed_bars.add(bp.id)
            continue

        bar = tqdm(
            total=total,
            desc=pad_desc(bp.id),
            unit="Sample",
            position=pos,
            leave=True,
            mininterval=0,
        )

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

    def _monitor():
        """Isolates background thread logic probing array resolution dynamically."""
        already_postprocessed = set()
        force_run = getattr(cfg, "force", False)

        while True:
            if len(closed_bars) == len(bp_bars):
                break

            for bp in bioprojects:
                if bp.id in closed_bars or bp.id not in bp_bars:
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

                if done_now >= len(bp.samples):
                    if bp.id not in already_postprocessed:
                        marker = bp.path / ".postprocessing_done"

                        if marker.exists() and not force_run:
                            log(f"[{bp.id}] Found .postprocessing_done marker. Skipping BP post-processing.", cfg.log)
                            bp.status = "done"
                            already_postprocessed.add(bp.id)
                        else:
                            bp.status = "done"
                            try:
                                _finalize_bioproject(bp, cfg)
                            except Exception as e:
                                if hasattr(bp, 'log_path'):
                                    log(f"[{bp.id}] Postprocessing failed: {e}", bp.log_path)
                            finally:
                                already_postprocessed.add(bp.id)

                    bp_bars[bp.id].close()
                    closed_bars.add(bp.id)

            time.sleep(poll_secs)

        for bid, bar in bp_bars.items():
            if bid not in closed_bars:
                bar.close()

    t = Thread(target=_monitor, daemon=True)
    t.start()

    return bp_bars, t


def _cfg(cfg, name, default=None):
    """Simple wrapper enforcing attribute safe resolution."""
    return getattr(cfg, name, default)


def _clamp(v, lo, hi):
    """Enforces mathematical boundary limits ensuring numbers stay within hard boundaries."""
    return max(lo, min(hi, v))


@dataclass
class ThreadPlan:
    """
    Data class structure representing mapped hardware allocation limits.
    """
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
    """
    Dynamic hardware inspection method converting system potential to threaded sub-allocations.

    Args:
        cfg (Config): Runtime setup parameters limiting ceiling potential.

    Returns:
        ThreadPlan: Data object explicitly tracking thread allocation blocks.
    """
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
    Primary parallel manager. Controls external fetching against physical processing pipelines.

    Args:
        dataset (Dataset): Reference target mapping structure.
        cfg (Config): Runtime setup constraints.
        cache_dir (Path): Bounded storage layer resolving prefetch requests.
        work_root (Path): Output structure location representing root boundaries.
        temp_dir (Path): Volatile temporary path routing pipeline artifacts.
        log_path (Path): Log structure location routing thread messages.
    """
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

    if getattr(dataset, "bioprojects", None):
        _bp_bars, t_bpmon = _start_bp_progress(dataset, cfg, start_position=2, poll_secs=0.5)
    else:
        _bp_bars, t_bpmon = {}, None

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

    if t_pref:
        t_pref.join()

    t_proc.join()

    if t_bpmon:
        t_bpmon.join()

    log("[orchestrator] all done", log_path)