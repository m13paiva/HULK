from __future__ import annotations

import math
import os
import time
from dataclasses import dataclass
from threading import Thread
from pathlib import Path
from typing import Optional

from tqdm.auto import tqdm

from .utils import log, pad_desc
from .prefetcher import prefetch
from .processor import process
from .cache_manager import CacheGate



def _start_bp_progress(bioprojects, cfg, *, start_position: int = 2, poll_secs: float = 0.5):
    """Create one tqdm bar per BioProject and return a monitor thread."""
    bp_bars = {}
    last_done = {}
    pos = start_position

    terminal = {"done", "failed", "skipped"}

    # -------------------------------
    # Create bars with initial offset
    # -------------------------------
    for bp in bioprojects:
        total = len(bp.samples)

        bar = tqdm(
            total=total,
            desc=pad_desc(bp.id),
            unit="Sample",
            position=pos,
            leave=True,
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

        while True:
            all_finished = True
            for bp in bioprojects:
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
                if done_now == len(bp.samples) and bp.id not in already_postprocessed:
                    # optional: keep status for humans
                    bp.status = "done"
                    try:
                        bp.run_postprocessing(cfg)
                    except Exception as e:
                        log(f"[{bp.id}] Postprocessing failed: {e}", bp.log_path)
                    finally:
                        already_postprocessed.add(bp.id)

                if done_now < len(bp.samples):
                    all_finished = False

            if all_finished:
                break

            time.sleep(poll_secs)

        # close at end
        for bar in bp_bars.values():
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

    dump_cap     = _cfg(cfg, "dump_cap", 1)
    fastp_cap    = _cfg(cfg, "fastp_cap", 8)
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
    # Decide cache mode using CacheGate
    # ----------------------------------------
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

    # ----------------------------------------
    # Thread planning
    # ----------------------------------------
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
        _bp_bars, t_bpmon = _start_bp_progress(dataset.bioprojects,cfg, start_position=2, poll_secs=0.5)
    else:
        _bp_bars, t_bpmon = {}, None

    # ----------------------------------------
    # Start PROCESSING daemon
    # ----------------------------------------
    # NOTE: processing daemon must receive cache mode
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
                mode="cache",          # <<< NEW
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

