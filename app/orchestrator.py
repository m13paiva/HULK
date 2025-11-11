# orchestrator.py
from __future__ import annotations
import os, math
from dataclasses import dataclass
from threading import Thread
from pathlib import Path

from .utils import log
from .prefetcher import prefetch
from .processor_daemon import run_processing_daemon

def _cfg(cfg, name, default):
    return getattr(cfg, name, default)

@dataclass
class ThreadPlan:
    bundle_concurrency: int
    bundle_threads: int
    prefetch_workers: int
    total_threads: int
    reserved_threads: int
    dump_cap: int | None
    fastp_cap: int | None
    kallisto_cap: int | None

def _plan_threads(cfg) -> ThreadPlan:
    total = int(_cfg(cfg, "threads", _cfg(cfg, "max_threads", os.cpu_count() or 8)))
    reserve = min(4, max(1, math.floor(total * 0.10)))
    avail = max(1, total - reserve)

    # how many bundles in parallel
    bundles = max(1, int(_cfg(cfg, "max_bundles", 2)))

    # threads per active tool inside the bundle
    forced_bt = _cfg(cfg, "bundle_threads", None)
    if forced_bt is not None:
        bundle_threads = max(1, int(forced_bt))
    else:
        bundle_threads = max(1, avail // bundles)

    # small, sensible caps (overridable in Config)
    dump_cap     = _cfg(cfg, "dump_cap", 1)         # fasterq-dump rarely benefits >1
    fastp_cap    = _cfg(cfg, "fastp_cap", 8)        # fastp scales but with diminishing returns
    kallisto_cap = _cfg(cfg, "kallisto_cap", 32)    # kallisto scales well up to ~32

    # prefetch is I/O-bound — give leftover lightweight workers (overrideable)
    default_pf = min(24, max(2, avail - bundles * bundle_threads))
    prefetch_workers = int(_cfg(cfg, "prefetch_workers", default_pf))

    return ThreadPlan(
        bundle_concurrency=bundles,
        bundle_threads=bundle_threads,
        prefetch_workers=prefetch_workers,
        total_threads=total,
        reserved_threads=reserve,
        dump_cap=None if dump_cap is None else int(dump_cap),
        fastp_cap=None if fastp_cap is None else int(fastp_cap),
        kallisto_cap=None if kallisto_cap is None else int(kallisto_cap),
    )

def run_download_and_process(
    dataset: "Dataset",
    *,
    cfg: "Config",
    cache_dir: Path,
    work_root: Path,
    temp_dir: Path,
    log_path: Path,
    cache_high_gb: int | None = None,
    cache_low_gb: int | None = None,
):
    cache_high_gb = cache_high_gb or _cfg(cfg, "cache_high_gb", 300)
    cache_low_gb  = cache_low_gb  or _cfg(cfg, "cache_low_gb", 250)

    plan = _plan_threads(cfg)
    log(
        (f"[orchestrator] threads={plan.total_threads} (reserve {plan.reserved_threads}) | "
         f"bundles={plan.bundle_concurrency} bundle_threads={plan.bundle_threads} | "
         f"caps(dump={plan.dump_cap},fastp={plan.fastp_cap},kallisto={plan.kallisto_cap}) | "
         f"prefetch_workers={plan.prefetch_workers} | "
         f"cache HIGH={cache_high_gb}GB LOW={cache_low_gb}GB"),
        log_path,
    )

    # Start processor first (consumer ready)
    t_proc = Thread(
        target=run_processing_daemon,
        kwargs=dict(
            dataset=dataset,
            cache_dir=cache_dir,
            work_root=work_root,
            temp_dir=temp_dir,
            log_path=log_path,
            max_bundles=plan.bundle_concurrency,
            bundle_threads=plan.bundle_threads,
            dump_cap=plan.dump_cap,
            fastp_cap=plan.fastp_cap,
            kallisto_cap=plan.kallisto_cap,
        ),
        daemon=True,
    )
    t_proc.start()

    # Prefetch (BP-prioritized list from dataset.to_do(); no-op in FASTQ mode)
    prefetch(
        dataset=dataset,
        cache_dir=cache_dir,
        log_path=log_path,
        max_workers=plan.prefetch_workers,
        cache_high_gb=cache_high_gb,
        cache_low_gb=cache_low_gb,
        overwrite=False,
    )

    t_proc.join()
    log("[orchestrator] all done", log_path)
