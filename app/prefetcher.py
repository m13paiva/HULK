# prefetcher.py
from __future__ import annotations
import shutil
import subprocess
import random
import time
from pathlib import Path
from typing import Iterable, Dict, Optional, List, Any
from concurrent.futures import ThreadPoolExecutor, as_completed

from .utils import log, run_cmd
from .cache_manager import CacheGate  # your polling gate

PREFETCH_MAX_DEFAULT = "200G"
POLL_SECS_DEFAULT    = 3.0
MAX_WORKERS_DEFAULT  = 16

def _normalize_sra_layout(cache_dir: Path, run_id: str) -> Path:
    cache_dir = cache_dir.expanduser().resolve()
    sra_file = cache_dir / f"{run_id}.sra"
    alt_dir  = cache_dir / run_id
    alt_file = alt_dir / f"{run_id}.sra"
    if (not sra_file.exists()) and alt_file.exists():
        alt_file.replace(sra_file)
        shutil.rmtree(alt_dir, ignore_errors=True)
    return sra_file

def prefetch_one(sample: "Sample", cache_dir: Path, log_path: Path, *,
                 overwrite=False, retries=3, prefetch_max=PREFETCH_MAX_DEFAULT,
                 cache_high_gb=300, cache_low_gb=250, poll_secs=POLL_SECS_DEFAULT) -> Dict[str,str]:

    cache_dir = cache_dir.expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    run_id = sample.id

    # pause/resume by watermarks
    CacheGate(cache_dir, cache_high_gb*(1024**3), cache_low_gb*(1024**3), poll_secs)\
        .wait_for_window(lambda m: log(m, log_path), run_id)

    sra_file = cache_dir / f"{run_id}.sra"
    if sra_file.exists() and not overwrite:
        sample.sra_path = sra_file
        sample.status = "prefetched"
        log(f"[{run_id}] SKIP prefetch: already cached.", log_path)
        return {"run_id": run_id, "status": "prefetched", "sra_path": str(sra_file)}

    if overwrite:
        try:
            if sra_file.exists(): sra_file.unlink()
            alt_dir = cache_dir / run_id
            if alt_dir.exists(): shutil.rmtree(alt_dir, ignore_errors=True)
        except Exception:
            pass

    strategies = [
        ["prefetch","--force","all","-X",prefetch_max,"--type","sra","--output-file",str(sra_file),run_id],
        ["prefetch","--force","all","-X",prefetch_max,"--output-directory",str(cache_dir),run_id],
    ]

    attempt, last_err = 0, None
    while attempt < max(1, retries):
        attempt += 1
        for i, cmd in enumerate(strategies, 1):
            try:
                log(f"[{run_id}] prefetch attempt {attempt}.{i} …", log_path)
                run_cmd(cmd, cache_dir, log_path)
                sra_file = _normalize_sra_layout(cache_dir, run_id)
                if not sra_file.exists():
                    raise FileNotFoundError(sra_file)
                sample.sra_path = sra_file
                sample.status = "prefetched"
                log(f"✅ Prefetched [{run_id}] → {sra_file.name}", log_path)
                return {"run_id": run_id, "status": "prefetched", "sra_path": str(sra_file)}
            except subprocess.CalledProcessError as e:
                last_err = f"prefetch exit {e.returncode}"
                log(f"[{run_id}] strategy {i} failed: {last_err}", log_path)
            except Exception as e:
                last_err = str(e)
                log(f"[{run_id}] strategy {i} error: {last_err}", log_path)

        if attempt < retries:
            delay = min(60, 2**attempt + random.random())
            log(f"[{run_id}] retrying in {delay:.1f}s …", log_path)
            time.sleep(delay)

    sample.status = "failed"
    log(f"❌ Prefetch FAILED [{run_id}]: {last_err or 'unknown error'}", log_path)
    return {"run_id": run_id, "status": "failed", "sra_path": str(sra_file)}

def prefetch(dataset: "Dataset", cache_dir: Path, log_path: Path, *,
             max_workers=MAX_WORKERS_DEFAULT, retries=3, prefetch_max=PREFETCH_MAX_DEFAULT,
             cache_high_gb=300, cache_low_gb=250, poll_secs=POLL_SECS_DEFAULT,
             overwrite=False) -> Dict[str, Dict[str,str]]:
    """
    Accepts Dataset. In SRR mode, prefetches pending SRR samples (BP-prioritized).
    In FASTQ mode, nothing to prefetch.
    """
    results: Dict[str, Dict[str,str]] = {}

    cache_dir = cache_dir.expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)

    candidates = [s for s in dataset.to_do() if s.type == "SRR"]
    total = len(candidates)
    log(f"[prefetch] start: {total} SRR | workers={max_workers} | cache={cache_dir}", log_path)

    def task(s: "Sample"):
        return prefetch_one(s, cache_dir, log_path,
                            overwrite=overwrite, retries=retries, prefetch_max=prefetch_max,
                            cache_high_gb=cache_high_gb, cache_low_gb=cache_low_gb, poll_secs=poll_secs)

    with ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="prefetch") as ex:
        futs = [ex.submit(task, s) for s in candidates]
        done = 0
        for fut in as_completed(futs):
            res = fut.result()
            results[res["run_id"]] = res
            done += 1
            ok = sum(1 for r in results.values() if r["status"] == "prefetched")
            fail = sum(1 for r in results.values() if r["status"] == "failed")
            log(f"[prefetch] progress {done}/{total} | ok={ok} fail={fail}", log_path)

    log(f"[prefetch] done. ok={sum(1 for r in results.values() if r['status']=='prefetched')}, "
        f"fail={sum(1 for r in results.values() if r['status']=='failed')}", log_path)
    return results
