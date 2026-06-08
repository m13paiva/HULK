from __future__ import annotations
import shutil
import subprocess
import random
import time
from pathlib import Path
from typing import Dict
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm.auto import tqdm
from .utils import log, run_cmd, pad_desc

try:
    from .cache_manager import CacheGate
except Exception:
    class CacheGate:
        """Fallback class mirroring the basic methods if cache gating fails to import."""

        def __init__(self, _cache_dir: Path, _high: int, _low: int, _poll: float) -> None:
            pass

        def wait_for_window(self, _logger, _run_id: str):
            return

PREFETCH_MAX_DEFAULT = "200G"
POLL_SECS_DEFAULT = 3.0
MAX_WORKERS_DEFAULT = 16


def _normalize_sra_layout(cache_dir: Path, run_id: str) -> Path:
    """
    Standardizes SRA output locations from nested structures created by `prefetch`.

    Args:
        cache_dir (Path): The designated target cache path.
        run_id (str): The accession identifier determining expected file naming.

    Returns:
        Path: The absolute path to the properly structured SRA file.
    """
    cache_dir = cache_dir.expanduser().resolve()
    sra_file = cache_dir / f"{run_id}.sra"
    alt_dir = cache_dir / run_id
    alt_file = alt_dir / f"{run_id}.sra"

    if (not sra_file.exists()) and alt_file.exists():
        alt_file.replace(sra_file)
        shutil.rmtree(alt_dir, ignore_errors=True)
    return sra_file


def prefetch_one(sample: "Sample", cache_dir: Path, *,
                 overwrite=False, retries=3, prefetch_max=PREFETCH_MAX_DEFAULT,
                 cache_high_gb=300, cache_low_gb=250, poll_secs=POLL_SECS_DEFAULT,
                 mode: str = "cache") -> Dict[str, str]:
    """
    Retrieves individual SRA archives utilizing SRA Toolkit logic while respecting capacity gates.

    Args:
        sample (Sample): The sample representing the target accession.
        cache_dir (Path): Base destination layer.
        overwrite (bool, optional): Overwrite existing valid cached files. Defaults to False.
        retries (int, optional): Max connection retries. Defaults to 3.
        prefetch_max (str, optional): Max download boundary string recognized by toolkit. Defaults to '200G'.
        cache_high_gb (int, optional): High watermark block trigger. Defaults to 300.
        cache_low_gb (int, optional): Low watermark unblock trigger. Defaults to 250.
        poll_secs (float, optional): Blocked polling interval. Defaults to 3.0.
        mode (str, optional): Cache mode operation flag. Defaults to "cache".

    Returns:
        Dict[str, str]: Dictionary mapping the target accession to completion status and path.
    """
    cache_dir = cache_dir.expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    run_id = sample.id

    if mode == "cache":
        CacheGate(
            cache_dir,
            cache_high_gb * (1024 ** 3),
            cache_low_gb * (1024 ** 3),
            poll_secs,
        ).wait_for_window(lambda m: log(m, sample.log_path), run_id)

    if mode == "cache":
        sra_file = cache_dir / f"{run_id}.sra"
    else:
        run_dir = sample.outdir
        run_dir.mkdir(parents=True, exist_ok=True)
        sra_file = run_dir / f"{run_id}.sra"

    if mode == "cache" and sra_file.exists() and not overwrite:
        sample.sra_path = sra_file
        sample.status = "prefetched"
        log(f"[{run_id}] SKIP prefetch: already cached.", sample.log_path)
        return {"run_id": run_id, "status": "prefetched", "sra_path": str(sra_file)}

    if overwrite:
        try:
            if sra_file.exists(): sra_file.unlink()
            alt_dir = cache_dir / run_id
            if alt_dir.exists(): shutil.rmtree(alt_dir, ignore_errors=True)
        except Exception:
            pass

    x_opt = ["-X", prefetch_max] if (mode == "cache" and prefetch_max) else []

    strategies = [
        ["prefetch", "--force", "all"] + x_opt +
        ["--type", "sra", "--output-file", str(sra_file), run_id]
    ]

    attempt, last_err = 0, None
    while attempt < max(1, retries):
        attempt += 1
        for i, cmd in enumerate(strategies, 1):
            try:
                log(f"[{run_id}] prefetch attempt {attempt}.{i} …", sample.log_path)
                run_cmd(cmd, cache_dir, sample.log_path)
                sra_file = _normalize_sra_layout(cache_dir, run_id)
                if not sra_file.exists():
                    raise FileNotFoundError(sra_file)
                sample.sra_path = sra_file
                sample.status = "prefetched"
                log(f"✅ Prefetched [{run_id}] → {sra_file.name}", sample.log_path)
                return {"run_id": run_id, "status": "prefetched", "sra_path": str(sra_file)}
            except subprocess.CalledProcessError as e:
                last_err = f"prefetch exit {e.returncode}"
                log(f"[{run_id}] strategy {i} failed: {last_err}", sample.log_path)
            except Exception as e:
                last_err = str(e)
                log(f"[{run_id}] strategy {i} error: {last_err}", sample.log_path)

        if attempt < retries:
            delay = min(60, 2 ** attempt + random.random())
            log(f"[{run_id}] retrying in {delay:.1f}s …", sample.log_path)
            time.sleep(delay)

    sample.status = "failed"
    log(f"❌ Prefetch FAILED [{run_id}]: {last_err or 'unknown error'}", sample.log_path)
    return {"run_id": run_id, "status": "failed", "sra_path": str(sra_file)}


def prefetch(dataset: "Dataset", cache_dir: Path, log_path: Path, *,
             max_workers=MAX_WORKERS_DEFAULT, retries=3, prefetch_max=PREFETCH_MAX_DEFAULT,
             cache_high_gb=300, cache_low_gb=250, poll_secs=POLL_SECS_DEFAULT,
             overwrite=False, mode="cache") -> Dict[str, Dict[str, str]]:
    """
    Manages concurrent retrieval scheduling mapped against dataset state.

    Args:
        dataset (Dataset): System state matrix containing target resolution objects.
        cache_dir (Path): Output caching layer.
        log_path (Path): Logging routing destination.
        max_workers (int, optional): Parallel download threads limit. Defaults to MAX_WORKERS_DEFAULT.
        retries (int, optional): Internal failure retries per worker block. Defaults to 3.
        prefetch_max (str, optional): Connection max parameter. Defaults to PREFETCH_MAX_DEFAULT.
        cache_high_gb (int, optional): Pause trigger threshold limit. Defaults to 300.
        cache_low_gb (int, optional): Unpause trigger threshold limit. Defaults to 250.
        poll_secs (float, optional): Block cycle interval checking. Defaults to 3.0.
        overwrite (bool, optional): Allow block overwrites if True. Defaults to False.
        mode (str, optional): Target operation structure. Defaults to "cache".

    Returns:
        Dict[str, Dict[str, str]]: Mapped aggregation of retrieval statuses across all queued samples.
    """
    results: Dict[str, Dict[str, str]] = {}

    cache_dir = cache_dir.expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)

    all_srrs = [
        s for bp in dataset.bioprojects
        for s in bp.samples
        if s.type == "SRR"
    ]
    total = len(all_srrs)
    if total == 0:
        return results

    candidates = [
        s for s in dataset.to_do()
        if s.type == "SRR" and s.status != "prefetched"
    ]

    already_done = sum(1 for s in all_srrs if s.status == "done")
    already_prefetched = sum(1 for s in all_srrs if s.status == "prefetched")
    initial_offset = already_done + already_prefetched

    log(
        f"[prefetch] start: {len(candidates)} SRR to fetch "
        f"({already_prefetched} already prefetched, {already_done} done) | total={total}",
        log_path,
    )

    def task(s: "Sample"):
        return prefetch_one(
            s, cache_dir,
            overwrite=overwrite, retries=retries,
            prefetch_max=prefetch_max,
            cache_high_gb=cache_high_gb, cache_low_gb=cache_low_gb,
            poll_secs=poll_secs,
            mode=mode,
        )

    with ThreadPoolExecutor(max_workers=max_workers,
                            thread_name_prefix="prefetch") as ex:
        futs = [ex.submit(task, s) for s in candidates]

        with tqdm(
                total=total,
                desc=pad_desc("Prefetch"),
                unit="Sample",
                position=0,
                leave=True,
                mininterval=0,
        ) as pbar:

            if initial_offset:
                pbar.update(initial_offset)
                pbar.refresh()

            for fut in as_completed(futs):
                res = fut.result()
                results[res["run_id"]] = res

                ok = sum(1 for r in results.values() if r["status"] == "prefetched")
                fail = sum(1 for r in results.values() if r["status"] == "failed")
                pbar.set_postfix(ok=ok, fail=fail)

                pbar.update(1)
                pbar.refresh()

    log(
        f"[prefetch] done. ok={sum(1 for r in results.values() if r['status'] == 'prefetched')}, "
        f"fail={sum(1 for r in results.values() if r['status'] == 'failed')}",
        log_path,
    )

    return results