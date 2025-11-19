from __future__ import annotations

import os
import shutil
import time
from pathlib import Path
from typing import Dict

from concurrent.futures import ThreadPoolExecutor
from threading import BoundedSemaphore

from tqdm.auto import tqdm

from .utils import log, run_cmd, pad_desc
from .prefetcher import prefetch_one  # ← used in no-cache mode

POLL_SECS = 3.0


def _temp_fs_free_bytes(temp_dir: Path) -> int:
    st = os.statvfs(str(temp_dir))
    # use f_bavail to respect reserved blocks
    return st.f_bavail * st.f_frsize


def _rm(p: Path):
    if not p.exists():
        return
    if p.is_dir():
        shutil.rmtree(p, ignore_errors=True)
    else:
        try:
            p.unlink()
        except IsADirectoryError:
            shutil.rmtree(p, ignore_errors=True)


def _bundle_srr(
    sample: "Sample",
    *,
    cfg: "Config",
    cache_dir: Path,
    work_root: Path,
    temp_dir: Path,
    dump_threads: int = 1,
    fastp_threads: int = 4,
    kallisto_threads: int = 12,
) -> Dict[str, str]:
    """
    One SRR end-to-end: prefetch (optionally) → fasterq-dump → fastp → kallisto quant.

    In cache mode:
      • SRA is expected in shared cache (prefetch thread).
    In no-cache mode:
      • prefetch_one() is called here, writing a temporary .sra into sample.outdir.
    """
    srr = sample.id

    # SRR output directory: should already be something like <outdir>/<BioProject>/<SRR>
    outdir = sample.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    # temp dir for fasterq-dump
    tmp = (temp_dir / f"{srr}_tmp").resolve()
    tmp.mkdir(parents=True, exist_ok=True)

    # how many bootstraps? come from cfg.align settings, default 100
    bootstraps = int(getattr(cfg, "kallisto_bootstrap", 100))

    # Determine .sra location depending on cache mode
    if getattr(cfg, "no_cache", False):
        # no-cache mode: prefetch directly into run directory if needed
        if getattr(sample, "sra_path", None) is None or not Path(sample.sra_path).exists():
            pre_res = prefetch_one(
                sample,
                outdir,       # use SRR directory as "cache_dir" → SRR.sra lives here
                sample.log_path,
                overwrite=True,
            )
            if pre_res.get("status") != "prefetched":
                raise RuntimeError(f"Prefetch failed for {srr} in no-cache mode")
            sra = Path(pre_res["sra_path"])
            sample.sra_path = sra
        else:
            sra = Path(sample.sra_path)
    else:
        # cache mode: use shared cache location
        sra = sample.sra_path or (cache_dir / f"{srr}.sra")

    try:
        sample.status = "processing"

        # -------------------------------
        # fasterq-dump → raw FASTQs in outdir
        # -------------------------------
        run_cmd(
            [
                "fasterq-dump",
                (str(sra) if sra.exists() else srr),
                "--split-files",
                "--skip-technical",
                "--threads",
                str(dump_threads),
                "--outdir",
                str(outdir),
                "--temp",
                str(tmp),
            ],
            outdir,
            sample.log_path,
        )

        # detect FASTQs
        fq1 = next(outdir.glob(f"{srr}_1.fastq"), None)
        fq2 = next(outdir.glob(f"{srr}_2.fastq"), None)
        if fq1 is None and fq2 is None:
            # single-end fasterq sometimes yields just <srr>.fastq when not split
            single = next(outdir.glob(f"{srr}.fastq"), None)
            if single is None:
                raise FileNotFoundError(f"No FASTQ files for {srr} in {outdir}")
            fqs = [single]
        else:
            fqs = [p for p in (fq1, fq2) if p is not None]

        _rm(sra)

        # -------------------------------
        # fastp → trimmed FASTQs
        # -------------------------------
        trimmed = [outdir / (p.stem + ".trim.fastq") for p in fqs]
        cmd = ["fastp", "-w", str(fastp_threads)]
        if len(fqs) == 2:
            cmd += [
                "-i", str(fqs[0]),
                "-I", str(fqs[1]),
                "-o", str(trimmed[0]),
                "-O", str(trimmed[1]),
            ]
        else:
            cmd += [
                "-i", str(fqs[0]),
                "-o", str(trimmed[0]),
            ]
        run_cmd(cmd, outdir, sample.log_path)

        # CLEANUP STEP 1 — remove *raw* FASTQs immediately after fastp
        for p in fqs:
            if p is not None and p.exists():
                _rm(p)

        # -------------------------------
        # kallisto quant
        # -------------------------------
        # Global index path from config: reference.fa → reference.idx
        idx = cfg.reference_path.with_suffix(".idx")

        # Sanity check (optional but helpful)
        if not idx.exists():
            raise FileNotFoundError(f"kallisto index not found at {idx}")

        if len(trimmed) == 2:
            # paired-end
            run_cmd(
                [
                    "kallisto", "quant",
                    "-t", str(kallisto_threads),
                    "-i", str(idx),
                    "-o", str(outdir),
                    "--plaintext",
                    "-b", str(bootstraps),
                    str(trimmed[0]),
                    str(trimmed[1]),
                ],
                outdir,
                sample.log_path,
            )
        else:
            # single-end; generic fragment params
            run_cmd(
                [
                    "kallisto", "quant",
                    "-t", str(kallisto_threads),
                    "-i", str(idx),
                    "-o", str(outdir),
                    "--single",
                    "-l", "200",
                    "-s", "20",
                    "--plaintext",
                    "-b", str(bootstraps),
                    str(trimmed[0]),
                ],
                outdir,
                sample.log_path,
            )

        # CLEANUP STEP 2 — remove trimmed FASTQs now that kallisto is done
        for p in trimmed:
            if p.exists():
                _rm(p)

        # CLEANUP STEP 3 — temp dir for fasterq-dump
        _rm(tmp)

        # CLEANUP STEP 4 — in no-cache mode, remove the local .sra as well
        if getattr(cfg, "no_cache", False) and sra.exists() and sra.parent == outdir:
            _rm(sra)

        sample.status = "done"
        log(f"✅ [{srr}] bundle done", sample.log_path)
        return {"run_id": srr, "status": "done"}
    except Exception as e:
        sample.status = "failed"
        log(f"❌ [{srr}] failed: {e}", sample.log_path)
        return {"run_id": srr, "status": "failed"}


def _process_fastq(
    sample: "Sample",
    *,
    cfg: "Config",
    work_root: Path,
    fastp_threads: int = 4,
    kallisto_threads: int = 12,
) -> Dict[str, str]:
    """Placeholder for FASTQ-only mode (not yet implemented)."""
    sid = sample.id
    outdir = (work_root / sid).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        sample.status = "processing"
        # TODO: implement FASTQ-only flow mirroring _bundle_srr
        sample.status = "done"
        return {"run_id": sid, "status": "done"}
    except Exception as e:
        sample.status = "failed"
        log(f"❌ [{sid}] failed: {e}", sample.log_path)
        return {"run_id": sid, "status": "failed"}


def process(
    dataset: "Dataset",
    *,
    cfg: "Config",
    cache_dir: Path,
    work_root: Path,
    temp_dir: Path,
    log_path: Path,
    max_bundles: int = 2,
    poll_secs: float = POLL_SECS,
    tmp_low_water_gb: int = 0,
    dump_threads: int = 1,
    fastp_threads: int = 4,
    kallisto_threads: int = 12,
) -> Dict[str, str]:
    """
    Daemon that launches SRR / FASTQ bundles with concurrency control.

    In cache mode:
      • waits for samples to be 'prefetched' by the prefetch thread.

    In no-cache mode (cfg.no_cache=True):
      • does NOT wait for 'prefetched' status — _bundle_srr() will call prefetch_one() itself.
    """
    cache_dir = cache_dir.resolve()
    work_root = work_root.resolve()
    temp_dir = temp_dir.resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    temp_dir.mkdir(parents=True, exist_ok=True)

    results: Dict[str, str] = {}
    sem = BoundedSemaphore(max_bundles)
    pool = ThreadPoolExecutor(max_workers=max_bundles, thread_name_prefix="bundle")
    futures = set()

    total = len(dataset)

    # Count already finished samples
    initial_offset = sum(
        1 for s in dataset.samples
        if getattr(s, "status", None) in {"done", "failed", "skipped"}
    )

    pbar = tqdm(
        total=total,
        desc=pad_desc("Processing"),
        unit="sample",
        position=1,
        leave=True,
        mininterval=0,
    )

    if initial_offset:
        pbar.update(initial_offset)
        pbar.refresh()

    def maybe_launch(sample: "Sample"):
        # Only launch when ready
        if dataset.mode == "SRR":
            if getattr(cfg, "no_cache", False):
                # no-cache mode → don't wait for prefetch thread; _bundle_srr will prefetch
                pass
            else:
                # cache mode → require prefetch thread to have set status="prefetched"
                if getattr(sample, "status", None) != "prefetched":
                    return False
        elif dataset.mode == "FASTQ":
            if getattr(sample, "status", None) not in {"pending", "ready"}:
                return False

        if tmp_low_water_gb > 0 and _temp_fs_free_bytes(temp_dir) / (1024 ** 3) < tmp_low_water_gb:
            return False
        if not sem.acquire(blocking=False):
            return False

        if dataset.mode == "SRR":
            target = _bundle_srr
            kwargs = dict(
                cfg=cfg,
                cache_dir=cache_dir,
                work_root=work_root,
                temp_dir=temp_dir,
                dump_threads=dump_threads,
                fastp_threads=fastp_threads,
                kallisto_threads=kallisto_threads,
            )
        else:
            target = _process_fastq
            kwargs = dict(
                cfg=cfg,
                work_root=work_root,
                fastp_threads=fastp_threads,
                kallisto_threads=kallisto_threads,
            )

        fut = pool.submit(target, sample, **kwargs)
        futures.add(fut)
        return True

    try:
        while True:
            # collect finished bundles
            for f in list(futures):
                if f.done():
                    futures.remove(f)
                    sem.release()
                    res = f.result()
                    results[res["run_id"]] = res["status"]
                    pbar.update(1)

            # launch new work
            for s in dataset.to_do():
                maybe_launch(s)

            # exit condition
            if not futures and all(getattr(s, "status", None) in {"done", "failed", "skipped"} for s in dataset.to_do()):
                break

            time.sleep(poll_secs)
    finally:
        pool.shutdown(wait=True)
        pbar.close()

    ok = sum(1 for v in results.values() if v == "done")
    fail = sum(1 for v in results.values() if v == "failed")
    log(f"[process] finished: ok={ok} fail={fail}", log_path)
    return results
