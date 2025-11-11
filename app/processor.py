# processor_daemon.py (core bits)
from __future__ import annotations
import os, shutil, time
from pathlib import Path
from typing import Dict
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import BoundedSemaphore
from .utils import log, run_cmd

POLL_SECS = 3.0

def _temp_fs_free_bytes(temp_dir: Path) -> int:
    st = os.statvfs(str(temp_dir))
    return st.f_bavail * st.f_frsize

def _rm(p: Path):
    if not p.exists(): return
    if p.is_dir(): shutil.rmtree(p, ignore_errors=True)
    else:
        try: p.unlink()
        except IsADirectoryError: shutil.rmtree(p, ignore_errors=True)

def _bundle_srr(sample: "Sample", *, cache_dir: Path, work_root: Path, temp_dir: Path,
                log_path: Path, dump_threads=1, fastp_threads=4, kallisto_threads=12) -> Dict[str,str]:
    srr = sample.id
    sra = sample.sra_path or (cache_dir / f"{srr}.sra")
    outdir = (work_root / srr).resolve(); outdir.mkdir(parents=True, exist_ok=True)
    tmp = (temp_dir / f"{srr}_tmp").resolve(); tmp.mkdir(parents=True, exist_ok=True)
    try:
        sample.status = "processing"
        run_cmd([
            "fasterq-dump", (str(sra) if sra.exists() else srr),
            "--split-files", "--skip-technical",
            "--threads", str(dump_threads),
            "--outdir", str(outdir),
            "--temp", str(tmp),
        ], outdir, log_path)

        # detect FASTQs
        fq1 = next(outdir.glob(f"{srr}_1.fastq"), None)
        fq2 = next(outdir.glob(f"{srr}_2.fastq"), None)
        fqs = [p for p in [fq1, fq2] if p] or [next(outdir.glob(f"{srr}.fastq"))]

        # fastp
        trimmed = [outdir / (p.stem + ".trim.fastq") for p in fqs]
        cmd = ["fastp","-w",str(fastp_threads)]
        if len(fqs) == 2:
            cmd += ["-i",str(fqs[0]),"-I",str(fqs[1]),"-o",str(trimmed[0]),"-O",str(trimmed[1])]
        else:
            cmd += ["-i",str(fqs[0]),"-o",str(trimmed[0])]
        run_cmd(cmd, outdir, log_path)

        # kallisto (adapt index source to your Config/metadata)
        if len(trimmed) == 2:
            run_cmd(["kallisto","quant","-t",str(kallisto_threads),"-o",str(outdir),"--plaintext",
                     sample.metadata["kallisto_index"], str(trimmed[0]), str(trimmed[1])], outdir, log_path)
        else:
            run_cmd(["kallisto","quant","-t",str(kallisto_threads),"-o",str(outdir),
                     "--single","-l","200","-s","20","--plaintext",
                     sample.metadata["kallisto_index"], str(trimmed[0])], outdir, log_path)

        # cleanup
        for p in outdir.glob("*.fastq"): _rm(p)
        _rm(tmp);  _rm(sra)
        sample.status = "done"
        log(f"✅ [{srr}] bundle done", log_path)
        return {"run_id": srr, "status": "done"}
    except Exception as e:
        sample.status = "failed"
        log(f"❌ [{srr}] failed: {e}", log_path)
        return {"run_id": srr, "status": "failed"}

def _process_fastq(sample: "Sample", *, work_root: Path, log_path: Path,
                   fastp_threads=4, kallisto_threads=12) -> Dict[str,str]:
    # sketch: adapt to your existing FASTQ-only pipeline
    sid = sample.id
    outdir = (work_root / sid).resolve(); outdir.mkdir(parents=True, exist_ok=True)
    try:
        sample.status = "processing"
        # fastp on sample.fastq_paths -> trimmed
        # kallisto quant -> outdir
        sample.status = "done"
        return {"run_id": sid, "status": "done"}
    except Exception as e:
        sample.status = "failed"
        log(f"❌ [{sid}] failed: {e}", log_path)
        return {"run_id": sid, "status": "failed"}

def run_processing_daemon(dataset: "Dataset", *, cache_dir: Path, work_root: Path, temp_dir: Path,
                          log_path: Path, max_bundles=2, poll_secs=POLL_SECS,
                          tmp_low_water_gb=0, dump_threads=1, fastp_threads=4, kallisto_threads=12) -> Dict[str,str]:
    cache_dir = cache_dir.resolve(); work_root = work_root.resolve(); temp_dir = temp_dir.resolve()
    work_root.mkdir(parents=True, exist_ok=True); temp_dir.mkdir(parents=True, exist_ok=True)

    results: Dict[str,str] = {}
    sem = BoundedSemaphore(max_bundles)
    pool = ThreadPoolExecutor(max_workers=max_bundles, thread_name_prefix="bundle")
    futures = set()

    def maybe_launch(sample: "Sample"):
        # SRR: only launch prefetched; FASTQ: launch pending
        if dataset.mode == "SRR" and sample.status != "prefetched":
            return False
        if dataset.mode == "FASTQ" and sample.status not in {"pending","ready"}:
            return False
        if tmp_low_water_gb > 0 and _temp_fs_free_bytes(temp_dir)/(1024**3) < tmp_low_water_gb:
            return False
        if not sem.acquire(blocking=False):
            return False
        fut = pool.submit(
            _bundle_srr if dataset.mode == "SRR" else _process_fastq,
            sample,
            **({"cache_dir":cache_dir,"work_root":work_root,"temp_dir":temp_dir,"log_path":log_path,
                "dump_threads":dump_threads,"fastp_threads":fastp_threads,"kallisto_threads":kallisto_threads}
               if dataset.mode=="SRR"
               else {"work_root":work_root,"log_path":log_path,
                     "fastp_threads":fastp_threads,"kallisto_threads":kallisto_threads})
        )
        futures.add(fut); return True

    try:
        while True:
            # collect done
            for f in list(futures):
                if f.done():
                    futures.remove(f); sem.release()
                    res = f.result(); results[res["run_id"]] = res["status"]

            # launch in priority order
            for s in dataset.to_do():
                maybe_launch(s)

            if not futures and all(s.status in {"done","failed","skipped"} for s in dataset.to_do()):
                break
            time.sleep(poll_secs)
    finally:
        pool.shutdown(wait=True)

    ok = sum(1 for v in results.values() if v == "done")
    fail = sum(1 for v in results.values() if v == "failed")
    log(f"[process] finished: ok={ok} fail={fail}", log_path)
    return results
