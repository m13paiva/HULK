# sra.py
import shutil
import subprocess
from pathlib import Path

from .utils import (
    log,
    run_cmd,
    is_sra_done
)


def download(run, threads: int, outdir: Path,
                     log_path: Path, overwrite: bool) -> None:
    """
        Download one SRR's FastQ file(s):
          prefetch → fasterq-dump
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
    prefetch_try1 = ["prefetch", "--force", "all", "-X", "200G", "--type", "sra", "--output-file", str(sra_file),
                     run_id]
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