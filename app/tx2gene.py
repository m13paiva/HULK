# tx2gene.py
import os
import sys
import pandas as pd
from pathlib import Path
from pytximport import tximport


from .utils import(
log,
log_err
)


def _gene_counts(
    abundance_files: list[Path],
    sample_names: list[str],
    tx2gene_path: Path,
    mode=None,
    ignore_tx_version: bool = False,
) -> pd.DataFrame:
    """
    Compute gene-level counts from kallisto outputs with pytximport (quiet mode).

    mode (user-facing) can be:
      - "raw_counts"          -> no transform (counts_from_abundance=None)
      - "scaled_tpm"          -> counts_from_abundance="scaled_tpm"
      - "length_scaled_tpm"   -> counts_from_abundance="length_scaled_tpm"
    """
    # ------------------------ normalize mode ------------------------
    mode_norm = (mode or "").strip().lower()

    if mode_norm in {"", "raw_counts", "no", "none"}:
        cfa = None
    elif mode_norm in {"scaled_tpm", "scaledtpm"}:
        cfa = "scaled_tpm"
    elif mode_norm in {"length_scaled_tpm", "lengthscaledtpm"}:
        cfa = "length_scaled_tpm"
    else:
        raise ValueError(
            "tximport_mode must be one of 'raw_counts', 'scaled_tpm', 'length_scaled_tpm' "
            f"(got {mode!r})"
        )

    prev_tqdm = os.environ.get("TQDM_DISABLE")
    os.environ["TQDM_DISABLE"] = "1"
    try:
        with open(os.devnull, "w") as devnull:
            sys_stdout, sys_stderr = sys.stdout, sys.stderr
            try:
                sys.stdout = devnull
                sys.stderr = devnull
                ds = tximport(
                    [str(p) for p in abundance_files],
                    data_type="kallisto",
                    transcript_gene_map=str(tx2gene_path),
                    counts_from_abundance=cfa,
                    output_type="xarray",
                    return_data=True,
                    existence_optional=False,
                    ignore_transcript_version=ignore_tx_version,
                )
            finally:
                sys.stdout, sys.stderr = sys_stdout, sys_stderr
    finally:
        if prev_tqdm is None:
            del os.environ["TQDM_DISABLE"]
        else:
            os.environ["TQDM_DISABLE"] = prev_tqdm

    df = ds["counts"].to_pandas()
    df.columns = sample_names
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Gene count merges
# ─────────────────────────────────────────────────────────────────────────────

def bp_gene_counts(bp,cfg) -> Path | None:
    """Write <bioproject_dir>/gene_counts.tsv (genes × samples). Logs missing SRAs but proceeds."""
    bioproject_dir = bp.path
    sra_ids = bp.get_sample_ids()
    tx2gene = cfg.tx2gene
    log_path = bp.log_path
    error_warnings = cfg.error_warnings
    tximport_opts = cfg.tximport_opts
    files, names = [], []
    for run_id in sra_ids:
        f = bioproject_dir / run_id / "abundance.tsv"
        if f.exists():
            files.append(f); names.append(run_id)
        else:
            log_err(error_warnings, log_path, f"[tximport] Missing file: {f}")

    if not files:
        log_err(error_warnings, log_path, f"[tximport] No abundance.tsv found in {bioproject_dir}, skipping gene table")
        return None
    _opts = {k: v for k, v in (tximport_opts or {}).items()
             if k in {"mode", "ignore_tx_version"} and v is not None}
    gdf = _gene_counts(files, names, tx2gene,**_opts)
    out = bioproject_dir / "gene_counts.tsv"
    gdf.to_csv(out, sep="\t")
    log(f"🧬 Gene counts (pytximport) written: {out}", log_path)
    return out

def global_gene_counts(output_root: Path, tx2gene: Path, log_path: Path, error_warnings: list[str], output_file: Path,tximport_opts=None) -> None:
    """
    Build a global gene count table across all BPs using pytximport directly.
    """
    files = list(output_root.glob("*/*/abundance.tsv"))
    if not files:
        log_err(error_warnings, log_path, f"[tximport] No abundance.tsv found under {output_root}")
        return
    names = [p.parent.name for p in files]
    _opts = {k: v for k, v in (tximport_opts or {}).items()
             if k in {"mode", "ignore_tx_version"} and v is not None}
    gdf = _gene_counts(files, names, tx2gene,**_opts)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_csv(output_file, sep="\t")
    log(f"🧬 Global gene counts (pytximport) written: {output_file}", log_path)
