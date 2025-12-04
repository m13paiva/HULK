# post_processing.py

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import List, Any

from .entities import Config
from .utils import log, log_err

R_SCRIPT_PATH = Path(__file__).with_name("post_processing.R")


def _map_counts_from_abundance(mode: str | None) -> str:
    """
    Map Config.tximport_mode (python side) to the R script's
    --counts-from-abundance argument.
    """
    m = (mode or "").strip().lower()
    if m in {"", "raw_counts", "no", "none"}:
        return "no"
    if m in {"scaled_tpm", "scaledtpm"}:
        return "scaledTPM"
    if m in {"length_scaled_tpm", "lengthscaledtpm"}:
        return "lengthScaledTPM"
    if m in {"dtu_scaled_tpm", "dtuscaledtpm"}:
        return "dtuScaledTPM"
    raise ValueError(
        "Invalid tximport_mode for R post-processing: "
        f"{mode!r}. Expected one of: raw_counts, scaled_tpm, "
        "length_scaled_tpm, dtu_scaled_tpm"
    )


def _build_r_args_for_bp(bp_id: str, cfg: Config, out_dir: Path) -> tuple[List[str], bool]:
    """
    Build argument list for a *per-BioProject* run of the R script.

    Returns:
      (args, tximport_only_mode_flag)
    """
    if cfg.tx2gene is None:
        raise ValueError("tx2gene is required for post-processing but cfg.tx2gene is None.")

    counts_from_abundance = _map_counts_from_abundance(getattr(cfg, "tximport_mode", None))

    args: List[str] = [
        "Rscript",
        str(R_SCRIPT_PATH),
        "--search-dir", str(cfg.outdir),
        "--tx2gene", str(cfg.tx2gene),
        "--out-dir", str(out_dir),
        "--prefix", bp_id,
        "--bioproject", bp_id,
        "--counts-from-abundance", counts_from_abundance,
    ]

    if getattr(cfg, "tximport_ignore_tx_version", False):
        args.append("--ignore-tx-version")

    # Matrix for Seidr exports: "vst" | "normalized"
    use_matrix = getattr(cfg, "expr_use_matrix", "vst")
    if use_matrix not in {"vst", "normalized"}:
        use_matrix = "vst"
    args.extend(["--use-matrix", use_matrix])

    # Drop non-varying genes?  (R flag is inverse: --no-drop-nonvarying)
    if not getattr(cfg, "drop_nonvarying_genes", True):
        args.append("--no-drop-nonvarying")

    # Optional target gene list (also applies in BP mode)
    if getattr(cfg, "target_genes_file", None):
        args.extend(["--target-genes", str(cfg.target_genes_file)])

    # Decide mode:
    # - plots-only       -> use existing VST, just plots
    # - tximport-only    -> gene counts only
    # - full DESeq2+VST  -> full treatment, incl. plots if requested
    deseq2_enabled = bool(getattr(cfg, "deseq2_vst_enabled", True))
    plots_only_mode = bool(getattr(cfg, "plots_only_mode", False))
    tximport_only_mode = bool(getattr(cfg, "tximport_only_mode", False))

    if plots_only_mode:
        args.append("--plots-only")
        # plots-only implicitly needs DESeq2/VST to have been run before.

    elif (not deseq2_enabled) or tximport_only_mode:
        # Pure tximport-only, this is where the R script writes tximport_counts.tsv
        args.append("--tximport-only")
        # Do NOT pass any plot flags; R quits right after writing counts.
        return args, True

    # Full DESeq2+VST pipeline (BP-specific "full treatment")
    if deseq2_enabled:
        if getattr(cfg, "plot_pca", False):
            args.append("--pca")
        if getattr(cfg, "plot_heatmap", False):
            args.append("--heatmap")
        if getattr(cfg, "plot_var_heatmap", False):
            # R script will automatically disable var-heatmap in --bioproject mode
            args.append("--var-heatmap")

    return args, False


def run_postprocessing_bp(bp, cfg: Config, *, r_script: Path | None = None) -> Path | None:
    """
    Run the R post-processing script for a single BioProject.

    Behaviour:
      - If DESeq2 is enabled (and not tximport-only):
          → full DESeq2 + VST + Seidr-style exports
          → BP-level PCA / expression heatmap (if plot flags enabled)
          → returns None (no gene_counts.tsv from R in this mode)
      - If DESeq2 is disabled OR tximport-only mode is set:
          → tximport-only; R writes <prefix>.tximport_counts.tsv
          → we move that to <BP>/gene_counts.tsv and return that Path.
    """
    log_path = bp.log_path
    errors = cfg.error_warnings

    if cfg.tx2gene is None:
        log_err(
            errors,
            log_path,
            f"[{bp.id}] [post-processing] tx2gene not provided; skipping R post-processing.",
        )
        return None

    script_path = Path(r_script) if r_script is not None else R_SCRIPT_PATH
    if not script_path.exists():
        log_err(
            errors,
            log_path,
            f"[{bp.id}] [post-processing] R script not found at {script_path}.",
        )
        return None

    out_dir = bp.path
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        args, is_txi_only = _build_r_args_for_bp(bp.id, cfg, out_dir)
    except Exception as e:
        log_err(
            errors,
            log_path,
            f"[{bp.id}] [post-processing] Failed to build R command: {e}",
        )
        return None

    # Ensure correct script path is set
    args[1] = str(script_path)

    try:
        mode_str = "tximport-only" if is_txi_only else "full DESeq2/VST"
        log(f"[{bp.id}] [post-processing] Running R script ({mode_str})…", log_path)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(f"\n[{bp.id}] [post-processing] R command:\n")
            fh.write("  " + " ".join(args) + "\n\n")
            fh.flush()
            result = subprocess.run(
                args,
                stdout=fh,
                stderr=fh,
                text=True,
                check=False,
            )
        if result.returncode != 0:
            log_err(
                errors,
                log_path,
                f"[{bp.id}] [post-processing] R script exited with code {result.returncode}",
            )
            return None
    except FileNotFoundError as e:
        log_err(
            errors,
            log_path,
            f"[{bp.id}] [post-processing] Failed to execute Rscript: {e}",
        )
        return None
    except Exception as e:
        log_err(
            errors,
            log_path,
            f"[{bp.id}] [post-processing] Unexpected error while running R script: {e}",
        )
        return None

    # If we were in tximport-only mode, the R script wrote tximport counts
    if is_txi_only:
        counts_file = out_dir / f"{bp.id}.tximport_counts.tsv"
        if not counts_file.exists():
            log_err(
                errors,
                log_path,
                f"[{bp.id}] [post-processing] Expected counts file not found: {counts_file}",
            )
            return None

        final_path = bp.path / "gene_counts.tsv"
        try:
            counts_file.replace(final_path)
        except Exception:
            import shutil
            try:
                shutil.copy2(counts_file, final_path)
            except Exception as e:
                log_err(
                    errors,
                    log_path,
                    f"[{bp.id}] [post-processing] Failed to move counts file: {e}",
                )
                return None

        log(f"[{bp.id}] Gene counts written: {final_path}", log_path)
        return final_path

    # Full DESeq2/VST mode: no gene_counts.tsv here, but all the juicy stuff
    # (VST, normalized counts, PCA/heatmaps, Seidr exports) is in out_dir.
    log(f"[{bp.id}] R post-processing (DESeq2/VST) completed.", log_path)
    return None

def _build_r_args_global(cfg: Config, out_dir: Path, mode: str | None = None) -> List[str]:

    """
    Build argument list for a *global* run of the R post-processing
    script (no --bioproject restriction).
    """
    if cfg.tx2gene is None:
        raise ValueError("tx2gene is required for post-processing but cfg.tx2gene is None.")

    counts_from_abundance = _map_counts_from_abundance(getattr(cfg, "tximport_mode", None))

    args: List[str] = [
        "Rscript",
        str(R_SCRIPT_PATH),
        "--search-dir", str(cfg.outdir),
        "--tx2gene", str(cfg.tx2gene),
        "--out-dir", str(out_dir),
        "--prefix", "hulk",  # arbitrary; R script default is 'seidr_input'
        "--counts-from-abundance", counts_from_abundance,
    ]

    # Mode hint for the R script (SRR vs FASTQ)
    if mode:
        args.extend(["--mode", mode])

    # tximport ignore version
    if getattr(cfg, "tximport_ignore_tx_version", False):
        args.append("--ignore-tx-version")

    # Which expression matrix to export for Seidr: vst | normalized
    use_matrix = getattr(cfg, "expr_use_matrix", "vst")
    if use_matrix not in {"vst", "normalized"}:
        use_matrix = "vst"
    args.extend(["--use-matrix", use_matrix])

    # Should we drop non-varying genes in the Seidr exports?
    # R flag is *inverse* semantics: --no-drop-nonvarying
    if not getattr(cfg, "drop_nonvarying_genes", True):
        args.append("--no-drop-nonvarying")

    # Optional target genes file (used for targeted plots, var heatmap, etc.)
    if getattr(cfg, "target_genes_file", None):
        args.extend(["--target-genes", str(cfg.target_genes_file)])

    # Determine high-level mode: tximport-only vs full DESeq2/plots.
    # NOTE: plotting only occurs when DESeq2 is enabled.
    deseq2_enabled = bool(getattr(cfg, "deseq2_vst_enabled", True))
    txi_only_mode = bool(getattr(cfg, "tximport_only_mode", False))
    plots_only_mode = bool(getattr(cfg, "plots_only_mode", False))

    if plots_only_mode:
        # Use existing VST file only (no tximport/DESeq2).
        args.append("--plots-only")
        # In this mode we still honour plot flags below.
    elif (not deseq2_enabled) or txi_only_mode:
        # Only run tximport and write counts table; no DESeq2, VST, or Seidr exports.
        args.append("--tximport-only")
        # We intentionally DO NOT pass any plot flags in this mode.
        return args

    # At this point: DESeq2/VST is enabled and we're not in tximport-only mode.

    # Plot flags (global-level)
    if deseq2_enabled:
        if getattr(cfg, "plot_pca", False):
            args.append("--pca")
        if getattr(cfg, "plot_heatmap", False):
            args.append("--heatmap")
        if getattr(cfg, "plot_var_heatmap", False):
            args.append("--var-heatmap")


    return args


def run_postprocessing(dataset: Dataset, cfg: Config, *, r_script: Path | None = None) -> None:
    """
    Run the R post-processing script once at the global level, using
    parameters stored in Config.

    - SEARCH_DIR is cfg.outdir (root that contains <BP>/<SRR>/abundance.tsv)
    - tx2gene comes from cfg.tx2gene
    - Output is written under <outdir>/shared/post_processing
    - tximport options come from cfg.tximport_mode and cfg.tximport_ignore_tx_version
    - DESeq2 + VST + Seidr exports/plots are controlled by:
        cfg.deseq2_vst_enabled
        cfg.expr_use_matrix
        cfg.drop_nonvarying_genes
        cfg.plot_pca
        cfg.plot_heatmap
        cfg.plot_var_heatmap
        cfg.plots_only_mode
        cfg.tximport_only_mode

    Plotting only happens if cfg.deseq2_vst_enabled is True.
    """
    log_path = cfg.log
    error_warnings: List[str] = cfg.error_warnings

    if cfg.tx2gene is None:
        log_err(
            error_warnings,
            log_path,
            "[post-processing] tx2gene not provided; skipping R post-processing.",
        )
        return

    script_path = Path(r_script) if r_script is not None else R_SCRIPT_PATH
    if not script_path.exists():
        log_err(
            error_warnings,
            log_path,
            f"[post-processing] R script not found at {script_path}. "
            "Place your R script there or pass an explicit path.",
        )
        return

    # Where to put all outputs from the R script
    out_dir = cfg.shared
    out_dir.mkdir(parents=True, exist_ok=True)

    # Build command-line arguments
    try:
        r_mode = getattr(dataset, "mode", None)  # "SRR" or "FASTQ"
        args = _build_r_args_global(cfg, out_dir, mode=r_mode)

    except Exception as e:
        log_err(
            error_warnings,
            log_path,
            f"[post-processing] Failed to build R command: {e}",
        )
        return

    # Replace script path in args (in case user overrode r_script argument)
    args[1] = str(script_path)

    # Run Rscript, streaming stdout/stderr into the main HULK log
    try:
        log("[post-processing] Running R script for global tximport/DESeq2/VST/plots…", log_path)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write("\n[post-processing] R command:\n")
            fh.write("  " + " ".join(args) + "\n\n")
            fh.flush()
            result = subprocess.run(
                args,
                stdout=fh,
                stderr=fh,
                text=True,
                check=False,
            )
        if result.returncode != 0:
            log_err(
                error_warnings,
                log_path,
                f"[post-processing] R script exited with code {result.returncode}",
            )
        else:
            log("[post-processing] R post-processing completed successfully.", log_path)
    except FileNotFoundError as e:
        # Typically Rscript not found
        log_err(
            error_warnings,
            log_path,
            f"[post-processing] Failed to execute Rscript: {e}",
        )
    except Exception as e:
        log_err(
            error_warnings,
            log_path,
            f"[post-processing] Unexpected error while running R script: {e}",
        )
