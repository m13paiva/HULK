from __future__ import annotations

import subprocess
from pathlib import Path
from typing import List, Any, Optional

from .entities import Config, Dataset
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


def _build_r_args_global(cfg: Config, out_dir: Path, mode: str | None = None) -> List[str]:
    """
    Build argument list for a *global* run of the R post-processing script.
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
        "--prefix", "hulk",
        "--counts-from-abundance", counts_from_abundance,
    ]

    if mode:
        args.extend(["--mode", mode])

    if getattr(cfg, "tximport_ignore_tx_version", False):
        args.append("--ignore-tx-version")

    use_matrix = getattr(cfg, "expr_use_matrix", "vst") or "vst"
    args.extend(["--use-matrix", use_matrix])

    if not getattr(cfg, "drop_nonvarying_genes", True):
        args.append("--no-drop-nonvarying")

    var_thresh = getattr(cfg, "deseq2_var_threshold", 0.1)
    args.extend(["--var-threshold", str(var_thresh)])

    # Top N
    top_n = getattr(cfg, "top_n_vars", 500)
    args.extend(["--top-n", str(top_n)])

    # --- Target Genes: Join multiple files with commas ---
    target_files = getattr(cfg, "target_genes_files", [])
    if target_files:
        joined_targets = ",".join(str(tf) for tf in target_files)
        args.extend(["--target-genes", joined_targets])

    # Plot Modes
    deseq2_enabled = bool(getattr(cfg, "deseq2_vst_enabled", True))
    txi_only_mode = bool(getattr(cfg, "tximport_only_mode", False))
    plots_only_mode = bool(getattr(cfg, "plots_only_mode", False))

    if plots_only_mode:
        args.append("--plots-only")

    elif (not deseq2_enabled) or txi_only_mode:
        args.append("--tximport-only")
        # FORCE RE-IMPORT even in txi-only mode to capture new samples
        args.append("--force-txi")
        return args

    else:
        # FULL RUN: Force re-import to avoid stale cache from partial runs
        args.append("--force-txi")

    if deseq2_enabled:
        if getattr(cfg, "plot_pca", False):
            args.append("--pca")
        if getattr(cfg, "plot_heatmap", False):
            args.append("--heatmap")
        if getattr(cfg, "plot_var_heatmap", False):
            args.append("--var-heatmap")

    return args


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

    use_matrix = getattr(cfg, "expr_use_matrix", "vst")
    if use_matrix not in {"vst", "normalized"}:
        use_matrix = "vst"
    args.extend(["--use-matrix", use_matrix])

    if not getattr(cfg, "drop_nonvarying_genes", True):
        args.append("--no-drop-nonvarying")

    var_thresh = getattr(cfg, "deseq2_var_threshold", 0.1)
    args.extend(["--var-threshold", str(var_thresh)])

    # --- Target Genes: Join multiple files with commas (FIXED) ---
    target_files = getattr(cfg, "target_genes_files", [])
    if target_files:
        joined_targets = ",".join(str(tf) for tf in target_files)
        args.extend(["--target-genes", joined_targets])

    deseq2_enabled = bool(getattr(cfg, "deseq2_vst_enabled", True))
    plots_only_mode = bool(getattr(cfg, "plots_only_mode", False))
    tximport_only_mode = bool(getattr(cfg, "tximport_only_mode", False))

    if plots_only_mode:
        args.append("--plots-only")

    elif (not deseq2_enabled) or tximport_only_mode:
        args.append("--tximport-only")
        # FORCE RE-IMPORT (FIXED)
        args.append("--force-txi")
        return args, True

    else:
        # FULL RUN: FORCE RE-IMPORT (FIXED)
        args.append("--force-txi")

    if deseq2_enabled:
        if getattr(cfg, "plot_pca", False):
            args.append("--pca")
        if getattr(cfg, "plot_heatmap", False):
            args.append("--heatmap")
        if getattr(cfg, "plot_var_heatmap", False):
            args.append("--var-heatmap")

    return args, False


def run_postprocessing_bp(bp, cfg: Config, *, r_script: Path | None = None) -> Path | None:
    """
    Run the R post-processing script for a single BioProject.
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
    log(f"[{bp.id}] R post-processing (DESeq2/VST) completed.", log_path)
    return None


def run_postprocessing(dataset: Dataset, cfg: Config, *, r_script: Path | None = None, skip_bp: bool = False) -> None:
    """
    Run the R post-processing script:
    1. Once at the GLOBAL level (all samples).
    2. Once per BioProject (unless skip_bp=True).
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

    # ------------------------------------------------------------------
    # 1. Global Run
    # ------------------------------------------------------------------
    out_dir = cfg.shared
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        r_mode = getattr(dataset, "mode", None)
        args = _build_r_args_global(cfg, out_dir, mode=r_mode)
    except Exception as e:
        log_err(error_warnings, log_path, f"[post-processing] Failed to build global R command: {e}")
        return

    args[1] = str(script_path)

    try:
        log("[post-processing] Running Global R analysis...", log_path)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write("\n[post-processing] Global R command:\n")
            fh.write("  " + " ".join(args) + "\n\n")
            fh.flush()
            result = subprocess.run(args, stdout=fh, stderr=fh, text=True, check=False)

        if result.returncode != 0:
            log_err(error_warnings, log_path, f"[post-processing] Global R script exited with code {result.returncode}")
        else:
            log("[post-processing] Global analysis completed.", log_path)

    except Exception as e:
        log_err(error_warnings, log_path, f"[post-processing] Global execution failed: {e}")

    # ------------------------------------------------------------------
    # 2. Per-BioProject Run
    # ------------------------------------------------------------------
    if skip_bp:
        log("[post-processing] Skipping per-BioProject analysis (--no-bp-postprocessing).", log_path)
        return

    # Dataset.reconstruct_from_output populates dataset.bioprojects with a list of IDs (strings)
    if not dataset.bioprojects:
        return

    log(f"[post-processing] Starting analysis for {len(dataset.bioprojects)} BioProjects...", log_path)

    for bp_id in dataset.bioprojects:
        # Create a mock BioProject object that matches what run_postprocessing_bp expects
        # (It expects an object with .id, .path, and .log_path)
        bp_obj = SimpleNamespace(
            id=bp_id,
            path=cfg.outdir / bp_id,
            log_path=cfg.outdir / bp_id / "log.txt"
        )

        # Ensure BP log exists or at least the folder exists
        bp_obj.path.mkdir(parents=True, exist_ok=True)

        # Run specific BP logic
        run_postprocessing_bp(bp_obj, cfg, r_script=script_path)