from __future__ import annotations
import subprocess
from pathlib import Path
from types import SimpleNamespace
from typing import List, Any, Optional, Set
from concurrent.futures import ThreadPoolExecutor, as_completed

from .entities import Config, Dataset
from .utils import log, log_err

R_SCRIPT_PATH = Path(__file__).with_name("post_processing.R")


def _map_counts_from_abundance(mode: str | None) -> str:
    if mode == "length_scaled_tpm":
        return "lengthScaledTPM"
    elif mode == "scaled_tpm":
        return "scaledTPM"
    elif mode == "dtu_scaled_tpm":
        return "dtuScaledTPM"
    return "no"


def _build_r_args_global(cfg: Config, out_dir: Path, mode: str | None = None) -> List[str]:
    """
    Build BASE argument list for a *global* run.
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

    top_n = getattr(cfg, "top_n_vars", 500)
    args.extend(["--top-n", str(top_n)])

    target_files = getattr(cfg, "target_genes_files", [])
    if target_files:
        joined_targets = ",".join(str(tf) for tf in target_files)
        args.extend(["--target-genes", joined_targets])

    deseq2_enabled = bool(getattr(cfg, "deseq2_vst_enabled", True))
    txi_only_mode = bool(getattr(cfg, "tximport_only_mode", False))

    if (not deseq2_enabled) or txi_only_mode:
        args.append("--tximport-only")

    # --- CRITICAL FIX ---
    # Dispersion plot MUST be generated during Phase 1 (Compute) because it needs the DESeq2 object.
    # It cannot be generated in Phase 2 (Plots Only) which only uses the VST matrix.
    if getattr(cfg, "plot_dispersion", False):
        args.append("--dispersion")

    # Always force update in global mode (Phase 1).
    args.append("--force-txi")

    return args


def _build_r_args_for_bp(bp_id: str, cfg: Config, out_dir: Path) -> tuple[List[str], bool]:
    """
    Build argument list for a *per-BioProject* run.
    """
    counts_from_abundance = _map_counts_from_abundance(getattr(cfg, "tximport_mode", None))
    args: List[str] = [
        "Rscript", str(R_SCRIPT_PATH),
        "--search-dir", str(cfg.outdir),
        "--tx2gene", str(cfg.tx2gene),
        "--out-dir", str(out_dir),
        "--prefix", bp_id,
        "--bioproject", bp_id,
        "--counts-from-abundance", counts_from_abundance,
        "--force-txi"
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

    top_n = getattr(cfg, "top_n_vars", 500)
    args.extend(["--top-n", str(top_n)])

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
        return args, True

    # Re-add plot flags for BP
    if deseq2_enabled:
        if getattr(cfg, "plot_pca", False): args.append("--pca")
        if getattr(cfg, "plot_heatmap", False): args.append("--heatmap")
        if getattr(cfg, "plot_var_heatmap", False): args.append("--var-heatmap")
        if getattr(cfg, "plot_sample_cor", False): args.append("--sample-cor")
        if getattr(cfg, "plot_dispersion", False): args.append("--dispersion")

    return args, False


def run_postprocessing_bp(bp: Any, cfg: Config, r_script: Path | None = None) -> None:
    log_path = bp.log_path
    out_dir = bp.path / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        args, is_txi_only = _build_r_args_for_bp(bp.id, cfg, out_dir)
        script_path = r_script if r_script else R_SCRIPT_PATH
        args[1] = str(script_path)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write(f"\n[post-processing] {bp.id} command:\n  {' '.join(args)}\n\n")
            fh.flush()
            subprocess.run(args, stdout=fh, stderr=fh, text=True, check=False)
    except Exception:
        pass

def run_postprocessing(
        dataset: Dataset,
        cfg: Config,
        *,
        r_script: Path | None = None,
        skip_bp: bool | None = None,  # Changed to Optional
        skip_global: bool | None = None  # New Argument
) -> None:
    log_path = cfg.log
    error_warnings: List[str] = cfg.error_warnings

    # Resolve flags (Argument overrides Config)
    do_global = not (skip_global if skip_global is not None else cfg.no_global_postprocessing)
    do_bp = not (skip_bp if skip_bp is not None else cfg.no_bp_postprocessing)

    if cfg.tx2gene is None:
        log_err(error_warnings, log_path, "[post-processing] tx2gene not provided; skipping.")
        return

    script_path = Path(r_script) if r_script is not None else R_SCRIPT_PATH
    if not script_path.exists():
        log_err(error_warnings, log_path, f"[post-processing] R script not found at {script_path}.")
        return

    # ------------------------------------------------------------------
    # 1. Global Analysis
    # ------------------------------------------------------------------
    if do_global:
        out_dir = cfg.shared
        out_dir.mkdir(parents=True, exist_ok=True)

        try:
            r_mode = getattr(dataset, "mode", None)
            base_args = _build_r_args_global(cfg, out_dir, mode=r_mode)
        except Exception as e:
            log_err(error_warnings, log_path, f"[post-processing] Failed to build command: {e}")
            return

        base_args[1] = str(script_path)

        # Determine Active Plots for Phase 2
        active_plots = []
        if getattr(cfg, "plot_pca", False): active_plots.append("--pca")
        if getattr(cfg, "plot_heatmap", False): active_plots.append("--heatmap")
        if getattr(cfg, "plot_var_heatmap", False): active_plots.append("--var-heatmap")
        if getattr(cfg, "plot_sample_cor", False): active_plots.append("--sample-cor")

        # Dispersion is handled in Phase 1 (Compute)

        # --- PHASE 1: COMPUTATION ---
        skip_compute = getattr(cfg, "plots_only_mode", False)

        if not skip_compute:
            log("[post-processing] Phase 1: Calculating Expression Matrices (Global)...", log_path)
            try:
                with open(log_path, "a", encoding="utf-8") as fh:
                    fh.write(f"\n[post-processing] Global Compute Command:\n  {' '.join(base_args)}\n\n")
                    fh.flush()
                    res = subprocess.run(base_args, stdout=fh, stderr=fh, text=True, check=False)
                if res.returncode != 0:
                    log_err(error_warnings, log_path,
                            f"[post-processing] Global Compute failed (code {res.returncode})")
                    # If compute fails, we probably shouldn't try plotting
                    active_plots = []
            except Exception as e:
                log_err(error_warnings, log_path, f"[post-processing] Compute execution failed: {e}")
                active_plots = []
        else:
            log("[post-processing] Phase 1: Skipped (Fast Mode).", log_path)

        # --- PHASE 2: PLOTTING ---
        if active_plots:
            SAFE_PLOT_WORKERS = 2
            log(f"[post-processing] Phase 2: Generating {len(active_plots)} plot types (Concurrency: {SAFE_PLOT_WORKERS})...",
                log_path)

            plot_base_args = [a for a in base_args if a not in ("--force-txi", "--dispersion")]
            plot_base_args.append("--plots-only")

            with ThreadPoolExecutor(max_workers=SAFE_PLOT_WORKERS) as executor:
                future_to_plot = {}
                for plot_flag in active_plots:
                    cmd = plot_base_args + [plot_flag]

                    def run_plot_job(c, name):
                        with open(log_path, "a", encoding="utf-8") as fh:
                            fh.write(f"\n[post-processing] Starting Plot: {name}\n")
                            fh.flush()
                            return subprocess.run(c, stdout=fh, stderr=fh, text=True, check=False)

                    f = executor.submit(run_plot_job, cmd, plot_flag)
                    future_to_plot[f] = plot_flag

                for future in as_completed(future_to_plot):
                    p_flag = future_to_plot[future]
                    try:
                        res = future.result()
                        if res.returncode != 0:
                            log_err(error_warnings, log_path,
                                    f"[post-processing] Plot {p_flag} failed (code {res.returncode})")
                    except Exception as e:
                        log_err(error_warnings, log_path, f"[post-processing] Plot {p_flag} execution error: {e}")

            log("[post-processing] Global plotting finished.", log_path)
    else:
        log("[post-processing] Skipping Global Analysis (--no-global-postprocessing).", log_path)

    # ------------------------------------------------------------------
    # 3. Per-BioProject Run
    # ------------------------------------------------------------------
    if not do_bp:
        log("[post-processing] Skipping per-BioProject analysis (--no-bp-postprocessing).", log_path)
        return

    if not dataset.bioprojects:
        return

    max_workers_bp = max(1, int(cfg.max_threads))
    log(f"[post-processing] Starting analysis for {len(dataset.bioprojects)} BioProjects (Concurrency: {max_workers_bp})...",
        log_path)

    with ThreadPoolExecutor(max_workers=max_workers_bp) as executor:
        # ... (rest of BP logic remains identical) ...
        future_to_bp = {}
        for bp_id in dataset.bioprojects:
            bp_obj = SimpleNamespace(
                id=bp_id,
                path=cfg.outdir / bp_id,
                log_path=cfg.outdir / bp_id / "log.txt"
            )
            bp_obj.path.mkdir(parents=True, exist_ok=True)
            f = executor.submit(run_postprocessing_bp, bp_obj, cfg, r_script=script_path)
            future_to_bp[f] = bp_id

        for future in as_completed(future_to_bp):
            try:
                future.result()
            except Exception:
                pass

    log("[post-processing] All analyses finished.", log_path)