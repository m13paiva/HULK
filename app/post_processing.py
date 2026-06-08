from __future__ import annotations
import subprocess
import shutil
from pathlib import Path
from types import SimpleNamespace
from typing import List, Any, Optional, Set
from concurrent.futures import ThreadPoolExecutor, as_completed

from .entities import Config, Dataset
from .utils import log

R_SCRIPT_PATH = Path(__file__).with_name("post_processing.R")


def _map_counts_from_abundance(mode: str | None) -> str:
    """
    Maps the specified quantification mode to its corresponding tximport parameter string.

    Args:
        mode (str | None): The quantification mode (e.g., 'length_scaled_tpm').

    Returns:
        str: The mapped string recognized by the tximport R package. Defaults to "no".
    """
    if mode == "length_scaled_tpm":
        return "lengthScaledTPM"
    elif mode == "scaled_tpm":
        return "scaledTPM"
    elif mode == "dtu_scaled_tpm":
        return "dtuScaledTPM"
    return "no"


def _build_r_args_global(cfg: Config, out_dir: Path, mode: str | None = None) -> List[str]:
    """
    Constructs the base argument list for a global execution of the post-processing R script.

    Args:
        cfg (Config): The configuration object containing runtime parameters.
        out_dir (Path): The designated output directory for the global run.
        mode (str | None): The operational mode (e.g., "SRR" or "FASTQ").

    Raises:
        ValueError: If the `tx2gene` mapping file is missing from the configuration.

    Returns:
        List[str]: A list of command-line arguments formatted for subprocess execution.
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

    if getattr(cfg, "plot_dispersion", False):
        args.append("--dispersion")

    # Enforce forced re-computation of tximport structures in global aggregation context.
    args.append("--force-txi")

    return args


def _build_r_args_for_bp(bp_id: str, cfg: Config, out_dir: Path) -> tuple[List[str], bool]:
    """
    Constructs the argument list for a per-BioProject execution of the post-processing R script.

    Args:
        bp_id (str): The specific BioProject identifier.
        cfg (Config): The configuration object containing runtime parameters.
        out_dir (Path): The designated output directory for the isolated run.

    Returns:
        tuple[List[str], bool]: A tuple containing the argument list and a boolean flag
        indicating whether the execution is operating in tximport-only mode.
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

    if deseq2_enabled:
        if getattr(cfg, "plot_pca", False): args.append("--pca")
        if getattr(cfg, "plot_heatmap", False): args.append("--heatmap")
        if getattr(cfg, "plot_var_heatmap", False): args.append("--var-heatmap")
        if getattr(cfg, "plot_sample_cor", False): args.append("--sample-cor")
        if getattr(cfg, "plot_dispersion", False): args.append("--dispersion")

    return args, False


def run_postprocessing_bp(bp: Any, cfg: Config, r_script: Path | None = None) -> Path | None:
    """
    Executes post-processing logic for a single, isolated BioProject block.

    Args:
        bp (Any): A BioProject instance or generic namespace mapping.
        cfg (Config): Runtime configuration.
        r_script (Path | None): Explicit override path to the R script.

    Returns:
        Path | None: The path to the generated plots directory if successful, otherwise None.
    """
    log_path = cfg.log
    if hasattr(bp, 'log_path'):
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
        return out_dir
    except Exception:
        return None


def run_postprocessing(
        dataset: Dataset,
        cfg: Config,
        *,
        r_script: Path | None = None,
        skip_bp: bool | None = None,
        skip_global: bool | None = None
) -> None:
    """
    Manages the overarching execution schedule for all post-processing analytics,
    distinguishing between global summaries and distinct BioProject reports.

    Args:
        dataset (Dataset): Loaded dataset tracking the execution context.
        cfg (Config): Configuration profile tracking active analytics flags.
        r_script (Path | None): Override script location.
        skip_bp (bool | None): Override boolean blocking per-project computations.
        skip_global (bool | None): Override boolean blocking global computations.
    """
    log_path = cfg.log

    do_global = not (skip_global if skip_global is not None else cfg.no_global_postprocessing)
    do_bp = not (skip_bp if skip_bp is not None else cfg.no_bp_postprocessing)

    if cfg.tx2gene is None:
        print("[post-processing] tx2gene not provided; skipping.")
        return

    script_path = Path(r_script) if r_script is not None else R_SCRIPT_PATH
    if not script_path.exists():
        print(f"[post-processing] R script not found at {script_path}.")
        return

    # ─────────────────────────────────────────────────────────────────────────────
    # 1. Global Analysis Phase
    # ─────────────────────────────────────────────────────────────────────────────
    if do_global:
        out_dir = cfg.shared
        out_dir.mkdir(parents=True, exist_ok=True)

        try:
            r_mode = getattr(dataset, "mode", None)
            base_args = _build_r_args_global(cfg, out_dir, mode=r_mode)
        except Exception as e:
            print(f"[post-processing] Failed to build command: {e}")
            return

        base_args[1] = str(script_path)

        active_plots = []
        if getattr(cfg, "plot_pca", False): active_plots.append("--pca")
        if getattr(cfg, "plot_heatmap", False): active_plots.append("--heatmap")
        if getattr(cfg, "plot_var_heatmap", False): active_plots.append("--var-heatmap")
        if getattr(cfg, "plot_sample_cor", False): active_plots.append("--sample-cor")

        skip_compute = getattr(cfg, "plots_only_mode", False)

        if not skip_compute:
            log("[post-processing] Phase 1: Calculating Expression Matrices (Global)...", log_path)
            try:
                with open(log_path, "a", encoding="utf-8") as fh:
                    fh.write(f"\n[post-processing] Global Compute Command:\n  {' '.join(base_args)}\n\n")
                    fh.flush()
                    res = subprocess.run(base_args, stdout=fh, stderr=fh, text=True, check=False)
                if res.returncode != 0:
                    print(f"[post-processing] Global Compute failed (code {res.returncode})")
                    active_plots = []
            except Exception as e:
                print(f"[post-processing] Compute execution failed: {e}")
                active_plots = []
        else:
            log("[post-processing] Phase 1: Skipped (Fast Mode).", log_path)

        # --- Plot Generation ---
        if active_plots:
            SAFE_PLOT_WORKERS = 2
            log(f"[post-processing] Phase 2: Generating {len(active_plots)} plot types (Concurrency: {SAFE_PLOT_WORKERS})...", log_path)

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
                            print(f"[post-processing] Plot {p_flag} failed (code {res.returncode})")
                    except Exception as e:
                        print(f"[post-processing] Plot {p_flag} execution error: {e}")

            log("[post-processing] Global plotting finished.", log_path)
    else:
        log("[post-processing] Skipping Global Analysis (--no-global-postprocessing).", log_path)

    # ─────────────────────────────────────────────────────────────────────────────
    # 2. Per-BioProject Phase
    # ─────────────────────────────────────────────────────────────────────────────
    if not do_bp:
        log("[post-processing] Skipping per-BioProject analysis (--no-bp-postprocessing).", log_path)
        return

    if not dataset.bioprojects:
        return

    max_workers_bp = max(1, int(getattr(cfg, "max_threads", 4)))
    log(f"[post-processing] Starting analysis for {len(dataset.bioprojects)} BioProjects (Concurrency: {max_workers_bp})...", log_path)

    with ThreadPoolExecutor(max_workers=max_workers_bp) as executor:
        future_to_bp = {}
        for bp_item in dataset.bioprojects:
            bp_id_str = str(bp_item)
            if hasattr(bp_item, "id"):
                bp_id_str = bp_item.id

            bp_obj = SimpleNamespace(
                id=bp_id_str,
                path=cfg.outdir / bp_id_str,
                log_path=cfg.outdir / bp_id_str / "log.txt"
            )
            bp_obj.path.mkdir(parents=True, exist_ok=True)
            f = executor.submit(run_postprocessing_bp, bp_obj, cfg, r_script=script_path)
            future_to_bp[f] = bp_id_str

        for future in as_completed(future_to_bp):
            try:
                future.result()
            except Exception:
                pass

    log("[post-processing] All analyses finished.", log_path)