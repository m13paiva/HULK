"""
Seidr Pipeline Integration for HULK.
"""
from __future__ import annotations
import os
import shutil
import sys
import io
import subprocess
import concurrent.futures
from pathlib import Path
from typing import Dict, List, Optional, Any
import pandas as pd

from .entities import Config
from .utils import log

# --- CONSTANTS ---
ENV_OVERRIDES = os.environ.copy()
ENV_OVERRIDES["LC_ALL"] = "C"
ENV_OVERRIDES["LC_NUMERIC"] = "C"
ENV_OVERRIDES["LANG"] = "C"

ALGO_MAP = {
    "PEARSON": ("correlation", "-m", "pearson", "lm"),
    "SPEARMAN": ("correlation", "-m", "spearman", "lm"),
    "PCOR": ("pcor", None, None, "lm"),
    "MI": ("mi", "-m", "RAW", "lm"),
    "CLR": ("mi", "-m", "CLR", "lm"),
    "ARACNE": ("mi", "-m", "ARACNE", "lm"),
    "GENIE3": ("genie3", None, None, "m"),
    "TIGRESS": ("tigress", None, None, "m"),
    "SVM": ("svm-ensemble", "-k", "POLY", "m"),
    "LLR": ("llr-ensemble", None, None, "m"),
    "PLSNET": ("plsnet", None, None, "m"),
    "NARROMI": ("narromi", "-m", "interior-point", "m"),
    "ELNET": ("el-ensemble", None, None, "m"),
}

PRESETS = {
    "FAST": ["PEARSON", "SPEARMAN", "PCOR"],
    "BALANCED": ["PEARSON", "SPEARMAN", "PCOR", "GENIE3", "CLR", "ARACNE"],
    "SLOW": ["PEARSON", "SPEARMAN", "PCOR", "GENIE3", "CLR", "ARACNE",
             "TIGRESS", "SVM", "LLR", "PLSNET", "NARROMI"],
    "ULTRA": ["PEARSON", "SPEARMAN", "PCOR", "GENIE3", "CLR", "ARACNE",
              "TIGRESS", "SVM", "LLR", "PLSNET", "NARROMI", "ELNET"]
}


def _resolve_binaries(needed_algos: List[str], bin_dir: Optional[Path]) -> Dict[str, str]:
    """
    Resolves paths for required Seidr binaries.
    Raises RuntimeError if any are missing.
    """
    if bin_dir:
        # Note: This modifies the process environment.
        # If you run this multiple times with different bin_dirs, it stacks.
        os.environ["PATH"] = str(bin_dir) + os.pathsep + os.environ.get("PATH", "")

    required_bins = {"seidr"}
    for alg in needed_algos:
        if alg in ALGO_MAP:
            required_bins.add(ALGO_MAP[alg][0])

    resolved = {}
    missing = []

    for b in required_bins:
        path = shutil.which(b)
        if path:
            resolved[b] = path
        else:
            missing.append(b)

    if missing:
        # This was previously being caught and swallowed in run_seidr.
        raise RuntimeError(f"Missing Seidr binaries: {', '.join(missing)}")

    return resolved


def _run_direct_visible(cmd: List[str], cwd: Path, log_path: Path, verbose: bool = True):
    cmd_str = " ".join(cmd)
    log(f"[EXEC] {cmd_str}", log_path)

    try:
        res = subprocess.run(
            cmd,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            env=ENV_OVERRIDES,
            check=True
        )
        if not verbose:
            with open(log_path, "a") as f:
                if res.stdout: f.write(res.stdout + "\n")
                if res.stderr: f.write(res.stderr + "\n")
    except subprocess.CalledProcessError as e:
        error_msg = f"\n[Seidr ERROR] Command failed: {cmd_str}\nSTDERR:\n{e.stderr}\nSTDOUT:\n{e.stdout}"
        print(error_msg, file=sys.stderr)
        log(error_msg, log_path)
        raise e


def _import_scores(seidr_bin: str, algo_name: str, outdir: Path, prefix: str,
                   tsv_in: Path, genes: Path, fmt: str, threads: int, log_path: Path, verbose: bool = True) -> Path:
    sf_path = outdir / f"{prefix}{algo_name.lower()}_scores.sf"

    # Check strict existence and size to avoid partial imports
    if sf_path.exists() and sf_path.stat().st_size > 0:
        return sf_path

    cmd = [seidr_bin, "import", "-n", algo_name, "-o", str(sf_path), "-F", fmt, "-i", str(tsv_in), "-g", str(genes)]

    # Ensure options match the algo map behavior
    if algo_name in ["PEARSON", "SPEARMAN", "PCOR"]:
        cmd.extend(["-A", "-r", "-u"])
    elif algo_name == "MI":
        cmd.extend(["-r", "-u", "-O", str(threads)])
    elif algo_name in ["CLR", "ARACNE"]:
        cmd.extend(["-r", "-u", "-z", "-O", str(threads)])
    else:
        cmd.extend(["-r", "-z", "-O", str(threads)])

    try:
        subprocess.run(cmd, cwd=str(outdir), capture_output=True, text=True, env=ENV_OVERRIDES, check=True)
    except subprocess.CalledProcessError as e:
        error_msg = f"\n[Seidr IMPORT ERROR] {algo_name} failed:\n{e.stderr}"
        print(error_msg, file=sys.stderr)
        log(error_msg, log_path)
        raise RuntimeError(f"Seidr import failed for {algo_name}")

    return sf_path


def _safe_float(x: Any) -> float:
    if isinstance(x, (float, int)): return float(x)
    if isinstance(x, str):
        try:
            val = x.split(";")[0]
            return float(val)
        except ValueError:
            return 0.0
    return 0.0


def _export_results(outdir: Path, algorithms: List[str], seidr: str, bb_sf: Path, label: str, no_full: bool,
                    log_path: Path, verbose: bool = True):
    msg = f"[Seidr] Exporting {label} results..."
    log(msg, log_path)
    if verbose: print(msg)

    cmd = [seidr, "view", "--column-headers", str(bb_sf)]
    # Added check=True so we don't silently fail on view errors
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, env=ENV_OVERRIDES, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[Seidr] Export failed: {e.stderr}", file=sys.stderr)
        return

    if not res.stdout.strip():
        print(f"[Seidr] Warning: {label} export resulted in empty output.")
        return

    try:
        df = pd.read_csv(io.StringIO(res.stdout), sep="\t", engine="python")
    except Exception as e:
        print(f"[Seidr] Failed to parse exported table: {e}", file=sys.stderr)
        return

    # Safety check for empty DF
    if df.empty:
        return

    numeric_cols = [c for c in df.columns if c not in ["Source", "Target"] and "interaction" not in c.lower()]
    for c in numeric_cols: df[c] = df[c].apply(_safe_float)

    # Use the last column as the aggregate weight
    agg_col = df.columns[-1]

    simple = df[["Source", "Target", agg_col]].copy().rename(columns={agg_col: "weight"})
    simple["direction"] = df["Interaction"] if "Interaction" in df.columns else "Undirected"
    simple.to_csv(outdir / f"network_{label}_edges.tsv", sep="\t", index=False)

    keepers = ["Source", "Target"]
    if "Interaction" in df.columns: keepers.append("Interaction")

    # Ensure we match algorithm columns case-insensitively
    upper_algos = [a.upper() for a in algorithms]
    for c in df.columns:
        # Check if column starts with any known algo name
        if any(c.upper().startswith(alg) for alg in upper_algos):
            if c not in keepers:
                keepers.append(c)

    df[keepers].to_csv(outdir / f"network_{label}_algs.tsv", sep="\t", index=False)

    if not no_full:
        df.to_csv(outdir / f"network_{label}_edges_full.tsv", sep="\t", index=False)


def _build_network_task(outdir: Path, genes_file: Path, expression_file: Path, threads: int, max_workers: int,
                        backbone: float, aggregate_mode: str, algorithms: List[str], tools: Dict[str, str],
                        label: str, targeted: bool, target_file: Optional[Path], no_full: bool, log_path: Path,
                        force: bool, verbose: bool = True) -> None:
    prefix = f"{label}_"
    seidr = tools["seidr"]
    final_done_marker = outdir / f".{label}.done"
    final_edge_file = outdir / f"network_{label}_edges.tsv"

    # Fixed: Check marker specifically to ensure complete run, not just file existence
    if not force and final_done_marker.exists() and final_edge_file.exists() and final_edge_file.stat().st_size > 0:
        if verbose: print(f"[Seidr] Found existing output for {label}. Skipping.")
        return

    sf_files = []

    def run_algo_task(algo: str, prerequisite_file: Optional[Path] = None) -> Optional[Path]:
        try:
            bin_name, m_flag, m_val, default_fmt = ALGO_MAP[algo]
            current_fmt = "el" if targeted else default_fmt

            # Construct paths
            out_tsv = outdir / f"{prefix}{algo.lower()}_scores.tsv"
            algo_done = outdir / f".{prefix}{algo.lower()}.done"
            out_sf = outdir / f"{prefix}{algo.lower()}_scores.sf"

            if not force and out_tsv.exists() and algo_done.exists() and out_sf.exists():
                return _import_scores(seidr, algo, outdir, prefix, out_tsv, genes_file, current_fmt, threads, log_path,
                                      verbose)

            # Cleanup previous fail
            for p in [out_tsv, algo_done, out_sf]:
                if p.exists(): p.unlink()

            cmd = [tools[bin_name], "-i", str(expression_file), "-g", str(genes_file), "-o", str(out_tsv)]
            if m_flag: cmd.extend([m_flag, m_val])

            # Correlation methods don't support -O threads usually, or handled by BLAS
            if bin_name not in ["correlation", "pcor"]:
                cmd.extend(["-O", str(threads)])
            else:
                cmd.append("--no-scale")

            if targeted:
                cmd.insert(1, "-t")
                cmd.insert(2, str(target_file))

            if algo in ["CLR", "ARACNE"]:
                if not prerequisite_file or not prerequisite_file.exists():
                    raise RuntimeError(f"Prerequisite file for {algo} missing: {prerequisite_file}")
                cmd.extend(["-M", str(prerequisite_file)])

            if verbose: print(f"[Seidr] Running {algo}...")
            _run_direct_visible(cmd, cwd=outdir, log_path=log_path, verbose=verbose)

            algo_done.touch()
            return _import_scores(seidr, algo, outdir, prefix, out_tsv, genes_file, current_fmt, threads, log_path,
                                  verbose)

        except Exception as e:
            # Caught exception inside thread to prevent immediate crash of main thread,
            # but we must print it. _run_direct_visible prints subprocess errors,
            # but we need to print other python errors.
            print(f"[Seidr THREAD ERROR] Algo {algo} failed: {e}", file=sys.stderr)
            return None

    # 1. Run MI first if needed by CLR/ARACNE
    mi_sf_path = None
    run_mi_first = False

    # Check if we need MI for dependencies (CLR/ARACNE) or if MI itself is requested
    if any(x in algorithms for x in ["MI", "CLR", "ARACNE"]):
        # MI is special. It runs if requested OR if needed as dependency.
        # If MI is not in 'algorithms' but CLR is, we run MI (but don't necessarily add to sf_files if user didn't ask for it?)
        # Actually standard practice: if CLR is run, MI matrix is intermediate.
        # If MI is in algorithms, we add it to sf_files.
        run_mi_first = True

    if run_mi_first:
        # Note: run_algo_task("MI") returns the SF file.
        # For CLR/ARACNE we need the TSV file usually, or the Seidr import?
        # Seidr tools (mi -m CLR) take the TSV matrix via -M.
        # run_algo_task creates TSV at: outdir / f"{prefix}mi_scores.tsv"
        mi_sf = run_algo_task("MI")
        if "MI" in algorithms and mi_sf:
            sf_files.append(mi_sf)

    # 2. Run everything else in parallel
    parallel_algos = [a for a in algorithms if a != "MI"]

    if parallel_algos:
        # Cast max_workers to int to avoid TypeErrors if config passed a string
        safe_workers = int(max_workers) if max_workers else 1

        with concurrent.futures.ThreadPoolExecutor(max_workers=safe_workers) as executor:
            future_map = {}
            for algo in parallel_algos:
                prereq = None
                if algo in ["CLR", "ARACNE"]:
                    # These need the MI TSV file, not the SF file
                    prereq = outdir / f"{prefix}mi_scores.tsv"

                future_map[executor.submit(run_algo_task, algo, prereq)] = algo

            for future in concurrent.futures.as_completed(future_map):
                res = future.result()
                if res: sf_files.append(res)

    if not sf_files:
        print("[Seidr] No algorithms finished successfully. Aborting aggregation.", file=sys.stderr)
        return

    net_sf = outdir / f"network_{label}.sf"
    if verbose: print(f"[Seidr] Aggregating {len(sf_files)} networks...")

    try:
        subprocess.run([seidr, "aggregate", "-o", str(net_sf), "-m", aggregate_mode] + [str(s) for s in sf_files],
                       cwd=str(outdir), capture_output=True, env=ENV_OVERRIDES, check=True)
    except subprocess.CalledProcessError as e:
        print(f"\n[Seidr AGGREGATE ERROR]\n{e.stderr}", file=sys.stderr)
        return

    bb_sf = outdir / f"network_{label}.bb.sf"
    try:
        # Backbone also needs explicit error handling
        _run_direct_visible([seidr, "backbone", "-F", str(backbone), str(net_sf)], cwd=outdir, log_path=log_path,
                            verbose=verbose)
    except Exception:
        # Error already logged by _run_direct_visible
        return

    if bb_sf.exists():
        _export_results(outdir, algorithms, seidr, bb_sf, label, no_full, log_path, verbose)

        # Cleanup
        for p in outdir.glob(f"{prefix}*_scores.*"): p.unlink()
        if net_sf.exists(): net_sf.unlink()
        if bb_sf.exists(): bb_sf.unlink()

        final_done_marker.touch()


def run_seidr_batch(cfg: Config, genes_file: Path, expression_file: Path, outdir: Path, preset: str = "FAST",
                    threads: Optional[int] = None) -> None:
    log_path = outdir / "seidr_batch.log"

    # Ensure algos are uppercase for ALGO_MAP lookup
    algos = [a.upper() for a in PRESETS.get(preset.upper(), PRESETS["FAST"])]

    # REMOVED: try-except block that swallowed errors. Now it will raise RuntimeError.
    tools = _resolve_binaries(algos, None)

    # max_workers=1 for batch usually, to avoid overloading since batch might be parallelized externally?
    # Or passed as arg. Here hardcoded to 1.
    _build_network_task(outdir, genes_file, expression_file, (threads or cfg.max_threads), 1, 1.28, "irp", algos, tools,
                        "saturation", False, None, True, log_path, True, False)

    for d in [p for p in outdir.glob("*-*-*-*") if p.is_dir()]:
        shutil.rmtree(d, ignore_errors=True)


def run_seidr(cfg: Config, force: bool = False) -> None:
    opts = cfg.get_tool_opts("seidr")

    # Check enabled status (handle string 'true' just in case yaml parser didn't)
    enabled = opts.get("enabled", False)
    if isinstance(enabled, str):
        enabled = enabled.lower() == "true"
    if not enabled:
        return

    # Determine algorithms
    preset_name = opts.get("preset", "BALANCED")
    raw_algos = opts.get("algorithms", PRESETS.get(preset_name, PRESETS["BALANCED"]))

    # Sanitize inputs: Uppercase algorithms
    algos = [a.upper() for a in raw_algos]

    outdir = Path(opts["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)

    # REMOVED: Silent try-except. If tools are missing, CRASH and tell the user.
    tools = _resolve_binaries(algos, None)

    # Ensure workers is int
    workers = opts.get("workers", 4)
    try:
        workers = int(workers)
    except ValueError:
        workers = 4

    task_args = {
        "outdir": outdir,
        "genes_file": Path(opts["genes_file"]),
        "expression_file": Path(opts["expression_file"]),
        "threads": cfg.max_threads,
        "max_workers": workers,
        "backbone": float(opts.get("backbone", 1.28)),
        "aggregate_mode": opts.get("aggregate", "irp"),
        "algorithms": algos,
        "tools": tools,
        "no_full": opts.get("no_full", False),
        "log_path": Path(getattr(cfg, "log", "seidr.log")),
        "force": (force or opts.get("force", False)),
        "verbose": True
    }

    _build_network_task(label="main", targeted=False, target_file=None, **task_args)