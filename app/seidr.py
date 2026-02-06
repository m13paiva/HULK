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
from tqdm import tqdm

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


# --- INTERNAL HELPERS ---

def _resolve_binaries(needed_algos: List[str], bin_dir: Optional[Path]) -> Dict[str, str]:
    if bin_dir:
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
        raise RuntimeError(f"Missing Seidr binaries: {', '.join(missing)}")

    return resolved


def _run_direct_visible(cmd: List[str], cwd: Path, log_path: Path, verbose: bool = True):
    """
    Executes a command.
    If verbose=False, redirects ALL output (stdout+stderr) to the log file to silence console.
    """
    cmd_str = " ".join(cmd)
    log(f"[EXEC] {cmd_str}", log_path)

    file_handle = None

    if verbose:
        # Standard behavior: console sees stderr (e.g., Seidr progress)
        stdout_target = subprocess.DEVNULL
        stderr_target = None
    else:
        # Silent behavior: Append stderr directly to log file
        try:
            file_handle = open(log_path, "a")
            stdout_target = file_handle
            stderr_target = file_handle
        except:
            stdout_target = subprocess.DEVNULL
            stderr_target = subprocess.DEVNULL

    try:
        subprocess.run(
            cmd,
            cwd=str(cwd),
            stdout=stdout_target,
            stderr=stderr_target,
            env=ENV_OVERRIDES,
            check=True
        )
    except subprocess.CalledProcessError as e:
        log(f"[ERROR] Command failed with exit code {e.returncode}", log_path)
        raise e
    finally:
        if file_handle:
            file_handle.close()


def _import_scores(seidr_bin: str, algo_name: str, outdir: Path, prefix: str,
                   tsv_in: Path, genes: Path, fmt: str, threads: int, log_path: Path, verbose: bool = True) -> Path:
    sf_path = outdir / f"{prefix}{algo_name.lower()}_scores.sf"

    if sf_path.exists() and sf_path.stat().st_size > 0:
        return sf_path

    cmd = [seidr_bin, "import", "-n", algo_name, "-o", str(sf_path),
           "-F", fmt, "-i", str(tsv_in), "-g", str(genes)]

    if algo_name in ["PEARSON", "SPEARMAN", "PCOR"]:
        cmd.extend(["-A", "-r", "-u"])
    elif algo_name == "MI":
        cmd.extend(["-r", "-u", "-O", str(threads)])
    elif algo_name in ["CLR", "ARACNE"]:
        cmd.extend(["-r", "-u", "-z", "-O", str(threads)])
    else:
        cmd.extend(["-r", "-z", "-O", str(threads)])

    cmd_str = " ".join(cmd)
    log(f"[EXEC] {cmd_str}", log_path)

    try:
        # Always capture import output to keep it clean (it's very verbose)
        subprocess.run(
            cmd,
            cwd=str(outdir),
            capture_output=True,
            text=True,
            env=ENV_OVERRIDES,
            check=True
        )
    except subprocess.CalledProcessError as e:
        log(f"[ERROR] Import failed: {e.stderr}", log_path)
        raise RuntimeError(f"Seidr import failed for {algo_name}: {e.stderr}")

    if not sf_path.exists():
        raise RuntimeError(f"Import failed for {algo_name}. Output {sf_path} missing.")

    return sf_path


def _safe_float(x: Any) -> float:
    if isinstance(x, (float, int)):
        return float(x)
    if isinstance(x, str):
        try:
            val = x.split(";")[0]
            return float(val)
        except (ValueError, TypeError):
            return 0.0
    return 0.0


def _export_results(outdir: Path, algorithms: List[str], seidr: str, bb_sf: Path, label: str, no_full: bool,
                    log_path: Path, verbose: bool = True):
    if verbose:
        msg = f"[Seidr] Exporting {label} results..."
        print(msg)
    else:
        msg = f"[Seidr] Exporting {label} results..."

    log(msg, log_path)

    cmd = [seidr, "view", "--column-headers", str(bb_sf)]
    log(f">> {' '.join(cmd)} (capturing output)", log_path)

    res = subprocess.run(cmd, capture_output=True, text=True, env=ENV_OVERRIDES)

    if not res.stdout.strip():
        msg = f"[Warn] No edges found in {bb_sf}. Backbone threshold too strict?"
        if verbose: print(msg)
        log(msg, log_path)
        return

    try:
        df = pd.read_csv(io.StringIO(res.stdout), sep="\t", engine="python")
    except Exception as e:
        msg = f"[Error] Failed to parse Seidr output: {e}"
        if verbose: print(msg)
        log(msg, log_path)
        return

    numeric_cols = [c for c in df.columns if c not in ["Source", "Target"] and "interaction" not in c.lower()]
    for c in numeric_cols:
        df[c] = df[c].apply(_safe_float)

    agg_col = df.columns[-1]
    simple = df[["Source", "Target", agg_col]].copy().rename(columns={agg_col: "weight"})
    simple["direction"] = df["Interaction"] if "Interaction" in df.columns else "Undirected"
    simple.to_csv(outdir / f"network_{label}_edges.tsv", sep="\t", index=False)

    keepers = ["Source", "Target"]
    if "Interaction" in df.columns: keepers.append("Interaction")
    for alg in algorithms:
        matches = [c for c in df.columns if c.upper().startswith(alg)]
        keepers.extend(matches)

    df[keepers].to_csv(outdir / f"network_{label}_algs.tsv", sep="\t", index=False)

    if not no_full:
        df.to_csv(outdir / f"network_{label}_edges_full.tsv", sep="\t", index=False)


def _build_network_task(
        outdir: Path,
        genes_file: Path,
        expression_file: Path,
        threads: int,
        max_workers: int,
        backbone: float,
        aggregate_mode: str,
        algorithms: List[str],
        tools: Dict[str, str],
        label: str,
        targeted: bool,
        target_file: Optional[Path],
        no_full: bool,
        log_path: Path,
        force: bool,
        verbose: bool = True
) -> None:
    prefix = f"{label}_"
    seidr = tools["seidr"]

    final_done_marker = outdir / f".{label}.done"
    final_edge_file = outdir / f"network_{label}_edges.tsv"

    # --- 0. RETROACTIVE COMPLETION CHECK ---
    # If edge file exists, we consider it done. Create marker if missing.
    if not force and final_edge_file.exists() and final_edge_file.stat().st_size > 0:
        if not final_done_marker.exists():
            final_done_marker.touch()
            # Also auto-clean just in case
            for p in outdir.glob(f"{prefix}*_scores.sf"): p.unlink()
            for p in outdir.glob(f"{prefix}*_scores.tsv"): p.unlink()
            bb = outdir / f"network_{label}.bb.sf"
            if bb.exists(): bb.unlink()

        if verbose: print(f"[Seidr] Found existing output for {label}. Skipping.")
        log(f"[Seidr] Found existing output for {label}. Marked done.", log_path)
        return

    sf_files = []

    def run_algo_task(algo: str, prerequisite_file: Optional[Path] = None) -> Optional[Path]:
        try:
            bin_name, m_flag, m_val, default_fmt = ALGO_MAP[algo]
            current_fmt = "el" if targeted else default_fmt

            out_tsv = outdir / f"{prefix}{algo.lower()}_scores.tsv"
            algo_done = outdir / f".{prefix}{algo.lower()}.done"
            out_sf = outdir / f"{prefix}{algo.lower()}_scores.sf"

            if not force and out_tsv.exists() and algo_done.exists() and out_sf.exists():
                if verbose: print(f"[Seidr] Found verified cache for {algo}. Skipping.")
                log(f"[Seidr] Found verified cache for {algo}", log_path)
                return _import_scores(seidr, algo, outdir, prefix, out_tsv, genes_file, current_fmt, threads, log_path,
                                      verbose)

            if out_tsv.exists(): out_tsv.unlink()
            if algo_done.exists(): algo_done.unlink()
            if out_sf.exists(): out_sf.unlink()

            cmd = [tools[bin_name]]
            cmd.extend(["-i", str(expression_file), "-g", str(genes_file), "-o", str(out_tsv)])

            if m_flag: cmd.extend([m_flag, m_val])
            if bin_name not in ["correlation", "pcor"]:
                cmd.extend(["-O", str(threads)])
            else:
                cmd.append("--no-scale")
            if targeted:
                if not target_file: raise ValueError("Targeted mode without file")
                cmd.insert(1, "-t")
                cmd.insert(2, str(target_file))
            if algo in ["CLR", "ARACNE"]:
                cmd.extend(["-M", str(prerequisite_file)])

            if verbose: print(f"[Seidr] Running {algo}...")

            # Use the silence-aware runner
            _run_direct_visible(cmd, cwd=outdir, log_path=log_path, verbose=verbose)

            if not out_tsv.exists():
                raise FileNotFoundError(f"{algo} reported success but {out_tsv} is missing.")

            algo_done.touch()
            return _import_scores(seidr, algo, outdir, prefix, out_tsv, genes_file, current_fmt, threads, log_path,
                                  verbose)

        except Exception as e:
            msg = f"[Error] {algo} failed: {e}"
            if verbose: print(msg)
            log(msg, log_path)
            d = outdir / f".{prefix}{algo.lower()}.done"
            if d.exists(): d.unlink()
            return None

    # 1. Dependency: MI
    mi_sf_path = None
    mi_tsv_path = outdir / f"{prefix}mi_scores.tsv"

    if any(x in algorithms for x in ["MI", "CLR", "ARACNE"]):
        mi_sf_path = run_algo_task("MI")
        if "MI" in algorithms and mi_sf_path:
            sf_files.append(mi_sf_path)

    # 2. Parallel Execution
    parallel_algos = [a for a in algorithms if a != "MI"]

    if parallel_algos:
        if verbose and max_workers > 1:
            print("\nWARNING: Workers > 1. Multiple progress bars will garble your console output.")

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {}
            for algo in parallel_algos:
                prereq = mi_tsv_path if algo in ["CLR", "ARACNE"] else None
                if prereq and not prereq.exists():
                    log(f"[Error] Cannot run {algo}: MI output missing.", log_path)
                    continue
                future = executor.submit(run_algo_task, algo, prereq)
                future_map[future] = algo

            if verbose and len(parallel_algos) > 1:
                iterator = tqdm(concurrent.futures.as_completed(future_map), total=len(parallel_algos), desc="Algos",
                                unit="algo", leave=False)
            else:
                iterator = concurrent.futures.as_completed(future_map)

            for future in iterator:
                res = future.result()
                if res: sf_files.append(res)

    # 3. Aggregate
    if not sf_files:
        return

    net_sf = outdir / f"network_{label}.sf"
    if force and net_sf.exists(): net_sf.unlink()

    if verbose: print(f"[Seidr] Aggregating {len(sf_files)} networks...")

    agg_cmd = [seidr, "aggregate", "-o", str(net_sf), "-m", aggregate_mode] + [str(s) for s in sf_files]

    try:
        cmd_str = " ".join(agg_cmd)
        log(f"[EXEC] {cmd_str}", log_path)

        # Capture stderr manually if verbose=False to ensure silence
        res = subprocess.run(agg_cmd, cwd=str(outdir), capture_output=True, text=True, check=True, env=ENV_OVERRIDES)

        # If verbose=True, we didn't print it above anyway (agg is fast), but let's log it
        with open(log_path, "a") as f:
            f.write(res.stderr + "\n")

    except Exception as e:
        log(f"[Error] Aggregation failed: {e}", log_path)
        return

    # 4. Backbone
    if verbose: print(f"[Seidr] Pruning (Backbone {backbone})...")

    bb_sf = outdir / f"network_{label}.bb.sf"
    if force and bb_sf.exists(): bb_sf.unlink()

    try:
        cmd = [seidr, "backbone", "-F", str(backbone), str(net_sf)]
        _run_direct_visible(cmd, cwd=outdir, log_path=log_path, verbose=verbose)
    except Exception as e:
        log(f"[Error] Backbone failed: {e}", log_path)
        return

    if bb_sf.exists():
        _export_results(outdir, algorithms, seidr, bb_sf, label, no_full, log_path, verbose)

        # --- CLEANUP & MARK DONE ---
        if verbose: print(f"[Seidr] Cleaning intermediate files...")

        for p in outdir.glob(f"{prefix}*_scores.sf"): p.unlink()
        for p in outdir.glob(f"{prefix}*_scores.tsv"): p.unlink()

        if bb_sf.exists(): bb_sf.unlink()

        # Mark as done
        final_done_marker.touch()


def run_seidr_batch(cfg: Config, genes_file: Path, expression_file: Path, outdir: Path, preset: str = "FAST",
                    threads: Optional[int] = None) -> None:
    log_path = outdir / "seidr_batch.log"
    algos = PRESETS.get(preset.upper(), PRESETS["FAST"])
    try:
        tools = _resolve_binaries(algos, None)
    except RuntimeError as e:
        with open(log_path, "a") as f:
            f.write(f"[FATAL] {e}\n")
        return
    effective_threads = threads if threads is not None else cfg.max_threads

    # CALLER WITH SILENCE
    _build_network_task(
        outdir=outdir,
        genes_file=genes_file,
        expression_file=expression_file,
        threads=effective_threads,
        max_workers=1,
        backbone=1.28,
        aggregate_mode="irp",
        algorithms=algos,
        tools=tools,
        label="saturation",
        targeted=False,
        target_file=None,
        no_full=True,
        log_path=log_path,
        force=True,
        verbose=False  # <--- Ensures Silence
    )

    junk_dirs = [p for p in outdir.glob("*-*-*-*") if p.is_dir()]
    for d in junk_dirs:
        try:
            shutil.rmtree(d)
        except:
            pass


def run_seidr(cfg: Config, force: bool = False) -> None:
    log_path = Path(getattr(cfg, "log", "seidr.log"))
    opts = cfg.get_tool_opts("seidr")
    if not opts.get("enabled", False): return
    active_force = force or opts.get("force", False)
    genes_file = Path(opts["genes_file"])
    expression_file = Path(opts["expression_file"])
    outdir = Path(opts["outdir"])
    targets = [Path(t) for t in opts.get("targets", [])]
    target_mode = opts.get("target_mode", "both")
    algos = opts.get("algorithms", [])
    if not algos: algos = PRESETS.get(opts.get("preset", "BALANCED"), PRESETS["BALANCED"])
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        tools = _resolve_binaries(algos, None)
    except RuntimeError as e:
        print(f"[Seidr] Error: {e}"); return

    print("=" * 60)
    print("STARTING SEIDR NETWORK INFERENCE" + (" (FORCE MODE)" if active_force else ""))

    task_args = {
        "outdir": outdir, "genes_file": genes_file, "expression_file": expression_file,
        "threads": cfg.max_threads, "max_workers": opts["workers"], "backbone": opts["backbone"],
        "aggregate_mode": opts.get("aggregate", "irp"), "algorithms": algos, "tools": tools,
        "no_full": opts.get("no_full", False), "log_path": log_path, "force": active_force, "verbose": True
    }
    if not targets or target_mode in ["both", "main_only"]:
        _build_network_task(label="main", targeted=False, target_file=None, **task_args)
    if targets and target_mode in ["both", "targeted_only"]:
        for t_file in targets:
            if t_file.exists(): _build_network_task(label=t_file.stem, targeted=True, target_file=t_file, **task_args)

    junk_dirs = [p for p in outdir.glob("*-*-*-*") if p.is_dir()]
    for d in junk_dirs:
        try:
            shutil.rmtree(d)
        except Exception:
            pass
    print("[Seidr] Analysis Finished.")