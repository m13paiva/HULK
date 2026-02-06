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
    cmd_str = " ".join(cmd)
    log(f"[EXEC] {cmd_str}", log_path)

    # Run and capture output to handle errors properly
    try:
        res = subprocess.run(
            cmd,
            cwd=str(cwd),
            capture_output=True,
            text=True,
            env=ENV_OVERRIDES,
            check=True
        )
        # If silent but successful, we just log the output
        if not verbose:
            with open(log_path, "a") as f:
                if res.stdout: f.write(res.stdout + "\n")
                if res.stderr: f.write(res.stderr + "\n")
    except subprocess.CalledProcessError as e:
        # ALWAYS print errors to console even if verbose=False
        error_msg = f"\n[Seidr ERROR] Command failed: {cmd_str}\n{e.stderr}"
        print(error_msg, file=sys.stderr)
        log(error_msg, log_path)
        raise e


def _import_scores(seidr_bin: str, algo_name: str, outdir: Path, prefix: str,
                   tsv_in: Path, genes: Path, fmt: str, threads: int, log_path: Path, verbose: bool = True) -> Path:
    sf_path = outdir / f"{prefix}{algo_name.lower()}_scores.sf"
    if sf_path.exists() and sf_path.stat().st_size > 0: return sf_path
    cmd = [seidr_bin, "import", "-n", algo_name, "-o", str(sf_path), "-F", fmt, "-i", str(tsv_in), "-g", str(genes)]
    if algo_name in ["PEARSON", "SPEARMAN", "PCOR"]:
        cmd.extend(["-A", "-r", "-u"])
    elif algo_name == "MI":
        cmd.extend(["-r", "-u", "-O", str(threads)])
    elif algo_name in ["CLR", "ARACNE"]:
        cmd.extend(["-r", "-u", "-z", "-O", str(threads)])
    else:
        cmd.extend(["-r", "-z", "-O", str(threads)])

    try:
        # Import is always silent unless it fails
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
        except:
            return 0.0
    return 0.0


def _export_results(outdir: Path, algorithms: List[str], seidr: str, bb_sf: Path, label: str, no_full: bool,
                    log_path: Path, verbose: bool = True):
    msg = f"[Seidr] Exporting {label} results..."
    log(msg, log_path)
    if verbose: print(msg)
    cmd = [seidr, "view", "--column-headers", str(bb_sf)]
    res = subprocess.run(cmd, capture_output=True, text=True, env=ENV_OVERRIDES)
    if not res.stdout.strip(): return
    try:
        df = pd.read_csv(io.StringIO(res.stdout), sep="\t", engine="python")
    except:
        return
    numeric_cols = [c for c in df.columns if c not in ["Source", "Target"] and "interaction" not in c.lower()]
    for c in numeric_cols: df[c] = df[c].apply(_safe_float)
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
    if not no_full: df.to_csv(outdir / f"network_{label}_edges_full.tsv", sep="\t", index=False)


def _build_network_task(outdir: Path, genes_file: Path, expression_file: Path, threads: int, max_workers: int,
                        backbone: float, aggregate_mode: str, algorithms: List[str], tools: Dict[str, str],
                        label: str, targeted: bool, target_file: Optional[Path], no_full: bool, log_path: Path,
                        force: bool, verbose: bool = True) -> None:
    prefix = f"{label}_"
    seidr = tools["seidr"]
    final_done_marker = outdir / f".{label}.done"
    final_edge_file = outdir / f"network_{label}_edges.tsv"

    if not force and final_edge_file.exists() and final_edge_file.stat().st_size > 0:
        if not final_done_marker.exists(): final_done_marker.touch()
        if verbose: print(f"[Seidr] Found existing output for {label}. Skipping.")
        return

    sf_files = []

    def run_algo_task(algo: str, prerequisite_file: Optional[Path] = None) -> Optional[Path]:
        try:
            bin_name, m_flag, m_val, default_fmt = ALGO_MAP[algo]
            current_fmt = "el" if targeted else default_fmt
            out_tsv, algo_done, out_sf = outdir / f"{prefix}{algo.lower()}_scores.tsv", outdir / f".{prefix}{algo.lower()}.done", outdir / f"{prefix}{algo.lower()}_scores.sf"

            if not force and out_tsv.exists() and algo_done.exists() and out_sf.exists():
                return _import_scores(seidr, algo, outdir, prefix, out_tsv, genes_file, current_fmt, threads, log_path,
                                      verbose)

            for p in [out_tsv, algo_done, out_sf]:
                if p.exists(): p.unlink()

            cmd = [tools[bin_name], "-i", str(expression_file), "-g", str(genes_file), "-o", str(out_tsv)]
            if m_flag: cmd.extend([m_flag, m_val])
            if bin_name not in ["correlation", "pcor"]:
                cmd.extend(["-O", str(threads)])
            else:
                cmd.append("--no-scale")
            if targeted: cmd.insert(1, "-t"); cmd.insert(2, str(target_file))
            if algo in ["CLR", "ARACNE"]: cmd.extend(["-M", str(prerequisite_file)])

            if verbose: print(f"[Seidr] Running {algo}...")
            _run_direct_visible(cmd, cwd=outdir, log_path=log_path, verbose=verbose)
            algo_done.touch()
            return _import_scores(seidr, algo, outdir, prefix, out_tsv, genes_file, current_fmt, threads, log_path,
                                  verbose)
        except Exception as e:
            # We don't catch here so the exception propagates and prints in _run_direct_visible
            return None

    mi_sf_path = None
    if any(x in algorithms for x in ["MI", "CLR", "ARACNE"]):
        mi_sf_path = run_algo_task("MI")
        if "MI" in algorithms and mi_sf_path: sf_files.append(mi_sf_path)

    parallel_algos = [a for a in algorithms if a != "MI"]
    if parallel_algos:
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {executor.submit(run_algo_task, algo, (
                outdir / f"{prefix}mi_scores.tsv" if algo in ["CLR", "ARACNE"] else None)): algo for algo in
                          parallel_algos}
            for future in concurrent.futures.as_completed(future_map):
                res = future.result()
                if res: sf_files.append(res)

    if not sf_files: return

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
        _run_direct_visible([seidr, "backbone", "-F", str(backbone), str(net_sf)], cwd=outdir, log_path=log_path,
                            verbose=verbose)
    except:
        return

    if bb_sf.exists():
        _export_results(outdir, algorithms, seidr, bb_sf, label, no_full, log_path, verbose)
        for p in outdir.glob(f"{prefix}*_scores.*"): p.unlink()
        if net_sf.exists(): net_sf.unlink()
        if bb_sf.exists(): bb_sf.unlink()
        final_done_marker.touch()


def run_seidr_batch(cfg: Config, genes_file: Path, expression_file: Path, outdir: Path, preset: str = "FAST",
                    threads: Optional[int] = None) -> None:
    log_path = outdir / "seidr_batch.log"
    algos = PRESETS.get(preset.upper(), PRESETS["FAST"])
    try:
        tools = _resolve_binaries(algos, None)
    except RuntimeError:
        return
    _build_network_task(outdir, genes_file, expression_file, (threads or cfg.max_threads), 1, 1.28, "irp", algos, tools,
                        "saturation", False, None, True, log_path, True, False)
    for d in [p for p in outdir.glob("*-*-*-*") if p.is_dir()]: shutil.rmtree(d, ignore_errors=True)


def run_seidr(cfg: Config, force: bool = False) -> None:
    opts = cfg.get_tool_opts("seidr")
    if not opts.get("enabled", False): return
    algos = opts.get("algorithms", PRESETS.get(opts.get("preset", "BALANCED"), PRESETS["BALANCED"]))
    outdir = Path(opts["outdir"]);
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        tools = _resolve_binaries(algos, None)
    except RuntimeError:
        return
    task_args = {"outdir": outdir, "genes_file": Path(opts["genes_file"]),
                 "expression_file": Path(opts["expression_file"]), "threads": cfg.max_threads,
                 "max_workers": opts["workers"], "backbone": opts["backbone"],
                 "aggregate_mode": opts.get("aggregate", "irp"), "algorithms": algos, "tools": tools,
                 "no_full": opts.get("no_full", False), "log_path": Path(getattr(cfg, "log", "seidr.log")),
                 "force": (force or opts.get("force", False)), "verbose": True}
    _build_network_task(label="main", targeted=False, target_file=None, **task_args)