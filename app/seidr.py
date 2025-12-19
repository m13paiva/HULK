"""
Seidr Pipeline Integration for HULK.
"""
from __future__ import annotations
import os
import shutil
import subprocess
import sys
import io
import concurrent.futures
from pathlib import Path
from typing import Dict, List, Optional, Any
import pandas as pd

from .entities import Config  # Import for type hinting

# --- CONSTANTS ---

ALGO_MAP = {
    "PEARSON": ("correlation", "-m", "pearson"),
    "SPEARMAN": ("correlation", "-m", "spearman"),
    "PCOR": ("pcor", None, None),
    "MI": ("mi", "-m", "RAW"),
    "CLR": ("mi", "-m", "CLR"),
    "ARACNE": ("mi", "-m", "ARACNE"),
    "GENIE3": ("genie3", None, None),
    "TIGRESS": ("tigress", None, None),
    "SVM": ("svm-ensemble", "-k", "POLY"),
    "LLR": ("llr-ensemble", None, None),
    "PLSNET": ("plsnet", None, None),
    "NARROMI": ("narromi", "-m", "interior-point"),
    "ELNET": ("el-ensemble", None, None),
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

def _run_cmd(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    cmd_str = [str(c) for c in cmd]
    print(f"[Seidr] $ {' '.join(cmd_str)}")
    sys.stdout.flush()
    return subprocess.run(cmd_str, check=check)


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


def _import_scores(seidr_bin: str, algo_name: str, outdir: Path, prefix: str,
                   tsv_in: Path, genes: Path, fmt: str, threads: int) -> Path:
    sf_path = outdir / f"{prefix}{algo_name.lower()}_scores.sf"

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

    _run_cmd(cmd)
    return sf_path


def _export_results(outdir: Path, algorithms: List[str], seidr: str, bb_sf: Path, label: str, no_full: bool):
    print(f"[Seidr] Exporting {label} results...")

    res = subprocess.run([seidr, "view", "--column-headers", str(bb_sf)], capture_output=True, text=True)
    if not res.stdout.strip():
        print(f"[Warn] No edges found in {bb_sf}")
        return

    try:
        df = pd.read_csv(io.StringIO(res.stdout), sep="\t", engine="c")
    except Exception:
        return

    def clean(x):
        return float(x.split(";")[0]) if isinstance(x, str) and ";" in x else x

    numeric_cols = [c for c in df.columns if c not in ["Source", "Target"] and "interaction" not in c.lower()]
    for c in numeric_cols:
        df[c] = df[c].apply(clean)

    # Edge list
    agg_col = df.columns[-1]
    simple = df[["Source", "Target", agg_col]].copy().rename(columns={agg_col: "weight"})
    simple["direction"] = df["Interaction"] if "Interaction" in df.columns else "Undirected"
    simple.to_csv(outdir / f"network_{label}_edges.tsv", sep="\t", index=False)

    # Algorithm subset
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
        no_full: bool
) -> None:
    prefix = f"{label}_"
    seidr = tools["seidr"]
    fmt = "el" if targeted else "m"

    sf_files = []

    def run_algo_task(algo: str, prerequisite_file: Optional[Path] = None) -> Optional[Path]:
        try:
            bin_name, m_flag, m_val = ALGO_MAP[algo]
            out_tsv = outdir / f"{prefix}{algo.lower()}_scores.tsv"

            if out_tsv.exists() and out_tsv.stat().st_size > 0:
                print(f"   [Skip] {algo} output exists.")
            else:
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
                    if not prerequisite_file: raise ValueError(f"{algo} needs MI output")
                    cmd.extend(["-M", str(prerequisite_file)])

                _run_cmd(cmd)

            return _import_scores(seidr, algo, outdir, prefix, out_tsv, genes_file, fmt, threads)

        except Exception as e:
            print(f"[Error] {algo} failed: {e}")
            return None

    # 1. Dependency: MI
    mi_sf_path = None
    mi_tsv_path = outdir / f"{prefix}mi_scores.tsv"

    if any(x in algorithms for x in ["MI", "CLR", "ARACNE"]):
        print(f"[Seidr] Running MI prereq for {label}...")
        mi_sf_path = run_algo_task("MI")
        if "MI" in algorithms and mi_sf_path:
            sf_files.append(mi_sf_path)

    # 2. Parallel Execution
    parallel_algos = [a for a in algorithms if a != "MI"]

    if parallel_algos:
        print(f"[Seidr] Launching {len(parallel_algos)} algorithms in parallel...")
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {}
            for algo in parallel_algos:
                prereq = mi_tsv_path if algo in ["CLR", "ARACNE"] else None
                future = executor.submit(run_algo_task, algo, prereq)
                future_map[future] = algo

            for future in concurrent.futures.as_completed(future_map):
                res = future.result()
                if res:
                    sf_files.append(res)

    # 3. Aggregate
    if not sf_files:
        print("[Error] No algorithms finished successfully.")
        return

    net_sf = outdir / f"network_{label}.sf"
    agg_cmd = [seidr, "aggregate", "-o", str(net_sf), "-m", aggregate_mode] + [str(s) for s in sf_files]
    _run_cmd(agg_cmd)

    # 4. Backbone
    print(f"[Seidr] Pruning (Backbone {backbone})...")
    _run_cmd([seidr, "backbone", "-F", str(backbone), str(net_sf)])

    bb_sf = outdir / f"network_{label}.bb.sf"
    if bb_sf.exists():
        _export_results(outdir, algorithms, seidr, bb_sf, label, no_full)


# --- MAIN RUNNER ---

def run_seidr(cfg: Config) -> None:
    """
    Main entry point called by the pipeline. Uses settings from the Config object.
    """
    opts = cfg.get_tool_opts("seidr")

    if not opts.get("enabled", False):
        return

    # 1. Locate inputs
    genes_file = Path(opts["genes_file"])
    expression_file = Path(opts["expression_file"])
    outdir = Path(opts["outdir"])
    targets = [Path(t) for t in opts.get("targets", [])]
    target_mode = opts.get("target_mode", "both")

    # Resolve Algorithms
    algos = opts.get("algorithms", [])
    if not algos:
        algos = PRESETS.get(opts.get("preset", "BALANCED"), PRESETS["BALANCED"])

    # Basic Validation
    if not genes_file.exists():
        print(f"[Seidr] SKIPPING: Genes file not found at {genes_file}")
        print("[Seidr] Ensure post-processing ran successfully with DESeq2 enabled.")
        return
    if not expression_file.exists():
        print(f"[Seidr] SKIPPING: Expression matrix not found at {expression_file}")
        return

    # 2. Setup Output
    outdir.mkdir(parents=True, exist_ok=True)

    # 3. Resolve Binaries
    try:
        tools = _resolve_binaries(algos, None)
    except RuntimeError as e:
        print(f"[Seidr] Error: {e}")
        return

    print("\n" + "=" * 60)
    print("STARTING SEIDR NETWORK INFERENCE")
    print(f"Algorithms: {', '.join(algos)}")
    print(f"Threads: {cfg.max_threads} | Workers: {opts['workers']}")
    print("=" * 60 + "\n")

    # Shared Args
    task_args = {
        "outdir": outdir,
        "genes_file": genes_file,
        "expression_file": expression_file,
        "threads": cfg.max_threads,
        "max_workers": opts["workers"],
        "backbone": opts["backbone"],
        "aggregate_mode": opts.get("aggregate", "irp"),
        "algorithms": algos,
        "tools": tools,
        "no_full": opts.get("no_full", False)
    }

    # Run Main (Global Network)
    if not targets or target_mode in ["both", "main_only"]:
        _build_network_task(label="main", targeted=False, target_file=None, **task_args)

    # Run Targeted
    if targets and target_mode in ["both", "targeted_only"]:
        for t_file in targets:
            if t_file.exists():
                #Use filename stem directly (e.g. 'apoptosis' -> 'network_apoptosis.sf')
                label = t_file.stem
                _build_network_task(label=label, targeted=True, target_file=t_file, **task_args)
            else:
                print(f"[Seidr] Warning: Target file {t_file} not found.")

    print("\n[Seidr] Analysis Finished.\n")