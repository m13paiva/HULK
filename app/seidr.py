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

# Force C locale to prevent std::stod errors in C++ binaries
# due to comma/dot decimal separators in different languages.
ENV_OVERRIDES = os.environ.copy()
ENV_OVERRIDES["LC_ALL"] = "C"
ENV_OVERRIDES["LC_NUMERIC"] = "C"
ENV_OVERRIDES["LANG"] = "C"

# ALGO_MAP Structure: (BinaryName, ModeFlag, ModeValue, ImportFormat)
# ImportFormat 'lm' = Lower Triangular (Symmetric)
# ImportFormat 'm'  = Dense Matrix (Asymmetric/Directed)
ALGO_MAP = {
    "PEARSON": ("correlation", "-m", "pearson", "lm"),
    "SPEARMAN": ("correlation", "-m", "spearman", "lm"),
    "PCOR": ("pcor", None, None, "lm"),
    "MI": ("mi", "-m", "RAW", "lm"),  # Reverted to RAW, marked as Lower Matrix
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


def _run_direct_visible(cmd: List[str], cwd: Path, log_path: Path):
    """
    Runs a subprocess letting stderr flow directly to the console.
    CRITICAL: Uses modified ENV_OVERRIDES to ensure LC_ALL=C.
    """
    cmd_str = " ".join(cmd)
    log(f"[EXEC] {cmd_str}", log_path)

    try:
        # We pass env=ENV_OVERRIDES to fix the 'stod' locale crashes
        subprocess.run(
            cmd,
            cwd=str(cwd),
            stderr=None,  # Inherit console stderr for progress bars
            stdout=subprocess.DEVNULL,
            env=ENV_OVERRIDES,
            check=True
        )
    except subprocess.CalledProcessError as e:
        log(f"[ERROR] Command failed with exit code {e.returncode}", log_path)
        raise e


def _import_scores(seidr_bin: str, algo_name: str, outdir: Path, prefix: str,
                   tsv_in: Path, genes: Path, fmt: str, threads: int, log_path: Path) -> Path:
    """
    Imports the TSV scores into Seidr binary format (.sf).
    """
    sf_path = outdir / f"{prefix}{algo_name.lower()}_scores.sf"

    if sf_path.exists() and sf_path.stat().st_size > 0:
        return sf_path

    cmd = [seidr_bin, "import", "-n", algo_name, "-o", str(sf_path),
           "-F", fmt, "-i", str(tsv_in), "-g", str(genes)]

    # Standardize import flags
    if algo_name in ["PEARSON", "SPEARMAN", "PCOR"]:
        cmd.extend(["-A", "-r", "-u"])
    elif algo_name == "MI":
        cmd.extend(["-r", "-u", "-O", str(threads)])
    elif algo_name in ["CLR", "ARACNE"]:
        cmd.extend(["-r", "-u", "-z", "-O", str(threads)])
    else:
        # Asymmetric methods
        cmd.extend(["-r", "-z", "-O", str(threads)])

    cmd_str = " ".join(cmd)
    log(f"[EXEC] {cmd_str}", log_path)
    try:
        # CRITICAL: env=ENV_OVERRIDES prevents locale issues during import parsing
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
    """Safely converts Seidr output strings to float, handling ';' rank info."""
    if isinstance(x, (float, int)):
        return float(x)
    if isinstance(x, str):
        try:
            # Handle '0.123;1' format from Seidr
            val = x.split(";")[0]
            return float(val)
        except (ValueError, TypeError):
            return 0.0
    return 0.0


def _export_results(outdir: Path, algorithms: List[str], seidr: str, bb_sf: Path, label: str, no_full: bool,
                    log_path: Path):
    msg = f"[Seidr] Exporting {label} results..."
    print(msg)
    log(msg, log_path)

    cmd = [seidr, "view", "--column-headers", str(bb_sf)]
    log(f">> {' '.join(cmd)} (capturing output)", log_path)

    # Use env overrides here too
    res = subprocess.run(cmd, capture_output=True, text=True, env=ENV_OVERRIDES)

    if not res.stdout.strip():
        msg = f"[Warn] No edges found in {bb_sf}. Backbone threshold too strict?"
        print(msg)
        log(msg, log_path)
        return

    try:
        # engine='python' is slower but more robust to bad lines than 'c'
        df = pd.read_csv(io.StringIO(res.stdout), sep="\t", engine="python")
    except Exception as e:
        msg = f"[Error] Failed to parse Seidr output: {e}"
        print(msg)
        log(msg, log_path)
        return

    # Clean numeric columns strictly
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
        log_path: Path
) -> None:
    prefix = f"{label}_"
    seidr = tools["seidr"]

    sf_files = []

    def run_algo_task(algo: str, prerequisite_file: Optional[Path] = None) -> Optional[Path]:
        try:
            # Unpack the 4-tuple from ALGO_MAP
            bin_name, m_flag, m_val, default_fmt = ALGO_MAP[algo]

            # Determine format: 'el' if targeted, otherwise the algo specific format ('lm' or 'm')
            current_fmt = "el" if targeted else default_fmt

            out_tsv = outdir / f"{prefix}{algo.lower()}_scores.tsv"
            done_marker = outdir / f".{prefix}{algo.lower()}.done"

            if out_tsv.exists() and done_marker.exists():
                msg = f"[Seidr] Found verified cache for {algo}. Skipping."
                print(msg)
                log(msg, log_path)
                return _import_scores(seidr, algo, outdir, prefix, out_tsv, genes_file, current_fmt, threads, log_path)

            if out_tsv.exists(): out_tsv.unlink()
            if done_marker.exists(): done_marker.unlink()

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

            # Run with direct visibility and ENV overrides
            print(f"[Seidr] Running {algo}...")
            _run_direct_visible(cmd, cwd=outdir, log_path=log_path)

            if not out_tsv.exists():
                raise FileNotFoundError(f"{algo} reported success but {out_tsv} is missing.")

            done_marker.touch()

            # Import using the correct format logic
            return _import_scores(seidr, algo, outdir, prefix, out_tsv, genes_file, current_fmt, threads, log_path)

        except Exception as e:
            msg = f"[Error] {algo} failed: {e}"
            print(msg)
            log(msg, log_path)

            done_marker_path = outdir / f".{prefix}{algo.lower()}.done"
            if done_marker_path.exists():
                done_marker_path.unlink()
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
        if max_workers > 1:
            print("\nWARNING: Workers > 1. Multiple progress bars will garble your console output.")
            print("         If you care about reading the bars, set workers: 1 in config.\n")

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {}
            for algo in parallel_algos:
                prereq = mi_tsv_path if algo in ["CLR", "ARACNE"] else None

                if prereq and not prereq.exists():
                    msg = f"[Error] Cannot run {algo}: MI output missing."
                    print(msg)
                    log(msg, log_path)
                    continue

                future = executor.submit(run_algo_task, algo, prereq)
                future_map[future] = algo

            for future in concurrent.futures.as_completed(future_map):
                res = future.result()
                if res:
                    sf_files.append(res)

    # 3. Aggregate
    if not sf_files:
        msg = f"[Error] No algorithms finished successfully for {label}."
        print(msg)
        log(msg, log_path)
        return

    net_sf = outdir / f"network_{label}.sf"
    if net_sf.exists(): net_sf.unlink()

    msg = f"[Seidr] Aggregating {len(sf_files)} networks..."
    print(msg)
    log(msg, log_path)

    agg_cmd = [seidr, "aggregate", "-o", str(net_sf), "-m", aggregate_mode] + [str(s) for s in sf_files]

    try:
        # Use subprocess manually to ensure ENV is passed (run_cmd might not support it)
        cmd_str = " ".join(agg_cmd)
        log(f"[EXEC] {cmd_str}", log_path)
        subprocess.run(agg_cmd, cwd=str(outdir), check=True, env=ENV_OVERRIDES, stdout=subprocess.DEVNULL)
    except Exception as e:
        msg = f"[Error] Aggregation failed: {e}"
        print(msg)
        log(msg, log_path)
        return

    # 4. Backbone
    msg = f"[Seidr] Pruning (Backbone {backbone})..."
    print(msg)
    log(msg, log_path)

    bb_sf = outdir / f"network_{label}.bb.sf"
    if bb_sf.exists(): bb_sf.unlink()

    try:
        cmd = [seidr, "backbone", "-F", str(backbone), str(net_sf)]
        log(f"[EXEC] {' '.join(cmd)}", log_path)
        subprocess.run(cmd, cwd=str(outdir), check=True, env=ENV_OVERRIDES, stdout=subprocess.DEVNULL)
    except Exception as e:
        msg = f"[Error] Backbone failed: {e}"
        print(msg)
        log(msg, log_path)
        return

    if bb_sf.exists():
        _export_results(outdir, algorithms, seidr, bb_sf, label, no_full, log_path)
    else:
        msg = "[Warn] Backbone file not created. Threshold might be too strict."
        print(msg)
        log(msg, log_path)


# --- MAIN RUNNER ---

def run_seidr(cfg: Config) -> None:
    log_path = Path(getattr(cfg, "log", "seidr.log"))
    opts = cfg.get_tool_opts("seidr")

    if not opts.get("enabled", False):
        return

    # 1. Locate inputs
    genes_file = Path(opts["genes_file"])
    expression_file = Path(opts["expression_file"])
    outdir = Path(opts["outdir"])
    targets = [Path(t) for t in opts.get("targets", [])]
    target_mode = opts.get("target_mode", "both")

    algos = opts.get("algorithms", [])
    if not algos:
        algos = PRESETS.get(opts.get("preset", "BALANCED"), PRESETS["BALANCED"])

    outdir.mkdir(parents=True, exist_ok=True)

    try:
        tools = _resolve_binaries(algos, None)
    except RuntimeError as e:
        msg = f"[Seidr] Error: {e}"
        print(msg)
        log(msg, log_path)
        return

    # Manual print+log for header
    sep = "=" * 60
    print(sep)
    log(sep, log_path)

    header = "STARTING SEIDR NETWORK INFERENCE"
    print(header)
    log(header, log_path)

    # Task execution
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
        "no_full": opts.get("no_full", False),
        "log_path": log_path
    }

    if not targets or target_mode in ["both", "main_only"]:
        _build_network_task(label="main", targeted=False, target_file=None, **task_args)

    if targets and target_mode in ["both", "targeted_only"]:
        for t_file in targets:
            if t_file.exists():
                _build_network_task(label=t_file.stem, targeted=True, target_file=t_file, **task_args)

    msg = "[Seidr] Analysis Finished."
    print(msg)
    log(msg, log_path)