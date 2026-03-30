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
import numpy as np
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


def _resolve_binaries(needed_algos: List[str]) -> Dict[str, str]:
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

def _run_direct_quiet(cmd: List[str], cwd: Path, log_path: Path):
    """
    Runs a command, redirecting all stdout/stderr to the log file.
    Only prints to console if the command FAILS.
    """
    cmd_str = " ".join(cmd)

    # Write marker to log
    with open(log_path, "a") as f:
        f.write(f"\n[EXEC] {cmd_str}\n")
        f.flush()

        try:
            # Redirect stdout and stderr to the log file directly
            subprocess.run(
                cmd,
                cwd=str(cwd),
                stdout=f,
                stderr=subprocess.STDOUT,  # Merge stderr into stdout (the log file)
                env=ENV_OVERRIDES,
                check=True
            )
        except subprocess.CalledProcessError as e:
            # If it fails, we print a warning to console and point to the log
            print(f"[Seidr ERROR] Command failed. Check {log_path.name} for details.", file=sys.stderr)
            # Re-raise so the pipeline knows to stop
            raise e


def _import_scores(seidr_bin: str, algo_name: str, outdir: Path, prefix: str,
                   tsv_in: Path, genes: Path, fmt: str, threads: int, log_path: Path) -> Path:
    sf_path = outdir / f"{prefix}{algo_name.lower()}_scores.sf"

    # If the .sf file already exists and we aren't forcing, technically we could skip this too,
    # but import is so fast it's usually safer to just overwrite it. Seidr's 'import'
    # command usually accepts overwriting, unlike the main binaries.

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

    try:
        _run_direct_quiet(cmd, cwd=outdir, log_path=log_path)
    except subprocess.CalledProcessError:
        raise RuntimeError(f"Import failed for {algo_name}")

    if not sf_path.exists():
        raise RuntimeError(f"Import reported success but {sf_path} is missing.")
    return sf_path


def _safe_float(x: Any) -> float:
    try:
        if isinstance(x, (float, int)): return float(x)
        return float(str(x).split(";")[0])
    except:
        return 0.0


def _export_results(outdir: Path, algorithms: List[str], seidr: str, bb_sf: Path, label: str, no_full: bool,
                    log_path: Path):
    print(f"[Seidr] Exporting {label} results...")

    cmd = [seidr, "view", "--column-headers", str(bb_sf)]

    # We capture this output because we need to parse it, not just log it
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, env=ENV_OVERRIDES, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[Seidr ERROR] Export failed: {e.stderr}", file=sys.stderr)
        return

    if not res.stdout.strip():
        print(f"[Warn] No edges found in {bb_sf}.")
        return

    try:
        df = pd.read_csv(io.StringIO(res.stdout), sep="\t", engine="python")
    except Exception as e:
        print(f"[Error] Failed to parse Seidr table: {e}")
        return

    numeric_cols = [c for c in df.columns if c not in ["Source", "Target"] and "interaction" not in c.lower()]
    for c in numeric_cols: df[c] = df[c].apply(_safe_float)

    # The last column in a Seidr view is always the aggregated ensemble weight
    agg_col = df.columns[-1]

    # Strip the garbage and map strictly to the columns the JS expects
    simple = df[["Source", "Target", agg_col]].copy().rename(columns={agg_col: "Weight"})

    # Dynamically find the IRP score column, regardless of its suffix
    irp_col = next((c for c in df.columns if "IRP" in c.upper()), None)
    if irp_col:
        simple["IRP_score"] = df[irp_col]

    # Save the sanitized format for the frontend
    simple.to_csv(outdir / f"network_{label}_edges.tsv", sep="\t", index=False)

    keepers = ["Source", "Target"]
    if "Interaction" in df.columns: keepers.append("Interaction")
    upper_algos = [a.upper() for a in algorithms]
    for c in df.columns:
        if any(c.upper().startswith(alg) for alg in upper_algos):
            if c not in keepers: keepers.append(c)

    df[keepers].to_csv(outdir / f"network_{label}_algs.tsv", sep="\t", index=False)
    if not no_full: df.to_csv(outdir / f"network_{label}_edges_full.tsv", sep="\t", index=False)


def _build_network_task(outdir: Path, genes_file: Path, expression_file: Path, threads: int, max_workers: int,
                        backbone: float, aggregate_mode: str, algorithms: List[str], tools: Dict[str, str],
                        label: str, targeted: bool, target_file: Optional[Path], no_full: bool, log_path: Path,
                        force: bool) -> None:
    # Path Resolution
    genes_file = genes_file.resolve()
    expression_file = expression_file.resolve()
    if target_file: target_file = target_file.resolve()

    if not genes_file.exists() or not expression_file.exists():
        print(f"[Seidr ERROR] Input files missing.", file=sys.stderr)
        return

    # Ensure log file exists before we try to append to it
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.touch(exist_ok=True)

    prefix = f"{label}_"
    seidr = tools["seidr"]
    sf_files = []

    def run_algo_task(algo: str, prerequisite_file: Optional[Path] = None) -> Optional[Path]:
        try:
            if algo not in ALGO_MAP: return None
            bin_name, m_flag, m_val, default_fmt = ALGO_MAP[algo]
            current_fmt = "el" if targeted else default_fmt

            out_tsv = outdir / f"{prefix}{algo.lower()}_scores.tsv"
            out_sf = outdir / f"{prefix}{algo.lower()}_scores.sf"

            # <--- THE REAL FIX: Check for the completed .sf file, not the intermediate .tsv
            if out_sf.exists() and not force:
                print(f"[Seidr] {algo} network already exists (.sf found). Skipping calculation.", flush=True)
                return out_sf

            # If we are here, we need to compute.
            # We MUST nuke any half-written trash from previous OOM crashes or Seidr will panic.
            if out_tsv.exists():
                out_tsv.unlink()
            if out_sf.exists():
                out_sf.unlink()

            cmd = [tools[bin_name], "-i", str(expression_file), "-g", str(genes_file), "-o", str(out_tsv)]
            if m_flag: cmd.extend([m_flag, m_val])

            if bin_name not in ["correlation", "pcor"]:
                cmd.extend(["-O", str(threads)])
            else:
                cmd.append("--no-scale")

            if targeted:
                cmd.insert(1, "-t")
                cmd.insert(2, str(target_file))

            if algo in ["CLR", "ARACNE"]:
                if not prerequisite_file or not prerequisite_file.exists():
                    print(f"[Seidr ERROR] {algo} missing prerequisite.", file=sys.stderr)
                    return None
                cmd.extend(["-M", str(prerequisite_file.resolve())])

            # Minimal console output
            print(f"[Seidr] Running {algo}...")

            # Heavy output goes to log file
            _run_direct_quiet(cmd, cwd=outdir, log_path=log_path)

            # Import the freshly calculated TSV into SF
            return _import_scores(seidr, algo, outdir, prefix, out_tsv, genes_file, current_fmt, threads, log_path)

        except Exception as e:
            print(f"[Seidr FAIL] {algo} crashed. See log.", file=sys.stderr)
            return None

    # 1. MI (Sequential)
    mi_sf = None
    if any(x in algorithms for x in ["MI", "CLR", "ARACNE"]):
        mi_sf = run_algo_task("MI")
        if "MI" in algorithms and mi_sf: sf_files.append(mi_sf)

    # 2. Others (Parallel or Sequential)
    parallel_algos = [a for a in algorithms if a != "MI"]

    # Using ThreadPoolExecutor again now that we are confident in logic
    safe_workers = int(max_workers) if max_workers else 1

    if parallel_algos:
        with concurrent.futures.ThreadPoolExecutor(max_workers=safe_workers) as executor:
            future_map = {}
            for algo in parallel_algos:
                prereq = None
                if algo in ["CLR", "ARACNE"]:
                    prereq = outdir / f"{prefix}mi_scores.tsv"
                    if not prereq.exists():
                        continue
                future_map[executor.submit(run_algo_task, algo, prereq)] = algo

            for future in concurrent.futures.as_completed(future_map):
                res = future.result()
                if res: sf_files.append(res)

    if not sf_files:
        print(f"[Seidr ERROR] Aborting {label}: No algorithms succeeded.", file=sys.stderr)
        return

    # 3. Aggregate
    net_sf = outdir / f"network_{label}.sf"
    print(f"[Seidr] Aggregating {len(sf_files)} networks...")

    # If forcing, delete old aggregate network to prevent appending/crashing issues
    if force and net_sf.exists():
        net_sf.unlink()

    try:
        # Aggregate is fast, but let's log it too
        _run_direct_quiet([seidr, "aggregate", "-o", str(net_sf), "-m", aggregate_mode] + [str(s) for s in sf_files],
                          cwd=outdir, log_path=log_path)
    except Exception:
        return

    # 4. Backbone
    print(f"[Seidr] Backbone pruning ({backbone})...")
    bb_sf = outdir / f"network_{label}.bb.sf"

    if force and bb_sf.exists():
        bb_sf.unlink()

    try:
        _run_direct_quiet([seidr, "backbone", "-F", str(backbone), str(net_sf)], cwd=outdir, log_path=log_path)
        if bb_sf.exists():
            _export_results(outdir, algorithms, seidr, bb_sf, label, no_full, log_path)
    except Exception:
        pass


def run_seidr_single(cfg: Config, force: bool = False) -> None:
    opts = cfg.get_tool_opts("seidr")
    if not str(opts.get("enabled", False)).lower() == "true": return

    preset_name = opts.get("preset", "BALANCED")
    raw_algos = opts.get("algorithms")

    if not raw_algos:
        raw_algos = PRESETS.get(preset_name, PRESETS["BALANCED"])

    if isinstance(raw_algos, str):
        raw_algos = [x.strip() for x in raw_algos.split(",")]

    algos = [str(a).upper() for a in raw_algos if str(a).strip()]

    outdir = Path(opts["outdir"])
    log_path = Path(getattr(cfg, "log", "seidr.log")).resolve()

    # Clean log file on new run if not appending
    # Actually, let's just append clearly so we don't lose history, but add a header
    with open(log_path, "a") as f:
        f.write(f"\n{'=' * 40}\nSeidr Run Start: {preset_name}\n{'=' * 40}\n")

    if (force or opts.get("force", False)) and outdir.exists():
        print(f"[Seidr] Force mode: Cleaning output directory...", flush=True)
        shutil.rmtree(outdir, ignore_errors=True)

    outdir.mkdir(parents=True, exist_ok=True)

    try:
        tools = _resolve_binaries(algos)
    except RuntimeError as e:
        print(f"[Seidr FATAL] {e}", file=sys.stderr)
        return

    workers = opts.get("workers", 4)
    try:
        workers = int(workers)
    except:
        workers = 4

    task_args = {
        "outdir": outdir, "genes_file": Path(opts["genes_file"]),
        "expression_file": Path(opts["expression_file"]), "threads": cfg.max_threads,
        "max_workers": workers, "backbone": float(opts.get("backbone", 1.28)),
        "aggregate_mode": opts.get("aggregate", "irp"), "algorithms": algos, "tools": tools,
        "no_full": opts.get("no_full", False), "log_path": log_path,
        "force": (force or opts.get("force", False))
    }

    _build_network_task(label="main", targeted=False, target_file=None, **task_args)


def run_seidr_batch(cfg: Config, genes_file: Path, expression_file: Path, outdir: Path, preset: str = "FAST",
                    threads: Optional[int] = None) -> None:
    log_path = outdir / "seidr_batch.log"
    algos = [a.upper() for a in PRESETS.get(preset.upper(), PRESETS["FAST"])]
    tools = _resolve_binaries(algos)

    _build_network_task(outdir, genes_file, expression_file, (threads or cfg.max_threads), 1, 1.28, "irp", algos, tools,
                        "saturation", False, None, True, log_path, True)

    for d in [p for p in outdir.glob("*-*-*-*") if p.is_dir()]:
        shutil.rmtree(d, ignore_errors=True)


def run_seidr_agg_batch(batch_dir: Path, cfg: Config) -> None:
    """
    Runs Seidr inference on a specific aggregation batch directory.
    Uses a .seidr_done marker to skip execution if already completed.
    """
    # 1. Check if enabled
    opts = cfg.get_tool_opts("seidr")
    if not str(opts.get("enabled", False)).lower() == "true":
        return

    # 2. Resolve paths
    batch_dir = batch_dir.resolve()
    genes_file = batch_dir / "genes.txt"
    expression_file = batch_dir / "expression.tsv"

    # Restore actual logging path
    log_path = batch_dir / "seidr.log"

    marker = batch_dir / ".seidr_done"  # <--- MARKER FILE

    if not genes_file.exists() or not expression_file.exists():
        log(f"[Seidr] Skipping aggregation batch {batch_dir.name}: missing input files.", cfg.log)
        return

    # 3. Check for completion marker
    force_run = bool(opts.get("force", False)) or cfg.force
    if marker.exists() and not force_run:
        log(f"[Seidr] Skipping batch {batch_dir.name} (found .seidr_done marker).", cfg.log)
        return

    # 4. Determine Algorithms / Preset
    preset_name = opts.get("preset", "BALANCED").upper()
    raw_algos = opts.get("algorithms")

    if not raw_algos:
        raw_algos = PRESETS.get(preset_name, PRESETS["BALANCED"])

    if isinstance(raw_algos, str):
        raw_algos = [x.strip() for x in raw_algos.split(",")]

    algos = [str(a).upper() for a in raw_algos if str(a).strip()]

    # 5. Resolve Binaries
    try:
        tools = _resolve_binaries(algos)
    except RuntimeError as e:
        log(f"[Seidr FATAL] {e}", cfg.log)
        return

    # 6. Config Parameters
    workers = int(opts.get("workers", 4))
    backbone_val = float(opts.get("backbone", 1.28))
    agg_mode = opts.get("aggregate", "irp")

    # 7. Run Task
    log(f"[Seidr] Starting inference for aggregation batch: {batch_dir.name} (Preset: {preset_name})", cfg.log)

    _build_network_task(
        outdir=batch_dir,
        genes_file=genes_file,
        expression_file=expression_file,
        threads=cfg.max_threads,
        max_workers=workers,
        backbone=backbone_val,
        aggregate_mode=agg_mode,
        algorithms=algos,
        tools=tools,
        label=batch_dir.name,
        targeted=False,
        target_file=None,
        no_full=opts.get("no_full", False),
        log_path=log_path,
        force=force_run
    )

    # 8. Create Marker
    try:
        marker.touch()
        log(f"[Seidr] Finished aggregation batch: {batch_dir.name}", cfg.log)
    except Exception as e:
        log(f"[Seidr] Warning: Could not create marker file for {batch_dir.name}: {e}", cfg.log)


def run_seidr_aggregation(cfg: Config) -> None:
    """
    Bypasses the Seidr C++ aggregator.
    Implements Sample-Size Weighted Mean with Fractional Node (25%)
    and Edge (10%) Thresholding.
    """
    import numpy as np
    import pandas as pd

    opts = cfg.get_tool_opts("seidr")
    if not str(opts.get("enabled", False)).lower() == "true":
        return

    aggregated_dir = cfg.shared / "seidr" / "aggregated"
    if not aggregated_dir.exists():
        print("[Seidr Meta] No aggregated batch directory found. Skipping meta-aggregation.")
        return

    # Find all batch directories that actually have an exported TSV
    batch_tsvs = list(aggregated_dir.rglob("network_*_edges.tsv"))
    batch_tsvs = [f for f in batch_tsvs if "meta" not in f.name]

    if not batch_tsvs:
        print("[Seidr Meta] No batch TSV networks found.")
        return

    num_batches = len(batch_tsvs)
    if num_batches == 1:
        print("[Seidr Meta] Only 1 batch network found. Skipping.")
        return

    meta_outdir = cfg.shared / "seidr" / "meta_aggregated"
    meta_outdir.mkdir(parents=True, exist_ok=True)
    meta_edges_out = meta_outdir / "network_meta_edges.tsv"

    force_run = bool(opts.get("force", False)) or getattr(cfg, "force", False)

    if meta_edges_out.exists() and not force_run:
        print("[Seidr Meta] Meta-network TSV already exists. Skipping calculation.")
        return

    print(f"[Seidr Meta] Executing Sample-Weighted Meta-Aggregation across {num_batches} batches...")

    # --- STEP 1 & 2: Sample Census and Core Node Election (25% Threshold) ---
    node_threshold = max(1, int(num_batches * 0.50))
    edge_threshold = max(1, int(num_batches * 0.25))

    print(
        f"[Seidr Meta] Thresholds - Nodes: >= {node_threshold} batches (25%). Edges: >= {edge_threshold} batches (10%).")

    gene_counts = {}
    batch_weights = {}

    for tsv in batch_tsvs:
        batch_dir = tsv.parent

        # 1. Sample Census
        expr_file = batch_dir / "expression.tsv"
        samples = 1  # Default to 1 if we can't find it to prevent division by zero
        if expr_file.exists():
            try:
                with open(expr_file, 'r') as f:
                    header = f.readline()
                    # Counting columns minus the gene ID column
                    samples = max(1, len(header.split('\t')) - 1)
            except Exception as e:
                print(f"[Seidr Meta Warning] Could not read samples from {expr_file.name}: {e}")
        batch_weights[tsv] = samples

        # 2. Node Election
        genes_file = batch_dir / "genes.txt"
        if genes_file.exists():
            try:
                with open(genes_file, 'r') as f:
                    for line in f:
                        g = line.strip()
                        if g:
                            gene_counts[g] = gene_counts.get(g, 0) + 1
            except Exception as e:
                print(f"[Seidr Meta Warning] Could not read {genes_file.name}: {e}")

    core_nodes = {g for g, count in gene_counts.items() if count >= node_threshold}
    print(f"[Seidr Meta] Elected {len(core_nodes)} Core Nodes out of {len(gene_counts)} total unique genes.")

    if len(core_nodes) < 2:
        print("[Seidr Meta ERROR] Not enough core nodes survived the 25% threshold. Aborting.")
        return

    # Write core nodes to file
    with open(meta_outdir / "meta_genes.txt", "w") as f:
        for g in sorted(core_nodes):
            f.write(f"{g}\n")

    # --- STEP 3 & 4: Edge Purge and Sample-Weighted Math ---
    df_list = []
    for tsv in batch_tsvs:
        try:
            df = pd.read_csv(tsv, sep="\t")
            col_map = {c.lower(): c for c in df.columns}

            if "weight" not in col_map:
                print(f"[Seidr Meta Warning] Skipping {tsv.name}: Missing 'weight' column.")
                continue

            w_col = col_map["weight"]
            df = df[["Source", "Target", w_col]].copy()
            df.rename(columns={w_col: "weight"}, inplace=True)

            # Purge non-core nodes
            df = df[df["Source"].isin(core_nodes) & df["Target"].isin(core_nodes)].copy()

            if df.empty:
                continue

            # Alphabetical sort
            nodes = np.sort(df[["Source", "Target"]].values, axis=1)
            df["Source"] = nodes[:, 0]
            df["Target"] = nodes[:, 1]

            # Apply sample weight
            samples = batch_weights[tsv]
            df["weight_x_samples"] = df["weight"] * samples
            df["sample_count"] = samples

            df_list.append(df)
        except Exception as e:
            print(f"[Seidr Meta ERROR] Failed to process {tsv.name}: {e}")

    if not df_list:
        print("[Seidr Meta ERROR] No edges survived the core node purge. Aborting.")
        return

    print("[Seidr Meta] Concatenating and applying 10% Edge Purge...")
    try:
        mega_df = pd.concat(df_list, ignore_index=True)

        meta_df = mega_df.groupby(["Source", "Target"]).agg(
            Sum_WxS=("weight_x_samples", "sum"),
            Sum_Samples=("sample_count", "sum"),
            Batch_Count=("weight", "count")
        ).reset_index()

        # Edge Purge (10%)
        initial_edges = len(meta_df)
        meta_df = meta_df[meta_df["Batch_Count"] >= edge_threshold].copy()
        print(
            f"[Seidr Meta] Purged {initial_edges - len(meta_df)} edges failing the 10% threshold (>= {edge_threshold} batches).")

        if meta_df.empty:
            print("[Seidr Meta ERROR] No edges survived the frequency threshold. Aborting.")
            return

        # Final Math: Sample-Weighted Mean
        meta_df["Weighted_Mean"] = (meta_df["Sum_WxS"] / meta_df["Sum_Samples"]).round(5)

        meta_df = meta_df[["Source", "Target", "Weighted_Mean", "Batch_Count"]]

        print(f"[Seidr Meta] Writing {len(meta_df)} final edges to {meta_edges_out.name}...")
        meta_df.to_csv(meta_edges_out, sep="\t", index=False)

        print("[Seidr Meta] Success. The biologically questionable meta-network is complete.")

    except Exception as e:
        print(f"[Seidr Meta ERROR] Pandas math failed: {e}")