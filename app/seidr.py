from __future__ import annotations
import os
import shutil
import sys
import io
import subprocess
import concurrent.futures
from pathlib import Path
from typing import Dict, List, Optional, Any
import re
import threading
import pandas as pd
from tqdm.auto import tqdm
from .entities import Config
from .utils import run_managed_subprocess, register_process, unregister_process, pad_desc

# ─────────────────────────────────────────────────────────────────────────────
# Environment Mapping
# ─────────────────────────────────────────────────────────────────────────────
ENV_OVERRIDES = os.environ.copy()
ENV_OVERRIDES["LC_ALL"] = "C"
ENV_OVERRIDES["LC_NUMERIC"] = "C"
ENV_OVERRIDES["LANG"] = "C"
# Prevent OpenBLAS/BLAS from spawning internal thread pools which compete or deadlock
# when called from within OpenMP loops in Seidr algorithms (e.g. GLMNET in TIGRESS, CLP in NARROMI).
ENV_OVERRIDES["OPENBLAS_NUM_THREADS"] = "1"
ENV_OVERRIDES["BLAS_NUM_THREADS"] = "1"
ENV_OVERRIDES["GOTO_NUM_THREADS"] = "1"
ENV_OVERRIDES["MKL_NUM_THREADS"] = "1"

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
    """
    Evaluates required execution binaries based on requested algorithms arrays,
    mapping localized environments.

    Args:
        needed_algos (List[str]): List of algorithms intended for execution.

    Raises:
        RuntimeError: If critical Seidr toolkit components are missing from PATH.

    Returns:
        Dict[str, str]: Absolute paths mapping component aliases.
    """
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


_seidr_pos_lock = threading.Lock()
_seidr_current_pos = 0

def _acquire_seidr_pos() -> int:
    global _seidr_current_pos
    with _seidr_pos_lock:
        pos = _seidr_current_pos
        _seidr_current_pos += 1
        return pos

def _reset_seidr_pos() -> None:
    global _seidr_current_pos
    with _seidr_pos_lock:
        _seidr_current_pos = 0

_SEIDR_PROGRESS_RE = re.compile(r"\|[=> -]+\|\s*([0-9.]+)%\s*\|\s*([^|]+)\|\s*([^|]+)\|\s*([0-9:? -]+)")

def _run_seidr_algo_stream(
    cmd: List[str],
    cwd: Path,
    env: Dict[str, str],
    log_path: Path,
    algo: str,
    threads: int
) -> None:
    """
    Executes a Seidr inference algorithm subprocess, capturing stdout/stderr in real-time.
    Progress updates (matching Seidr's progress bar format) are routed directly to a live
    Python tqdm progress bar on the console and stripped from the log file.
    Informational and diagnostic log lines are written cleanly to the log file.
    """
    cmd_str = " ".join(cmd)

    with open(log_path, "a", buffering=1) as flog:
        flog.write(f"\n[EXEC] {cmd_str}\n")
        flog.flush()

        pos = _acquire_seidr_pos()
        desc_str = pad_desc(f"[Seidr] {algo}", 18)
        pbar = tqdm(
            total=100.0,
            desc=desc_str,
            unit="%",
            position=pos,
            leave=True,
            bar_format="{desc}: {percentage:3.0f}%|{bar}| {n:.1f}% [{elapsed}{postfix}]",
            mininterval=0.1,
        )

        def _execute_process(exec_cmd: List[str], exec_env: Dict[str, str]) -> int:
            kwargs = {}
            if hasattr(os, 'setsid'):
                kwargs['start_new_session'] = True

            proc = subprocess.Popen(
                exec_cmd,
                cwd=str(cwd),
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                env=exec_env,
                text=True,
                bufsize=1,
                **kwargs
            )
            register_process(proc)

            try:
                assert proc.stdout is not None
                for raw_line in proc.stdout:
                    matches = list(_SEIDR_PROGRESS_RE.finditer(raw_line))
                    if matches:
                        last_m = matches[-1]
                        pct, rate, elapsed, eta = last_m.groups()
                        try:
                            val = min(100.0, float(pct))
                            pbar.n = val
                            pbar.set_postfix_str(f"{rate.strip()} | ETA: {eta.strip()}")
                            pbar.refresh()
                        except Exception:
                            pass

                        cleaned = _SEIDR_PROGRESS_RE.sub("", raw_line).strip()
                        if cleaned:
                            flog.write(cleaned + "\n")
                            flog.flush()
                    else:
                        flog.write(raw_line)
                        flog.flush()
            finally:
                unregister_process(proc)

            return proc.wait()

        try:
            ret = _execute_process(cmd, env)

            # Defensive fallback: If tool failed with -O, retry with -O 1
            if ret != 0 and "-O" in cmd and str(threads) != "1":
                o_idx = cmd.index("-O")
                fallback_cmd = list(cmd)
                fallback_cmd[o_idx + 1] = "1"
                fallback_env = env.copy()
                fallback_env["OMP_NUM_THREADS"] = "1"
                retry_msg = f"[Seidr] {algo} failed with -O {threads}. Retrying with single thread (-O 1)..."
                flog.write(f"\n{retry_msg}\n[EXEC] {' '.join(fallback_cmd)}\n")
                flog.flush()
                pbar.set_postfix_str("Retrying -O 1...")
                pbar.refresh()
                ret = _execute_process(fallback_cmd, fallback_env)

            if ret != 0:
                pbar.set_postfix_str("FAILED")
                pbar.refresh()
                raise subprocess.CalledProcessError(ret, cmd)
            else:
                pbar.n = 100.0
                pbar.set_postfix_str("Done")
                pbar.refresh()

        finally:
            pbar.close()


def _run_direct_quiet(cmd: List[str], cwd: Path, log_path: Path, env: Optional[Dict[str, str]] = None):
    """
    Executes a direct subprocess hook routing standard output securely against logs,
    surfacing errors globally.

    Args:
        cmd (List[str]): Executable sequence.
        cwd (Path): Expected working directory.
        log_path (Path): File destination routing output.
        env (Optional[Dict[str, str]]): Environment variables override.
    """
    cmd_str = " ".join(cmd)
    run_env = env if env is not None else ENV_OVERRIDES

    with open(log_path, "a") as f:
        f.write(f"\n[EXEC] {cmd_str}\n")
        f.flush()

        try:
            run_managed_subprocess(
                cmd,
                cwd=str(cwd),
                stdout=f,
                stderr=subprocess.STDOUT,
                env=run_env,
                check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"[Seidr ERROR] Command failed. Check {log_path.name} for details.", file=sys.stderr)
            raise e


def _import_scores(seidr_bin: str, algo_name: str, outdir: Path, prefix: str,
                   tsv_in: Path, genes: Path, fmt: str, threads: int, log_path: Path) -> Path:
    """
    Imports raw target scores mapped via algorithm execution into standard structured frames.

    Args:
        seidr_bin (str): Tool path pointing directly to standard component.
        algo_name (str): Identifier defining algorithm behavior formats.
        outdir (Path): Output mapping structure.
        prefix (str): Text grouping indicator label.
        tsv_in (Path): Input TSV file generated by prior tool run.
        genes (Path): Input node identification constraint mapping file.
        fmt (str): Expected formatting string parsing flag.
        threads (int): Subsystem thread allocation bound limits.
        log_path (Path): Expected destination for emitted messages.

    Raises:
        RuntimeError: On failure indicating malformed input components or internal toolkit failure.

    Returns:
        Path: Evaluated output `.sf` path mapping.
    """
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

    import_env = ENV_OVERRIDES.copy()
    import_env["OMP_NUM_THREADS"] = str(threads)

    try:
        _run_direct_quiet(cmd, cwd=outdir, log_path=log_path, env=import_env)
    except subprocess.CalledProcessError:
        raise RuntimeError(f"Import failed for {algo_name}")

    if not sf_path.exists():
        raise RuntimeError(f"Import reported success but {sf_path} is missing.")
    return sf_path


def _safe_float(x: Any) -> float:
    """Safe evaluation translating edge values securely."""
    try:
        if isinstance(x, (float, int)): return float(x)
        return float(str(x).split(";")[0])
    except:
        return 0.0


def _export_results(outdir: Path, algorithms: List[str], seidr: str, bb_sf: Path, label: str, no_full: bool,
                    log_path: Path):
    """
    Processes final aggregated outputs directly into flat file components natively ready for ingestion
    into downstream standard analytical R packages (like EGAD).

    Args:
        outdir (Path): Standard component location map.
        algorithms (List[str]): Input lists identifying aggregated properties.
        seidr (str): Active standard tool execution block.
        bb_sf (Path): Input target defining aggregated and backboned properties natively.
        label (str): Operational descriptor indicating output grouping string limits.
        no_full (bool): Flag defining truncation of heavy matrices output arrays limits.
        log_path (Path): Logging routing indicator output format.
    """
    print(f"[Seidr] Exporting {label} results...")

    cmd = [seidr, "view", "--column-headers", str(bb_sf)]

    try:
        proc = run_managed_subprocess(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=ENV_OVERRIDES, check=True)
        res_stdout = proc.stdout.read() if proc.stdout else ""
    except subprocess.CalledProcessError as e:
        print(f"[Seidr ERROR] Export failed: {e}", file=sys.stderr)
        return

    if not res_stdout.strip():
        print(f"[Warn] No edges found in {bb_sf}.")
        return

    try:
        df = pd.read_csv(io.StringIO(res_stdout), sep="\t", engine="python")
    except Exception as e:
        print(f"[Error] Failed to parse Seidr table: {e}")
        return

    numeric_cols = [c for c in df.columns if c not in ["Source", "Target"] and "interaction" not in c.lower()]
    for c in numeric_cols: df[c] = df[c].apply(_safe_float)

    agg_col = df.columns[-1]
    simple = df[["Source", "Target", agg_col]].copy().rename(columns={agg_col: "Weight"})

    irp_col = next((c for c in df.columns if "IRP" in c.upper()), None)
    if irp_col:
        simple["IRP_score"] = df[irp_col]

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
    """
    Subroutine managing standard analytical task sequence building dynamic network prediction components
    from available matrices parameters via mapped Seidr components.

    Args:
        outdir (Path): Directory layout for processing.
        genes_file (Path): Reference constraints file array indicator.
        expression_file (Path): Raw variance inputs matrices file indicator.
        threads (int): Concurrency sub-limits.
        max_workers (int): Processing pool worker thread bounds limitation.
        backbone (float): Discard significance floor threshold limit.
        aggregate_mode (str): Seidr aggregation sequence indicator properties block.
        algorithms (List[str]): Required calculation methods mappings limit constraints.
        tools (Dict[str, str]): Binaries environment lookup array string.
        label (str): Operational descriptor indicating output grouping string limits.
        targeted (bool): Boolean boolean indicator mapping edge output requirements parameters target array.
        target_file (Optional[Path]): Override specific target limitations string property variable.
        no_full (bool): Flag defining truncation of heavy matrices output arrays limits.
        log_path (Path): Default indicator structure limiting property mapping arrays parameters limits.
        force (bool): Toggle hard disk overwriting logic sequences indicator parameter maps.
    """
    genes_file = genes_file.resolve()
    expression_file = expression_file.resolve()
    if target_file: target_file = target_file.resolve()

    if not genes_file.exists() or not expression_file.exists():
        print(f"[Seidr ERROR] Input files missing.", file=sys.stderr)
        return

    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.touch(exist_ok=True)

    prefix = f"{label}_"
    seidr = tools["seidr"]
    _reset_seidr_pos()

    sf_files: List[Path] = []
    to_run: List[str] = []

    # 1. Pre-check existing .sf files in workflow order
    for algo in algorithms:
        out_sf = outdir / f"{prefix}{algo.lower()}_scores.sf"
        if out_sf.exists() and not force:
            print(f"[Seidr] {algo} network already exists (.sf found). Skipping calculation.", flush=True)
            sf_files.append(out_sf)
        else:
            to_run.append(algo)

    def run_algo_task(algo: str, prerequisite_file: Optional[Path] = None) -> Optional[Path]:
        try:
            if algo not in ALGO_MAP: return None
            bin_name, m_flag, m_val, default_fmt = ALGO_MAP[algo]
            current_fmt = "el" if targeted else default_fmt

            out_tsv = outdir / f"{prefix}{algo.lower()}_scores.tsv"
            out_sf = outdir / f"{prefix}{algo.lower()}_scores.sf"

            if out_sf.exists() and not force:
                return out_sf

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
                    print(f"[Seidr ERROR] {algo} missing prerequisite {prerequisite_file}.", file=sys.stderr)
                    return None
                cmd.extend(["-M", str(prerequisite_file.resolve())])

            task_env = ENV_OVERRIDES.copy()
            task_env["OPENBLAS_NUM_THREADS"] = "1"
            task_env["BLAS_NUM_THREADS"] = "1"
            task_env["GOTO_NUM_THREADS"] = "1"
            task_env["MKL_NUM_THREADS"] = "1"
            task_env["OMP_NUM_THREADS"] = str(threads)

            _run_seidr_algo_stream(cmd, outdir, task_env, log_path, algo, threads)
            return _import_scores(seidr, algo, outdir, prefix, out_tsv, genes_file, current_fmt, threads, log_path)

        except Exception as e:
            err_msg = f"[Seidr FAIL] {algo} crashed: {type(e).__name__}: {e}"
            print(err_msg, file=sys.stderr)
            with open(log_path, "a") as f:
                f.write(f"\n{err_msg}\n")
            return None

    # 2. Sequential MI if required as dependency or target
    need_mi = ("MI" in to_run) or (any(a in to_run for a in ["CLR", "ARACNE"]) and not (outdir / f"{prefix}mi_scores.tsv").exists())
    if need_mi:
        mi_sf = run_algo_task("MI")
        if mi_sf and "MI" in algorithms and mi_sf not in sf_files:
            sf_files.append(mi_sf)
        if "MI" in to_run:
            to_run.remove("MI")

    # 3. Parallel processing for remaining algorithms
    parallel_algos = [a for a in to_run if a != "MI"]
    safe_workers = int(max_workers) if max_workers else 1

    if parallel_algos:
        with concurrent.futures.ThreadPoolExecutor(max_workers=safe_workers) as executor:
            future_map = {}
            for algo in parallel_algos:
                prereq = None
                if algo in ["CLR", "ARACNE"]:
                    prereq = outdir / f"{prefix}mi_scores.tsv"
                    if not prereq.exists():
                        print(f"[Seidr ERROR] {algo} missing prerequisite {prereq.name}.", file=sys.stderr)
                        continue
                future_map[executor.submit(run_algo_task, algo, prereq)] = algo

            for future in concurrent.futures.as_completed(future_map):
                res = future.result()
                if res and res not in sf_files:
                    sf_files.append(res)

        print("", flush=True)

    if not sf_files:
        print(f"[Seidr ERROR] Aborting {label}: No algorithms succeeded.", file=sys.stderr)
        return

    # 3. Aggregate
    net_sf = outdir / f"network_{label}.sf"
    print(f"[Seidr] Aggregating {len(sf_files)} networks...")

    if force and net_sf.exists():
        net_sf.unlink()

    try:
        _run_direct_quiet([seidr, "aggregate", "-o", str(net_sf), "-m", aggregate_mode] + [str(s) for s in sf_files],
                          cwd=outdir, log_path=log_path)
    except Exception:
        return

    # 4. Backbone
    print(f"[Seidr] Backbone pruning ({backbone})...")
    bb_sf = outdir / f"network_{label}.bb.sf"

    if force and bb_sf.exists():
        bb_sf.unlink()

    out_edges = outdir / f"network_{label}_edges.tsv"

    try:
        _run_direct_quiet([seidr, "backbone", "-F", str(backbone), str(net_sf)], cwd=outdir, log_path=log_path)
        if bb_sf.exists():
            _export_results(outdir, algorithms, seidr, bb_sf, label, no_full, log_path)
    except Exception as e:
        print(f"[Warn] Backbone failed for {label}: {e}", file=sys.stderr)

    # Fallback: if backbone pruning produced no valid edge file, export aggregate network
    if not out_edges.exists() or out_edges.stat().st_size == 0:
        if net_sf.exists():
            print(f"[Seidr] Backbone yielded no edges for {label}. Exporting aggregate network as fallback...", flush=True)
            _export_results(outdir, algorithms, seidr, net_sf, label, no_full, log_path)


def run_seidr(cfg: Config, force: bool = False) -> None:
    """
    Subroutine managing execution logic boundaries extracting environment variables natively.
    """
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
    """
    Executes specific component task processing parameters directly handling saturation batch variables target strings limits mapping outputs target arrays bounds variable mapping limit properties constraint strings mapping limitation variable strings variable array limits targets variables indicator limit array targets map parameter property constraints bounds map indicator mappings property bounds properties target.
    """
    log_path = outdir / "seidr_batch.log"
    algos = [a.upper() for a in PRESETS.get(preset.upper(), PRESETS["FAST"])]
    tools = _resolve_binaries(algos)

    _build_network_task(outdir, genes_file, expression_file, (threads or cfg.max_threads), 1, 1.28, "irp", algos, tools,
                        "saturation", False, None, True, log_path, True)

    allowed_files = {"expression.tsv", "genes.txt", "network_saturation_edges.tsv", "network_edges.tsv", "seidr_batch.log", ".seidr.done"}

    for item in outdir.iterdir():
        if item.is_dir():
            shutil.rmtree(item, ignore_errors=True)
        elif item.name not in allowed_files:
            item.unlink(missing_ok=True)

    expected_edges = outdir / "network_saturation_edges.tsv"
    if expected_edges.exists() and expected_edges.stat().st_size > 0:
        (outdir / ".seidr.done").touch()
    else:
        (outdir / ".seidr.done").unlink(missing_ok=True)
        print(f"[Seidr ERROR] Batch in {outdir.name} did not yield valid edge table.", file=sys.stderr)