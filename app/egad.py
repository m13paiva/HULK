import os
import sys
import subprocess
import pandas as pd
import click
from pathlib import Path
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from .entities import Config

def get_egad_script_path() -> Path:
    """
    Locates the egad.R script.
    """
    base_path = Path(__file__).parent
    script_path = base_path / "scripts" / "egad.R"
    if not script_path.exists():
        script_path = base_path / "egad.R"
    if not script_path.exists():
        raise FileNotFoundError(f"EGAD R script not found at: {script_path}")
    return script_path


def run_egad_task(
        network_file: Path,
        mapman_file: Path,
        out_file: Path,
        script_path: Path,
        log_path: Path
) -> Optional[float]:
    """
    Wraps the R script execution. Silent execution for batch/saturation processing.
    """
    cmd = [
        "Rscript", str(script_path),
        "--network", str(network_file),
        "--mapman", str(mapman_file),
        "--output", str(out_file)
    ]

    with open(log_path, "a") as log:
        log.write(f"\n{'=' * 40}\n[EXEC] {' '.join(cmd)}\n")

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )
        with open(log_path, "a") as log:
            log.write(result.stderr)
            log.write("[SUCCESS]\n")

    except subprocess.CalledProcessError as e:
        with open(log_path, "a") as log:
            log.write(f"[ERROR] Output:\n{e.stderr}\n")
            log.write(f"[FAILED] Exit Code: {e.returncode}\n")
        return None

    if out_file.exists():
        try:
            df = pd.read_csv(out_file, sep="\t")
            if df.empty: return None
            return df["AUC"].mean()
        except Exception:
            return None
    return None


def run_vocal_evaluation(cfg: "Config", mapman_path: Path):
    """
    CLI-friendly wrapper for EGAD. Always runs the evaluation and rewrites output.
    """
    opts = cfg.get_tool_opts("egad")
    net_path = opts["network_file"]
    out_path = opts["out_file"]
    log_path = opts["log_path"]

    # Basic check for input existence
    if not net_path.exists():
        click.secho(f"[Error] Network file not found: {net_path}", fg="red")
        return

    click.secho(f"[Evaluate] Network: {net_path.name}", fg="cyan")
    click.secho(f"[Evaluate] Annotation: {mapman_path.name}", fg="cyan")
    click.secho(f"[Evaluate] Rewriting results to: {out_path.name}", fg="yellow")

    # Call the core task - logic is now "always execute"
    mean_auc = run_egad_task(
        network_file=net_path,
        mapman_file=mapman_path,
        out_file=out_path,
        script_path=get_egad_script_path(),
        log_path=log_path
    )

    if mean_auc is not None:
        click.secho(f"[Success] Evaluation complete!", fg="green", bold=True)
        click.secho(f"[Results] Mean AUROC: {mean_auc:.4f}", fg="green")
        click.secho(f"[File] {out_path}", fg="blue")
    else:
        click.secho(f"[Error] EGAD failed. Check logs at {log_path}", fg="red")