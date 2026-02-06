import os
import sys
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Optional

# Fix Matplotlib cache permission issues in containers
os.environ["MPLCONFIGDIR"] = "/tmp/mpl_cache"
import matplotlib.pyplot as plt


def get_egad_script_path() -> Path:
    """
    Locates the egad.R script.
    Assumes it is located at: app/scripts/egad.R
    """
    base_path = Path(__file__).parent
    # Try 'scripts' subdirectory first
    script_path = base_path / "scripts" / "egad.R"

    if not script_path.exists():
        # Fallback to checking the same directory (if user flattened it)
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
    Wraps the R script execution.
    CAPTURES and PRINTS output to console for debugging.
    """
    cmd = [
        "Rscript", str(script_path),
        "--network", str(network_file),
        "--mapman", str(mapman_file),
        "--output", str(out_file)
    ]

    # Header for the log file
    with open(log_path, "a") as log:
        log.write(f"\n{'=' * 40}\n[EXEC] {' '.join(cmd)}\n")

    try:
        # Capture both stdout and stderr
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )

        # 1. Print to Console (Standard Error)
        # We print a header so you know which task this output belongs to
        task_id = f"[{network_file.parent.parent.name}/{network_file.parent.name}]"
        print(f"\n{task_id} R Output:\n{result.stderr}", file=sys.stderr)

        # 2. Write to Log File
        with open(log_path, "a") as log:
            log.write(result.stderr)
            log.write("[SUCCESS]\n")

    except subprocess.CalledProcessError as e:
        # ON ERROR: Print everything we captured
        task_id = f"[{network_file.parent.parent.name}/{network_file.parent.name}]"
        print(f"\n{task_id} CRASHED:\n{e.stderr}", file=sys.stderr)

        with open(log_path, "a") as log:
            log.write(f"[ERROR] Output:\n{e.stderr}\n")
            log.write(f"[FAILED] Exit Code: {e.returncode}\n")
        return None

    # Parse Result
    if out_file.exists():
        try:
            df = pd.read_csv(out_file, sep="\t")
            if df.empty: return None
            return df["AUC"].mean()
        except Exception:
            return None
    return None


def aggregate_and_plot_egad(results: List[dict], out_dir: Path):
    """
    Aggregates AUROC results and generates the saturation plot.
    """
    df = pd.DataFrame(results)
    if df.empty:
        print("[EGAD] No valid results to plot.")
        return

    # Save raw table
    raw_path = out_dir / "saturation_auroc_raw.tsv"
    df.to_csv(raw_path, sep="\t", index=False)

    # Aggregate stats
    agg = df.groupby("pct")["auc"].agg(["mean", "std", "count"]).reset_index()
    agg["se"] = agg["std"] / np.sqrt(agg["count"])

    agg_path = out_dir / "saturation_auroc_summary.tsv"
    agg.to_csv(agg_path, sep="\t", index=False)

    # Plotting
    try:
        plt.figure(figsize=(8, 6))
        plt.errorbar(
            agg["pct"], agg["mean"], yerr=agg["std"],
            fmt='-o', capsize=5, color='#2c3e50', ecolor='#e74c3c',
            linewidth=2, markersize=8, label="Mean AUROC ± Std Dev"
        )

        plt.title("Network Reconstruction Performance vs Dataset Size", fontsize=14)
        plt.xlabel("Dataset Size (% of Total Samples)", fontsize=12)
        plt.ylabel("Mean AUROC (EGAD)", fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.ylim(0.4, 1.0)

        plot_path = out_dir / "saturation_plot.pdf"
        plt.savefig(plot_path)
        plt.close()

        print(f"\n[EGAD] Results saved:")
        print(f"       Table: {agg_path}")
        print(f"       Plot:  {plot_path}")
    except Exception as e:
        print(f"[EGAD] Plotting failed: {e}")