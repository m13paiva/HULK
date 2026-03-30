import subprocess
import pandas as pd
import click
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Optional, Tuple, TYPE_CHECKING

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
        log_path: Path,
        curves_prefix: Optional[Path] = None,
        do_auroc: bool = True,
        do_aupr: bool = False
) -> Tuple[Optional[float], Optional[float]]:
    """
    Wraps the R script execution. Conditionally computes AUROC/AUPR based on flags.
    """
    cmd = [
        "Rscript", str(script_path),
        "--network", str(network_file),
        "--mapman", str(mapman_file),
        "--output", str(out_file)
    ]

    if do_auroc:
        cmd.append("--auroc")
    if do_aupr:
        cmd.append("--aupr")
    if curves_prefix:
        cmd.extend(["--curves", str(curves_prefix)])

    with open(log_path, "a") as log:
        log.write(f"\n{'=' * 40}\n[EXEC] {' '.join(cmd)}\n")

    try:
        result = subprocess.run(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True
        )
        with open(log_path, "a") as log:
            log.write(result.stderr)
            log.write("[SUCCESS]\n")
    except subprocess.CalledProcessError as e:
        with open(log_path, "a") as log:
            log.write(f"[ERROR] Output:\n{e.stderr}\n")
            log.write(f"[FAILED] Exit Code: {e.returncode}\n")
        return None, None

    if out_file.exists():
        try:
            df = pd.read_csv(out_file, sep="\t")
            if df.empty: return None, None

            mean_auc = df["AUC"].mean() if "AUC" in df.columns else None
            mean_aupr = df["AUPR"].mean() if "AUPR" in df.columns else None
            return mean_auc, mean_aupr

        except Exception:
            return None, None
    return None, None


def run_vocal_evaluation(cfg: "Config", mapman_path: Path, metrics: str = "both",
                         custom_network: Optional[Path] = None):
    """
    CLI-friendly wrapper for EGAD. Computes AUROC/AUPR based on user choice.
    Allows for a custom network input to override the default consensus.
    """
    opts = cfg.get_tool_opts("egad")

    # Override default network if a custom one is provided
    net_path = custom_network if custom_network else opts["network_file"]
    out_path = opts["out_file"]
    log_path = opts["log_path"]

    curves_prefix = out_path.parent / "vocal_curves"

    if not net_path.exists():
        click.secho(f"[Error] Network file not found: {net_path}", fg="red")
        return

    # Parse the metrics choice
    do_auroc = metrics in ["both", "auroc"]
    do_aupr = metrics in ["both", "aupr"]

    click.secho(f"[Evaluate] Network: {net_path.name}", fg="cyan")
    click.secho(f"[Evaluate] Annotation: {mapman_path.name}", fg="cyan")
    click.secho(f"[Evaluate] Metrics requested: {metrics.upper()}", fg="cyan")
    click.secho(f"[Evaluate] Generating metrics and curves at: {out_path.parent}", fg="yellow")

    mean_auc, mean_aupr = run_egad_task(
        network_file=net_path,
        mapman_file=mapman_path,
        out_file=out_path,
        script_path=get_egad_script_path(),
        log_path=log_path,
        curves_prefix=curves_prefix,
        do_auroc=do_auroc,
        do_aupr=do_aupr
    )

    if mean_auc is not None or mean_aupr is not None:
        click.secho(f"[Success] Evaluation complete!", fg="green", bold=True)

        res_str = []
        if do_auroc and mean_auc is not None: res_str.append(f"Mean AUROC: {mean_auc:.4f}")
        if do_aupr and mean_aupr is not None: res_str.append(f"Mean AUPR: {mean_aupr:.4f}")
        click.secho(f"[Results] {' | '.join(res_str)}", fg="green")

        try:
            # Load dataframes if they exist
            roc_df = pd.read_csv(f"{curves_prefix}_roc.tsv", sep="\t") if Path(
                f"{curves_prefix}_roc.tsv").exists() else None
            prc_df = pd.read_csv(f"{curves_prefix}_prc.tsv", sep="\t") if Path(
                f"{curves_prefix}_prc.tsv").exists() else None

            baseline_val = 0.0
            baseline_file = Path(f"{curves_prefix}_baseline.tsv")
            if baseline_file.exists():
                b_df = pd.read_csv(baseline_file, sep="\t")
                if not b_df.empty and "Baseline" in b_df.columns:
                    baseline_val = b_df["Baseline"].iloc[0]

            # --- COMBINED PLOT ---
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

            if do_auroc and mean_auc is not None and roc_df is not None:
                ax1.plot(roc_df["FPR"], roc_df["TPR"], color='darkblue', lw=2, label="Model")
                ax1.plot([0, 1], [0, 1], color='grey', linestyle='--', label="Random (AUROC: 0.5000)")
                ax1.set_title(f"ROC Curve (Mean AUROC: {mean_auc:.4f})")
                ax1.set_xlabel("False Positive Rate")
                ax1.set_ylabel("True Positive Rate")
                ax1.legend(loc="lower right")
                ax1.grid(True, linestyle='--', alpha=0.6)
            else:
                ax1.text(0.5, 0.5, "ROC Data Unavailable\n(or not requested)", ha='center', va='center')

            if do_aupr and mean_aupr is not None and prc_df is not None:
                ax2.plot(prc_df["Recall"], prc_df["Precision"], color='darkred', lw=2, label="Model")
                ax2.axhline(y=baseline_val, color='grey', linestyle='--', label=f"Random (AUPR: {baseline_val:.4f})")
                ax2.set_title(f"PR Curve (Mean AUPR: {mean_aupr:.4f})")
                ax2.set_xlabel("Recall")
                ax2.set_ylabel("Precision")
                ax2.legend(loc="upper right")
                ax2.grid(True, linestyle='--', alpha=0.6)
            else:
                ax2.text(0.5, 0.5, "PR Data Unavailable\n(or not requested)", ha='center', va='center')

            combined_plot_path = out_path.parent / "vocal_evaluation_curves.pdf"
            plt.tight_layout()
            plt.savefig(combined_plot_path)
            plt.close()

            # --- SEPARATE ROC PLOT ---
            if do_auroc and mean_auc is not None and roc_df is not None:
                fig_roc, ax_roc = plt.subplots(figsize=(6, 5))
                ax_roc.plot(roc_df["FPR"], roc_df["TPR"], color='darkblue', lw=2, label="Model")
                ax_roc.plot([0, 1], [0, 1], color='grey', linestyle='--', label="Random (AUROC: 0.5000)")
                ax_roc.set_title(f"ROC Curve (Mean AUROC: {mean_auc:.4f})")
                ax_roc.set_xlabel("False Positive Rate")
                ax_roc.set_ylabel("True Positive Rate")
                ax_roc.legend(loc="lower right")
                ax_roc.grid(True, linestyle='--', alpha=0.6)
                roc_plot_path = out_path.parent / "vocal_evaluation_roc.pdf"
                fig_roc.tight_layout()
                fig_roc.savefig(roc_plot_path)
                plt.close(fig_roc)

            # --- SEPARATE PR PLOT ---
            if do_aupr and mean_aupr is not None and prc_df is not None:
                fig_pr, ax_pr = plt.subplots(figsize=(6, 5))
                ax_pr.plot(prc_df["Recall"], prc_df["Precision"], color='darkred', lw=2, label="Model")
                ax_pr.axhline(y=baseline_val, color='grey', linestyle='--', label=f"Random (AUPR: {baseline_val:.4f})")
                ax_pr.set_title(f"PR Curve (Mean AUPR: {mean_aupr:.4f})")
                ax_pr.set_xlabel("Recall")
                ax_pr.set_ylabel("Precision")
                ax_pr.legend(loc="upper right")
                ax_pr.grid(True, linestyle='--', alpha=0.6)
                pr_plot_path = out_path.parent / "vocal_evaluation_prc.pdf"
                fig_pr.tight_layout()
                fig_pr.savefig(pr_plot_path)
                plt.close(fig_pr)

            click.secho(f"[Plots] Saved combined and individual PDFs to {out_path.parent}", fg="blue")

        except Exception as e:
            click.secho(f"[Warn] Plot generation failed. Files might be corrupt: {e}", fg="yellow")
    else:
        click.secho(f"[Error] EGAD failed entirely. Check logs at {log_path}", fg="red")