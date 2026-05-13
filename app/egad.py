import subprocess
import pandas as pd
import click
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import Optional, Tuple, TYPE_CHECKING, Dict, Any

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
        out_file: Path,
        script_path: Path,
        log_path: Path,
        mapman_file: Optional[Path] = None,
        go_file: Optional[Path] = None,
        curves_prefix: Optional[Path] = None,
        do_auroc: bool = True,
        do_aupr: bool = False
) -> Dict[str, Dict[str, Optional[float]]]:
    """
    Wraps the R script execution. Dynamically handles multiple annotation sources.
    Returns a dictionary mapping Annotation_Source to its AUC/AUPR scores.
    """
    cmd = [
        "Rscript", str(script_path),
        "--network", str(network_file),
        "--output", str(out_file)
    ]

    if mapman_file:
        cmd.extend(["--mapman", str(mapman_file)])
    if go_file:
        cmd.extend(["--go_file", str(go_file)])

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
        return {}

    results = {}
    if out_file.exists():
        try:
            df = pd.read_csv(out_file, sep="\t")
            if df.empty: return {}

            if "Annotation_Source" in df.columns:
                for source, group in df.groupby("Annotation_Source"):
                    valid_group = group[group["Term"] != "None"]
                    if valid_group.empty: continue

                    results[source] = {
                        "auc": valid_group["AUC"].mean() if "AUC" in valid_group.columns else None,
                        "aupr": valid_group["AUPR"].mean() if "AUPR" in valid_group.columns else None
                    }
            else:
                results["Legacy"] = {
                    "auc": df["AUC"].mean() if "AUC" in df.columns else None,
                    "aupr": df["AUPR"].mean() if "AUPR" in df.columns else None
                }
        except Exception:
            return {}

    return results


def run_vocal_evaluation(cfg: "Config", mapman_path: Optional[Path] = None, go_file_path: Optional[Path] = None,
                         metrics: str = "both", custom_network: Optional[Path] = None):
    """
    CLI-friendly wrapper for EGAD. Computes AUROC/AUPR based on user choice.
    Allows for a custom network input to override the default consensus.
    """
    opts = cfg.get_tool_opts("egad")

    net_path = custom_network if custom_network else opts["network_file"]
    out_path = opts["out_file"]
    log_path = opts["log_path"]

    curves_prefix = out_path.parent / "vocal_curves"

    if not net_path.exists():
        click.secho(f"[Error] Network file not found: {net_path}", fg="red")
        return

    do_auroc = metrics in ["both", "auroc"]
    do_aupr = metrics in ["both", "aupr"]

    click.secho(f"[Evaluate] Network: {net_path.name}", fg="cyan")
    if mapman_path: click.secho(f"[Evaluate] MapMan: {mapman_path.name}", fg="cyan")
    if go_file_path: click.secho(f"[Evaluate] GO: {go_file_path.name}", fg="cyan")
    click.secho(f"[Evaluate] Metrics requested: {metrics.upper()}", fg="cyan")
    click.secho(f"[Evaluate] Generating metrics and curves at: {out_path.parent}", fg="yellow")

    results = run_egad_task(
        network_file=net_path,
        mapman_file=mapman_path,
        go_file=go_file_path,
        out_file=out_path,
        script_path=get_egad_script_path(),
        log_path=log_path,
        curves_prefix=curves_prefix,
        do_auroc=do_auroc,
        do_aupr=do_aupr
    )

    if results:
        click.secho(f"\n[Success] Evaluation complete!", fg="green", bold=True)

        for src, mets in results.items():
            res_str = []
            if do_auroc and mets["auc"] is not None: res_str.append(f"Mean AUROC: {mets['auc']:.4f}")
            if do_aupr and mets["aupr"] is not None: res_str.append(f"Mean AUPR: {mets['aupr']:.4f}")
            click.secho(f"[Results - {src}] {' | '.join(res_str)}", fg="green")

        try:
            colors = {
                "MapMan": {"roc": "darkblue", "prc": "darkred", "bar_auc": "royalblue", "bar_prc": "firebrick"},
                "GO": {"roc": "teal", "prc": "darkorange", "bar_auc": "lightseagreen", "bar_prc": "orange"}
            }
            default_colors = {"roc": "purple", "prc": "brown", "bar_auc": "mediumpurple", "bar_prc": "sienna"}

            # --- LOOP TO GENERATE BOTH WITH AND WITHOUT BASELINE ---
            for show_baseline in [True, False]:
                suffix = "_with_baseline" if show_baseline else "_no_baseline"

                fig_comb, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
                fig_roc, ax_roc = plt.subplots(figsize=(6, 5))
                fig_pr, ax_pr = plt.subplots(figsize=(6, 5))

                if show_baseline:
                    ax1.plot([0, 1], [0, 1], color='grey', linestyle='--', label="Random (AUROC: 0.5000)")
                    ax_roc.plot([0, 1], [0, 1], color='grey', linestyle='--', label="Random (AUROC: 0.5000)")

                roc_plotted = False
                prc_plotted = False

                for src, mets in results.items():
                    src_color = colors.get(src, default_colors)
                    mean_auc = mets["auc"]
                    mean_aupr = mets["aupr"]

                    roc_file = Path(f"{curves_prefix}_{src}_roc.tsv")
                    prc_file = Path(f"{curves_prefix}_{src}_prc.tsv")
                    base_file = Path(f"{curves_prefix}_{src}_baseline.tsv")

                    roc_df = pd.read_csv(roc_file, sep="\t") if roc_file.exists() else None
                    prc_df = pd.read_csv(prc_file, sep="\t") if prc_file.exists() else None

                    baseline_val = 0.0
                    if base_file.exists():
                        b_df = pd.read_csv(base_file, sep="\t")
                        if not b_df.empty and "Baseline" in b_df.columns:
                            baseline_val = b_df["Baseline"].iloc[0]

                    if do_auroc and mean_auc is not None and roc_df is not None:
                        roc_plotted = True
                        label = f"{src} Model"
                        ax1.plot(roc_df["FPR"], roc_df["TPR"], color=src_color["roc"], lw=2, label=label)
                        ax_roc.plot(roc_df["FPR"], roc_df["TPR"], color=src_color["roc"], lw=2, label=label)

                    if do_aupr and mean_aupr is not None and prc_df is not None:
                        prc_plotted = True
                        label = f"{src} Model"

                        if show_baseline:
                            base_label = f"{src} Random ({baseline_val:.4f})"
                            ax2.axhline(y=baseline_val, color=src_color["roc"], linestyle=':', label=base_label,
                                        alpha=0.7)
                            ax_pr.axhline(y=baseline_val, color=src_color["roc"], linestyle=':', label=base_label,
                                          alpha=0.7)

                        ax2.plot(prc_df["Recall"], prc_df["Precision"], color=src_color["prc"], lw=2, label=label)
                        ax_pr.plot(prc_df["Recall"], prc_df["Precision"], color=src_color["prc"], lw=2, label=label)

                # Finalize plots for this baseline iteration
                if roc_plotted:
                    for ax in [ax1, ax_roc]:
                        ax.set_title("Combined ROC Curves")
                        ax.set_xlabel("False Positive Rate")
                        ax.set_ylabel("True Positive Rate")
                        ax.legend(loc="lower right")
                        ax.grid(True, linestyle='--', alpha=0.6)
                    fig_roc.tight_layout()
                    fig_roc.savefig(out_path.parent / f"vocal_evaluation_roc{suffix}.pdf")
                else:
                    ax1.text(0.5, 0.5, "ROC Data Unavailable", ha='center', va='center')

                if prc_plotted:
                    for ax in [ax2, ax_pr]:
                        ax.set_title("Combined PR Curves")
                        ax.set_xlabel("Recall")
                        ax.set_ylabel("Precision")
                        ax.legend(loc="upper right", fontsize='small')
                        ax.grid(True, linestyle='--', alpha=0.6)
                    fig_pr.tight_layout()
                    fig_pr.savefig(out_path.parent / f"vocal_evaluation_prc{suffix}.pdf")
                else:
                    ax2.text(0.5, 0.5, "PR Data Unavailable", ha='center', va='center')

                fig_comb.tight_layout()
                fig_comb.savefig(out_path.parent / f"vocal_evaluation_curves{suffix}.pdf")

                plt.close(fig_comb)
                plt.close(fig_roc)
                plt.close(fig_pr)

            # --- COMPARISON BAR CHART (Generated once) ---
            if len(results) > 1:
                fig_bar, ax_bar = plt.subplots(figsize=(8, 6))

                sources = list(results.keys())
                auc_vals = [results[s]["auc"] if results[s]["auc"] is not None else 0 for s in sources]
                aupr_vals = [results[s]["aupr"] if results[s]["aupr"] is not None else 0 for s in sources]

                x = np.arange(len(sources))

                if do_auroc and do_aupr:
                    width = 0.35
                    bars1 = ax_bar.bar(x - width / 2, auc_vals, width, label='Mean AUROC',
                                       color=[colors.get(s, default_colors)["bar_auc"] for s in sources])
                    bars2 = ax_bar.bar(x + width / 2, aupr_vals, width, label='Mean AUPR',
                                       color=[colors.get(s, default_colors)["bar_prc"] for s in sources])
                    ax_bar.bar_label(bars1, fmt='%.3f', padding=3)
                    ax_bar.bar_label(bars2, fmt='%.3f', padding=3)
                elif do_auroc:
                    width = 0.5
                    bars1 = ax_bar.bar(x, auc_vals, width, label='Mean AUROC',
                                       color=[colors.get(s, default_colors)["bar_auc"] for s in sources])
                    ax_bar.bar_label(bars1, fmt='%.3f', padding=3)
                elif do_aupr:
                    width = 0.5
                    bars2 = ax_bar.bar(x, aupr_vals, width, label='Mean AUPR',
                                       color=[colors.get(s, default_colors)["bar_prc"] for s in sources])
                    ax_bar.bar_label(bars2, fmt='%.3f', padding=3)

                ax_bar.set_ylabel('Score')
                ax_bar.set_title('EGAD Performance Comparison')
                ax_bar.set_xticks(x)
                ax_bar.set_xticklabels(sources)
                ax_bar.set_ylim(0, 1.1)

                from matplotlib.patches import Patch
                legend_elements = []
                if do_auroc: legend_elements.append(Patch(facecolor='grey', label='Mean AUROC'))
                if do_aupr: legend_elements.append(Patch(facecolor='black', label='Mean AUPR'))
                ax_bar.legend(handles=legend_elements, loc='upper left')

                fig_bar.tight_layout()
                fig_bar.savefig(out_path.parent / "vocal_evaluation_comparison_bars.pdf")
                plt.close(fig_bar)

            click.secho(f"[Plots] Saved ALL plot variations to {out_path.parent}", fg="blue")

        except Exception as e:
            click.secho(f"[Warn] Plot generation failed. Files might be corrupt: {e}", fg="yellow")
    else:
        click.secho(f"[Error] EGAD failed entirely or returned no valid data. Check logs at {log_path}", fg="red")