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
    Locates the underlying R execution script required for EGAD analytics.

    Returns:
        Path: Absolute path to `egad.R`.

    Raises:
        FileNotFoundError: If the internal file structure has been disrupted.
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
    Delegates to the external R script to generate evaluation distributions.

    Args:
        network_file (Path): The network adjacency input.
        out_file (Path): File destination for results TSV.
        script_path (Path): Resolved path to `egad.R`.
        log_path (Path): Path for appending execution output.
        mapman_file (Path, optional): Path to MapMan mappings.
        go_file (Path, optional): Path to BioMart GO exports.
        curves_prefix (Path, optional): Base destination for storing graph interpolation data.
        do_auroc (bool): Toggle AUROC analysis.
        do_aupr (bool): Toggle AUPR analysis.

    Returns:
        Dict: Extracted performance dictionaries keyed by annotation source.
    """
    cmd = [
        "Rscript", str(script_path),
        "--network", str(network_file),
        "--output", str(out_file)
    ]

    if mapman_file: cmd.extend(["--mapman", str(mapman_file)])
    if go_file: cmd.extend(["--go_file", str(go_file)])
    if do_auroc: cmd.append("--auroc")
    if do_aupr: cmd.append("--aupr")
    if curves_prefix: cmd.extend(["--curves", str(curves_prefix)])

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
                        "macro_auc": valid_group["AUC"].mean() if "AUC" in valid_group.columns else None,
                        "macro_aupr": valid_group["AUPR"].mean() if "AUPR" in valid_group.columns else None
                    }
            else:
                results["Legacy"] = {
                    "macro_auc": df["AUC"].mean() if "AUC" in df.columns else None,
                    "macro_aupr": df["AUPR"].mean() if "AUPR" in df.columns else None
                }
        except Exception:
            return {}
    return results


def calculate_auc_from_df(df, x_col, y_col):
    """
    Calculates Area Under Curve using standard geometric summation.

    Args:
        df (DataFrame): DataFrame holding curve points.
        x_col (str): DataFrame key for x-axis points.
        y_col (str): DataFrame key for y-axis points.

    Returns:
        float | None: Geometric AUC approximation or None.
    """
    if df is None or df.empty: return None
    x = df[x_col].values
    y = df[y_col].values
    idx = np.argsort(x)
    return np.abs(np.trapezoid(y[idx], x[idx]))


def create_1x2_plot(source_data_list, show_baseline, title_prefix, out_file):
    """
    Builds standard dual-plot layout (ROC left, PR right).

    Args:
        source_data_list (list): Configuration arrays defining the curve content.
        show_baseline (bool): Toggle for rendering random chance markers.
        title_prefix (str): Text tag for the title header.
        out_file (Path): Destination output path.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    if show_baseline:
        ax1.plot([0, 1], [0, 1], color='grey', linestyle='--', label="Random (0.5000)")

    roc_plotted, prc_plotted = False, False

    for data in source_data_list:
        src = data['src']
        color = data['color']
        roc_df = data['roc_df']
        prc_df = data['prc_df']
        baseline = data['baseline']

        micro_auc = calculate_auc_from_df(roc_df, "FPR", "TPR")
        micro_aupr = calculate_auc_from_df(prc_df, "Recall", "Precision")

        if roc_df is not None and micro_auc is not None:
            roc_plotted = True
            ax1.plot(roc_df["FPR"], roc_df["TPR"], color=color, lw=2, label=f"{src} Micro (AUC: {micro_auc:.4f})")

        if prc_df is not None and micro_aupr is not None:
            prc_plotted = True
            if show_baseline:
                ax2.axhline(y=baseline, color=color, linestyle=':', alpha=0.7, label=f"{src} Random ({baseline:.4f})")
            ax2.plot(prc_df["Recall"], prc_df["Precision"], color=color, lw=2,
                     label=f"{src} Micro (AUPR: {micro_aupr:.4f})")

    if roc_plotted:
        ax1.set_title(f"{title_prefix} Micro-Averaged ROC")
        ax1.set_xlabel("False Positive Rate")
        ax1.set_ylabel("True Positive Rate")
        ax1.legend(loc="lower right", fontsize='small')
        ax1.grid(True, linestyle='--', alpha=0.6)
    else:
        ax1.text(0.5, 0.5, "ROC Data Unavailable", ha='center', va='center')

    if prc_plotted:
        ax2.set_title(f"{title_prefix} Micro-Averaged PR")
        ax2.set_xlabel("Recall")
        ax2.set_ylabel("Precision")
        ax2.legend(loc="upper right", fontsize='small')
        ax2.grid(True, linestyle='--', alpha=0.6)
    else:
        ax2.text(0.5, 0.5, "PR Data Unavailable", ha='center', va='center')

    fig.tight_layout()
    fig.savefig(out_file)
    plt.close(fig)


def create_macro_boxplots(df, out_file, colors_dict):
    """
    Renders boxplot layouts mapping term distribution profiles based on macro values.

    Args:
        df (DataFrame): The merged results DataFrame containing raw performance numbers.
        out_file (Path): Export path.
        colors_dict (dict): Dictionary mapping source types to display colors.
    """
    if df is None or df.empty: return
    sources = [s for s in ["GO", "MapMan"] if s in df["Annotation_Source"].unique()]
    if not sources: return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    auc_data = [df[(df["Annotation_Source"] == s) & (df["AUC"].notnull())]["AUC"].values for s in sources]
    aupr_data = [df[(df["Annotation_Source"] == s) & (df["AUPR"].notnull())]["AUPR"].values for s in sources]

    if any(len(d) > 0 for d in auc_data):
        bplot1 = ax1.boxplot(auc_data, labels=sources, patch_artist=True)
        for patch, src in zip(bplot1['boxes'], sources):
            patch.set_facecolor(colors_dict.get(src, "grey"))
            patch.set_alpha(0.8)
        ax1.set_title("Macro-Average Distribution (AUROC)")
        ax1.set_ylabel("Term AUROC")
        ax1.grid(True, linestyle='--', alpha=0.6)

    if any(len(d) > 0 for d in aupr_data):
        bplot2 = ax2.boxplot(aupr_data, labels=sources, patch_artist=True)
        for patch, src in zip(bplot2['boxes'], sources):
            patch.set_facecolor(colors_dict.get(src, "grey"))
            patch.set_alpha(0.8)
        ax2.set_title("Macro-Average Distribution (AUPR)")
        ax2.set_ylabel("Term AUPR")
        ax2.grid(True, linestyle='--', alpha=0.6)

    fig.tight_layout()
    fig.savefig(out_file)
    plt.close(fig)


def run_vocal_evaluation(cfg: "Config", mapman_path: Optional[Path] = None, go_file_path: Optional[Path] = None,
                         metrics: str = "both", custom_network: Optional[Path] = None):
    """
    Main execution wrapper orchestrating R execution and downstream plot rendering for analytical feedback.

    Args:
        cfg (Config): Runtime setup state.
        mapman_path (Path, optional): Explicit override for MapMan mappings.
        go_file_path (Path, optional): Explicit override for GO exports.
        metrics (str): Select which metrics to process ('both', 'auroc', 'aupr').
        custom_network (Path, optional): Network edges file override.
    """
    opts = cfg.get_tool_opts("egad")

    net_path = custom_network if custom_network else opts["network_file"]

    seidr_dir = opts["out_file"].parent
    egad_dir = seidr_dir.parent / "egad"
    egad_dir.mkdir(parents=True, exist_ok=True)

    out_path = egad_dir / opts["out_file"].name
    log_path = egad_dir / opts["log_path"].name
    curves_prefix = egad_dir / "curves"

    if not net_path.exists():
        click.secho(f"[Error] Network file not found: {net_path}", fg="red")
        return

    do_auroc = metrics in ["both", "auroc"]
    do_aupr = metrics in ["both", "aupr"]

    click.secho(f"[Evaluate] Network: {net_path.name}", fg="cyan")
    if mapman_path: click.secho(f"[Evaluate] MapMan: {mapman_path.name}", fg="cyan")
    if go_file_path: click.secho(f"[Evaluate] GO: {go_file_path.name}", fg="cyan")
    click.secho(f"[Evaluate] Generating metrics and plots at: {egad_dir}", fg="yellow")

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

        colors = {
            "MapMan": "firebrick",
            "GO": "darkorange"
        }

        source_data_list = []

        try:
            full_df = pd.read_csv(out_path, sep="\t")
        except:
            full_df = None

        for src, mets in results.items():
            res_str = []
            if do_auroc and mets["macro_auc"] is not None: res_str.append(f"Macro AUROC: {mets['macro_auc']:.4f}")
            if do_aupr and mets["macro_aupr"] is not None: res_str.append(f"Macro AUPR: {mets['macro_aupr']:.4f}")
            click.secho(f"[Results - {src}] {' | '.join(res_str)}", fg="green")

            roc_file = Path(f"{curves_prefix}_{src}_roc.tsv")
            prc_file = Path(f"{curves_prefix}_{src}_prc.tsv")
            base_file = Path(f"{curves_prefix}_{src}_baseline.tsv")

            data_pack = {
                'src': src,
                'color': colors.get(src, "purple"),
                'roc_df': pd.read_csv(roc_file, sep="\t") if roc_file.exists() else None,
                'prc_df': pd.read_csv(prc_file, sep="\t") if prc_file.exists() else None,
                'baseline': 0.0
            }

            if base_file.exists():
                b_df = pd.read_csv(base_file, sep="\t")
                if not b_df.empty and "Baseline" in b_df.columns:
                    data_pack['baseline'] = b_df["Baseline"].iloc[0]

            source_data_list.append(data_pack)

        try:
            # Generate unconstrained line graphs
            for data in source_data_list:
                src = data['src']
                create_1x2_plot([data], show_baseline=False, title_prefix=src,
                                out_file=egad_dir / f"{src}_no_baseline.pdf")
                create_1x2_plot([data], show_baseline=True, title_prefix=src,
                                out_file=egad_dir / f"{src}_with_baseline.pdf")

            # Generate cross-comparison aggregates if viable
            if len(source_data_list) > 1:
                create_1x2_plot(source_data_list, show_baseline=False, title_prefix="Combined",
                                out_file=egad_dir / "combined_no_baseline.pdf")
                create_1x2_plot(source_data_list, show_baseline=True, title_prefix="Combined",
                                out_file=egad_dir / "combined_with_baseline.pdf")

                fig_bar, ax_bar = plt.subplots(figsize=(8, 6))
                sources = list(results.keys())
                auc_vals = [results[s]["macro_auc"] if results[s]["macro_auc"] is not None else 0 for s in sources]
                aupr_vals = [results[s]["macro_aupr"] if results[s]["macro_aupr"] is not None else 0 for s in sources]

                x = np.arange(len(sources))
                bar_colors = [colors.get(s, "purple") for s in sources]

                if do_auroc and do_aupr:
                    width = 0.35
                    bars1 = ax_bar.bar(x - width / 2, auc_vals, width, label='Macro AUROC', color=bar_colors, alpha=0.8)
                    bars2 = ax_bar.bar(x + width / 2, aupr_vals, width, label='Macro AUPR', color=bar_colors,
                                       hatch='//')
                    ax_bar.bar_label(bars1, fmt='%.3f', padding=3)
                    ax_bar.bar_label(bars2, fmt='%.3f', padding=3)
                elif do_auroc:
                    width = 0.5
                    bars1 = ax_bar.bar(x, auc_vals, width, label='Macro AUROC', color=bar_colors)
                    ax_bar.bar_label(bars1, fmt='%.3f', padding=3)
                elif do_aupr:
                    width = 0.5
                    bars2 = ax_bar.bar(x, aupr_vals, width, label='Macro AUPR', color=bar_colors, hatch='//')
                    ax_bar.bar_label(bars2, fmt='%.3f', padding=3)

                ax_bar.set_ylabel('Score')
                ax_bar.set_title('Macro-Averaged Performance')
                ax_bar.set_xticks(x)
                ax_bar.set_xticklabels(sources)
                ax_bar.set_ylim(0, 1.1)

                from matplotlib.patches import Patch
                legend_elements = []
                if do_auroc: legend_elements.append(Patch(facecolor='grey', alpha=0.8, label='Macro AUROC'))
                if do_aupr: legend_elements.append(Patch(facecolor='grey', hatch='//', label='Macro AUPR'))
                ax_bar.legend(handles=legend_elements, loc='upper left')

                fig_bar.tight_layout()
                fig_bar.savefig(egad_dir / "bar_comparison.pdf")
                plt.close(fig_bar)

                create_macro_boxplots(full_df, egad_dir / "macro_boxplot_comparison.pdf", colors)

            click.secho(f"[Plots] Saved strictly formatted PDF array to {egad_dir}", fg="blue")

        except Exception as e:
            click.secho(f"[Warn] Plot generation failed. {e}", fg="yellow")
    else:
        click.secho(f"[Error] EGAD failed entirely. Check logs at {log_path}", fg="red")