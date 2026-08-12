import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, List, TYPE_CHECKING
import concurrent.futures
from tqdm import tqdm
import random
import os
import shutil
import json

os.environ["MPLCONFIGDIR"] = "/tmp"
import matplotlib.pyplot as plt

from .seidr import run_seidr_batch
from .egad import run_egad_task, get_egad_script_path
from .utils import setup_interrupt_handlers, kill_all_active_processes

if TYPE_CHECKING:
    from .entities import Dataset, Config


class BatchOrchestrator:
    """
    Manages the overall execution state processing saturation loops and mapping batch subset intervals
    to dynamically evaluate network prediction performance bounds.
    """

    def __init__(self, dataset: "Dataset", config: "Config",
                 seed: Optional[int] = None,
                 workers: int = 4,
                 max_threads: Optional[int] = None,
                 mapman_file: Optional[Path] = None,
                 go_file: Optional[Path] = None,
                 force: bool = False,
                 num_steps: int = 10,
                 metrics: str = "both",
                 target_genes_files: Optional[List[Path]] = None):

        self.dataset = dataset
        self.config = config
        self.seed = seed
        self.workers = workers
        self.mapman_file = mapman_file
        self.go_file = go_file
        self.force = force
        self.num_steps = num_steps

        if target_genes_files is not None:
            self.target_genes_files = target_genes_files
        else:
            self.target_genes_files = getattr(self.config, "target_genes_files", []) or []

        self.do_auroc = metrics.lower() in ["both", "auroc"]
        self.do_aupr = metrics.lower() in ["both", "aupr"]

        if max_threads is not None:
            self.total_thread_budget = max_threads
        else:
            self.total_thread_budget = self.config.max_threads

        if self.seed is not None:
            np.random.seed(self.seed)
            random.seed(self.seed)
            print(f"[Saturation] Random seed: {self.seed}")

        self.base_outdir = self.config.shared / "saturation"
        self.base_outdir.mkdir(parents=True, exist_ok=True)

        self.iterations = 10
        self.deseq2_dir = self.config.shared / "deseq2"
        self.vst_path = self.deseq2_dir / "vst.tsv"
        self.master_df = None

    def _get_master_df(self) -> pd.DataFrame:
        if self.master_df is not None: return self.master_df
        try:
            df = pd.read_csv(self.vst_path, sep="\t", index_col=0, engine="pyarrow")
            self.master_df = df
            return self.master_df
        except:
            df = pd.read_csv(self.vst_path, sep="\t", index_col=0)
            self.master_df = df
            return df

    def _write_seidr_files(self, df: pd.DataFrame, out_dir: Path):
        df = df.fillna(0)
        g_file, e_file = out_dir / "genes.txt", out_dir / "expression.tsv"
        with open(g_file, "w") as f:
            f.write("\n".join(df.columns.astype(str)))
        df.to_csv(e_file, sep="\t", header=False, index=False, float_format='%.6f')
        return g_file, e_file

    def _count_samples_from_file(self, expr_file: Path) -> int:
        if not expr_file.exists(): return 0
        try:
            with open(expr_file, 'rb') as f:
                return sum(1 for _ in f)
        except:
            return 0

    def _calculate_steps(self, total_bps: int) -> List[int]:
        if total_bps <= self.num_steps: return list(range(1, total_bps + 1))
        base = total_bps // self.num_steps
        remainder = total_bps % self.num_steps
        increments = [base] * self.num_steps
        for i in range(1, remainder + 1): increments[-i] += 1
        steps, current = [], 0
        for inc in increments:
            current += inc
            steps.append(current)
        return sorted(list(set(steps)))

    def _read_target_genes(self, target_file: Path) -> set:
        if not target_file.exists():
            return set()
        genes = set()
        with open(target_file, "r") as f:
            for line in f:
                line_str = line.strip()
                if line_str and not line_str.startswith("#"):
                    genes.add(line_str)
        return genes

    def _filter_network(self, orig_net_path: Path, target_genes: set, filtered_net_path: Path) -> bool:
        if not orig_net_path.exists():
            return False
        try:
            df = pd.read_csv(orig_net_path, sep="\t")
            if df.empty or "Source" not in df.columns or "Target" not in df.columns:
                filtered_net_path.parent.mkdir(parents=True, exist_ok=True)
                df.to_csv(filtered_net_path, sep="\t", index=False)
                return True

            df["Source"] = df["Source"].astype(str)
            df["Target"] = df["Target"].astype(str)

            targets_set = set(target_genes)

            source_in_targets = df["Source"].isin(targets_set)
            target_in_targets = df["Target"].isin(targets_set)

            neighbors = set(df.loc[source_in_targets, "Target"]).union(
                set(df.loc[target_in_targets, "Source"])
            ) - targets_set

            remaining_genes = targets_set.union(neighbors)

            filtered_df = df[df["Source"].isin(remaining_genes) & df["Target"].isin(remaining_genes)].copy()

            filtered_net_path.parent.mkdir(parents=True, exist_ok=True)
            filtered_df.to_csv(filtered_net_path, sep="\t", index=False)
            return True
        except Exception as e:
            print(f"[Error] Failed to filter network {orig_net_path} with target genes: {e}")
            return False

    def _get_unique_target_names(self, target_files: List[Path]) -> List[tuple[Path, str]]:
        unique_names = []
        seen = {}
        for tf in target_files:
            tf_path = Path(tf).expanduser().resolve()
            base = tf_path.stem
            if base in seen:
                seen[base] += 1
                name = f"{base}_{seen[base]}"
            else:
                seen[base] = 0
                name = base
            unique_names.append((tf_path, name))
        return unique_names

    def _create_1x2_saturation_plot(self, source_data_list, samp_df, show_samples, title_prefix, out_file):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        ax1_twin = ax1.twinx() if show_samples else None
        ax2_twin = ax2.twinx() if show_samples else None

        roc_plotted, prc_plotted = False, False
        lines_ax1, labels_ax1, lines_ax2, labels_ax2 = [], [], [], []

        for data in source_data_list:
            src, color, df = data['src'], data['color'], data['df']

            if self.do_auroc and "auc_mean" in df.columns:
                roc_plotted = True
                l1 = ax1.plot(df["n_bps"], df["auc_mean"], color=color, linewidth=2, marker='o', label=f"{src} AUROC")
                ax1.errorbar(df["n_bps"], df["auc_mean"], yerr=df["auc_std"], fmt='none', ecolor=color, capsize=3,
                             alpha=0.6)
                lines_ax1.extend(l1)
                labels_ax1.append(f"{src} AUROC")

            if self.do_aupr and "aupr_mean" in df.columns:
                prc_plotted = True
                l2 = ax2.plot(df["n_bps"], df["aupr_mean"], color=color, linewidth=2, marker='^', label=f"{src} AUPR")
                ax2.errorbar(df["n_bps"], df["aupr_mean"], yerr=df["aupr_std"], fmt='none', ecolor=color, capsize=3,
                             alpha=0.6)
                lines_ax2.extend(l2)
                labels_ax2.append(f"{src} AUPR")

        if show_samples and not samp_df.empty:
            color_samp = 'gray'
            ls1 = ax1_twin.plot(samp_df["n_bps"], samp_df["n_samples_mean"], color=color_samp, linewidth=2,
                                linestyle='--', marker='s', label="Samples")
            ax1_twin.errorbar(samp_df["n_bps"], samp_df["n_samples_mean"], yerr=samp_df["n_samples_std"], fmt='none',
                              ecolor=color_samp, capsize=3, alpha=0.5)
            lines_ax1.extend(ls1)
            labels_ax1.append("Samples")
            ax1_twin.set_ylabel('Number of Samples', color=color_samp)
            ax1_twin.tick_params(axis='y', labelcolor=color_samp)

            ls2 = ax2_twin.plot(samp_df["n_bps"], samp_df["n_samples_mean"], color=color_samp, linewidth=2,
                                linestyle='--', marker='s', label="Samples")
            ax2_twin.errorbar(samp_df["n_bps"], samp_df["n_samples_mean"], yerr=samp_df["n_samples_std"], fmt='none',
                              ecolor=color_samp, capsize=3, alpha=0.5)
            lines_ax2.extend(ls2)
            labels_ax2.append("Samples")
            ax2_twin.set_ylabel('Number of Samples', color=color_samp)
            ax2_twin.tick_params(axis='y', labelcolor=color_samp)

        if roc_plotted:
            ax1.set_title(f"AUROC Saturation")
            ax1.set_xlabel("Number of BioProjects")
            ax1.set_ylabel("Macro-averaged AUROC")
            ax1.legend(lines_ax1, labels_ax1, loc="best", fontsize='small')
            ax1.grid(True, linestyle='--', alpha=0.6)
        else:
            ax1.text(0.5, 0.5, "AUROC Data Unavailable", ha='center', va='center')
            if show_samples: ax1_twin.set_yticks([])

        if prc_plotted:
            ax2.set_title(f"AUPR Saturation")
            ax2.set_xlabel("Number of BioProjects")
            ax2.set_ylabel("Macro-averaged AUPR")
            ax2.legend(lines_ax2, labels_ax2, loc="best", fontsize='small')
            ax2.grid(True, linestyle='--', alpha=0.6)
        else:
            ax2.text(0.5, 0.5, "AUPR Data Unavailable", ha='center', va='center')
            if show_samples: ax2_twin.set_yticks([])

        fig.tight_layout()
        fig.savefig(out_file)
        plt.close(fig)

    def _generate_plots(self, results: List[dict], outdir: Optional[Path] = None):
        if outdir is None:
            outdir = self.base_outdir

        if not results:
            print("[Saturation] No results to plot.")
            return

        df = pd.DataFrame(results)
        df.to_csv(outdir / "saturation_results_raw.tsv", sep="\t", index=False)

        samp_agg = df[["n_bps", "iter", "n_samples"]].drop_duplicates().groupby("n_bps").agg(
            n_samples_mean=("n_samples", "mean"),
            n_samples_std=("n_samples", "std")
        ).reset_index().fillna(0)

        if "source" not in df.columns: df["source"] = "Legacy"

        agg_dict = {}
        if self.do_auroc and "auc" in df.columns: agg_dict["auc"] = ["mean", "std"]
        if self.do_aupr and "aupr" in df.columns: agg_dict["aupr"] = ["mean", "std"]

        if not agg_dict:
            print("[Saturation] No metric data to plot curves for.")
            return

        metric_agg = df.groupby(["n_bps", "source"]).agg(agg_dict).reset_index().fillna(0)
        metric_agg.columns = ['_'.join(col).strip('_') for col in metric_agg.columns.values]

        summary_df = pd.merge(metric_agg, samp_agg, on="n_bps", how="left")
        summary_df.to_csv(outdir / "saturation_results_summary.tsv", sep="\t", index=False)

        colors = {"MapMan": "firebrick", "GO": "darkorange", "Legacy": "purple"}
        sources = metric_agg["source"].unique()
        source_data_list = [
            {'src': src, 'color': colors.get(src, "purple"), 'df': metric_agg[metric_agg["source"] == src]} for src in
            sources]

        try:
            for data in source_data_list:
                src = data['src']
                self._create_1x2_saturation_plot([data], samp_agg, False, src,
                                                 outdir / f"{src}_no_samples.pdf")
                self._create_1x2_saturation_plot([data], samp_agg, True, src,
                                                 outdir / f"{src}_with_samples.pdf")

            if len(source_data_list) > 1:
                self._create_1x2_saturation_plot(source_data_list, samp_agg, False, "",
                                                 outdir / "combined_no_samples.pdf")
                self._create_1x2_saturation_plot(source_data_list, samp_agg, True, "",
                                                 outdir / "combined_with_samples.pdf")

            print(f"\n[Saturation] Plots updated in {outdir}")
        except Exception as e:
            print(f"[Plot Error] {e}")

    def regenerate_plots_from_raw(self):
        raw_file = self.base_outdir / "saturation_results_raw.tsv"
        if raw_file.exists():
            try:
                print(f"[Saturation] Plot-only mode active. Reading {raw_file}...")
                df = pd.read_csv(raw_file, sep="\t")
                self._generate_plots(df.to_dict('records'), outdir=self.base_outdir)
            except Exception as e:
                print(f"[Error] Failed to read {raw_file}: {e}")
        else:
            print(f"[Error] Cannot plot main saturation. Missing {raw_file}.")

        if self.target_genes_files:
            target_specs = self._get_unique_target_names(self.target_genes_files)
            for tf, target_name in target_specs:
                target_outdir = self.base_outdir / target_name
                t_raw = target_outdir / "saturation_results_raw.tsv"
                if t_raw.exists():
                    try:
                        print(f"[Saturation] [Target Genes: {target_name}] Reading {t_raw}...")
                        tdf = pd.read_csv(t_raw, sep="\t")
                        self._generate_plots(tdf.to_dict('records'), outdir=target_outdir)
                    except Exception as e:
                        print(f"[Error] Failed to read {t_raw}: {e}")
                else:
                    print(f"[Warn] Cannot plot target saturation for {target_name}. Missing {t_raw}.")

    def run(self):
        setup_interrupt_handlers()
        state_file = self.base_outdir / "run_state.json"
        target_specs = self._get_unique_target_names(self.target_genes_files) if self.target_genes_files else []

        current_state = {
            "seed": self.seed,
            "iterations": self.iterations,
            "num_steps": self.num_steps,
            "seidr_preset": getattr(self.config, "seidr_preset", "FAST"),
            "mapman_file": str(self.mapman_file) if self.mapman_file else None,
            "go_file": str(self.go_file) if self.go_file else None,
            "do_auroc": self.do_auroc,
            "do_aupr": self.do_aupr,
            "target_genes_files": [str(tf) for tf, _ in target_specs]
        }

        force_wipe_level = 0

        if self.force:
            force_wipe_level = 3
        elif state_file.exists():
            try:
                with open(state_file, "r") as f:
                    old_state = json.load(f)

                if old_state.get("seed") != current_state["seed"] or \
                        old_state.get("iterations") != current_state["iterations"] or \
                        old_state.get("num_steps") != current_state["num_steps"]:
                    print("[Saturation] Structural parameters changed. Forcing Level 3 wipe (Complete Restart).")
                    force_wipe_level = 3
                elif old_state.get("seidr_preset") != current_state["seidr_preset"]:
                    print("[Saturation] Seidr preset changed. Forcing Level 2 wipe (Restarting from Inference).")
                    force_wipe_level = 2
                elif old_state.get("mapman_file") != current_state["mapman_file"] or \
                        old_state.get("go_file") != current_state["go_file"] or \
                        old_state.get("do_auroc") != current_state["do_auroc"] or \
                        old_state.get("do_aupr") != current_state["do_aupr"]:
                    print("[Saturation] Evaluation metrics changed. Forcing Level 1 wipe (Restarting EGAD Analysis).")
                    force_wipe_level = 1
                else:
                    print("[Saturation] Parameters identical. Resuming execution using phase markers.")

            except Exception:
                print("[Warn] Corrupted state file. Forcing Level 3 complete restart.")
                force_wipe_level = 3
        else:
            force_wipe_level = 3

        if force_wipe_level == 3 and self.base_outdir.exists():
            print(f"[Saturation] Purging entire directory: {self.base_outdir}...")
            for item in self.base_outdir.iterdir():
                try:
                    if item.is_dir():
                        shutil.rmtree(item)
                    else:
                        item.unlink()
                except Exception:
                    pass
        elif force_wipe_level == 2:
            print("[Saturation] Purging previous Network and Evaluation outputs...")
            for iter_dir in self.base_outdir.glob("step*/iter*"):
                for f in [".seidr.done", ".egad.done", "network_saturation_edges.tsv", "egad_results.tsv",
                          "seidr_batch.log", "egad.log"]:
                    (iter_dir / f).unlink(missing_ok=True)
            (self.base_outdir / ".saturation.done").unlink(missing_ok=True)
        elif force_wipe_level == 1:
            print("[Saturation] Purging previous Evaluation outputs...")
            for iter_dir in self.base_outdir.glob("step*/iter*"):
                for f in [".egad.done", "egad_results.tsv", "egad.log"]:
                    (iter_dir / f).unlink(missing_ok=True)
            (self.base_outdir / ".saturation.done").unlink(missing_ok=True)

        self.base_outdir.mkdir(parents=True, exist_ok=True)
        with open(state_file, "w") as f:
            json.dump(current_state, f, indent=4)

        bp_meta = []
        for bp in self.dataset.bioprojects:
            s_ids = bp.get_sample_ids()
            if s_ids: bp_meta.append({'obj': bp, 'size': len(s_ids), 'ids': s_ids})

        total_bps = len(bp_meta)
        if total_bps == 0: return
        avg_samples_per_bp = sum(x['size'] for x in bp_meta) / total_bps
        bp_steps = self._calculate_steps(total_bps)
        step_name_map = {n: f"step{i + 1}" for i, n in enumerate(bp_steps)}

        print(f"\n[Saturation] BioProject Steps (Total {total_bps} BPs):")
        for n, sname in step_name_map.items(): print(f"  {sname}: {n} BPs")

        work_definitions = []
        for n_bps in bp_steps:
            iters = 1 if n_bps == total_bps else self.iterations
            for i in range(1, iters + 1):
                work_definitions.append((n_bps, i, step_name_map[n_bps]))

        seidr_queue, egad_queue = [], []
        global_genes, global_expr = self.deseq2_dir / "genes.txt", self.deseq2_dir / "expression.tsv"
        matrix = self._get_master_df()

        step_pools = {n_bps: list(bp_meta) for n_bps in bp_steps}

        print(f"\n[Phase 1] Checking/Generating Batches...")
        with tqdm(total=len(work_definitions), desc="Batch Gen", unit="batch") as pbar:
            for n_bps, iter_num, s_name in work_definitions:
                iter_dir = self.base_outdir / s_name / f"iter{iter_num}"
                iter_dir.mkdir(parents=True, exist_ok=True)
                g_path = global_genes if n_bps == total_bps else iter_dir / "genes.txt"
                e_path = global_expr if n_bps == total_bps else iter_dir / "expression.tsv"

                seidr_queue.append((g_path, e_path, iter_dir, s_name, iter_num))
                egad_queue.append({
                    "n_bps": n_bps, "iter": iter_num, "dir": iter_dir,
                    "net": iter_dir / "network_saturation_edges.tsv",
                    "out": iter_dir / "egad_results.tsv", "expr": e_path
                })

                if not (iter_dir / ".batch.done").exists():
                    if n_bps < total_bps:

                        if len(step_pools[n_bps]) < n_bps:
                            step_pools[n_bps] = list(bp_meta)

                        target_samples = n_bps * avg_samples_per_bp
                        current_selection = random.sample(step_pools[n_bps], n_bps)

                        step_pools[n_bps] = [x for x in step_pools[n_bps] if x not in current_selection]

                        current_n = sum(x['size'] for x in current_selection)
                        current_diff = abs(current_n - target_samples)

                        remaining_pool = list(step_pools[n_bps])
                        improved, max_swaps, swaps = True, 50, 0

                        while improved and swaps < max_swaps:
                            improved = False
                            for i in range(len(current_selection)):
                                for j in range(len(remaining_pool)):
                                    new_n = current_n - current_selection[i]['size'] + remaining_pool[j]['size']
                                    new_diff = abs(new_n - target_samples)
                                    if new_diff < current_diff:
                                        current_n, current_diff = new_n, new_diff
                                        current_selection[i], remaining_pool[j] = remaining_pool[j], current_selection[
                                            i]
                                        improved, swaps = True, swaps + 1
                                        break
                                if improved: break

                        step_pools[n_bps] = list(remaining_pool)

                        selected_samples = []
                        for item in current_selection: selected_samples.extend(item['ids'])
                        final_s = [s for s in selected_samples if not matrix.empty and s in matrix.index]
                        if final_s:
                            sub = matrix.loc[final_s].copy().fillna(0)
                            sub = sub.loc[:, sub.var(axis=0) > 0.1]
                            if not sub.empty: self._write_seidr_files(sub, iter_dir)

                    (iter_dir / ".batch.done").touch()

                pbar.update(1)

        print(f"\n[Phase 2] Seidr Inference...")
        seidr_tasks = [i for i in seidr_queue if not (i[2] / ".seidr.done").exists()]
        if seidr_tasks:
            threads = max(1, self.total_thread_budget // self.workers)
            preset = getattr(self.config, "seidr_preset", "FAST")
            with tqdm(total=len(seidr_queue), initial=len(seidr_queue) - len(seidr_tasks), desc="Seidr Inf") as p_inf:
                # Group tasks by step (n_bps) so lower steps complete sequentially
                step_groups = {}
                for task in seidr_tasks:
                    # task format: (g_path, e_path, iter_dir, s_name, iter_num)
                    step_key = task[3]  # s_name (e.g. step1, step2)
                    step_groups.setdefault(step_key, []).append(task)

                for step_name in sorted(step_groups.keys(), key=lambda x: int(x.replace('step', '')) if x.replace('step', '').isdigit() else x):
                    group = step_groups[step_name]
                    with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                        f_map = {executor.submit(run_seidr_batch, self.config, i[0], i[1], i[2], preset, threads): i for i in group}
                        for f in concurrent.futures.as_completed(f_map):
                            it = f_map[f]
                            try:
                                f.result()
                            except Exception as e:
                                print(f"[Warn] Seidr failed for {it[3]} iter {it[4]}: {e}")
                            p_inf.update(1)

        final_results = []
        if self.mapman_file or self.go_file:
            print(f"\n[Phase 3] EGAD Analysis...")
            r_script = get_egad_script_path()
            egad_tasks = []

            for item in egad_queue:
                if (item["dir"] / ".egad.done").exists() and item["out"].exists():
                    try:
                        df_check = pd.read_csv(item["out"], sep="\t")
                        if not df_check.empty:
                            if "Annotation_Source" in df_check.columns:
                                for src, group in df_check.groupby("Annotation_Source"):
                                    valid_group = group[group["Term"] != "None"]
                                    if valid_group.empty: continue
                                    res_dict = {
                                        'n_bps': item['n_bps'], 'iter': item['iter'],
                                        'n_samples': self._count_samples_from_file(item["expr"]),
                                        'source': src
                                    }
                                    if self.do_auroc and "AUC" in valid_group.columns: res_dict['auc'] = valid_group["AUC"].mean()
                                    if self.do_aupr and "AUPR" in valid_group.columns: res_dict['aupr'] = valid_group["AUPR"].mean()
                                    final_results.append(res_dict)
                            else:
                                res_dict = {
                                    'n_bps': item['n_bps'], 'iter': item['iter'],
                                    'n_samples': self._count_samples_from_file(item["expr"]),
                                    'source': "Legacy"
                                }
                                if self.do_auroc and "AUC" in df_check.columns: res_dict['auc'] = df_check["AUC"].mean()
                                if self.do_aupr and "AUPR" in df_check.columns: res_dict['aupr'] = df_check["AUPR"].mean()
                                final_results.append(res_dict)
                    except Exception as e:
                        print(f"[Warn] Failed reading {item['out']}: {e}")
                else:
                    if (item["dir"] / ".seidr.done").exists():
                        egad_tasks.append(item)

            if final_results:
                self._generate_plots(final_results, outdir=self.base_outdir)

            if egad_tasks:
                with tqdm(total=len(egad_queue), initial=len(egad_queue) - len(egad_tasks), desc="EGAD Eval") as p_egad:
                    # Group tasks by n_bps to process steps sequentially and save incrementally
                    step_egad_groups = {}
                    for item in egad_tasks:
                        step_egad_groups.setdefault(item["n_bps"], []).append(item)

                    for n_bps in sorted(step_egad_groups.keys()):
                        group = step_egad_groups[n_bps]
                        with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                            f_map = {
                                executor.submit(
                                    run_egad_task, network_file=i["net"], mapman_file=self.mapman_file,
                                    go_file=self.go_file, out_file=i["out"], script_path=r_script,
                                    log_path=i["dir"] / "egad.log", curves_prefix=None,
                                    do_auroc=self.do_auroc, do_aupr=self.do_aupr
                                ): i for i in group
                            }

                            for f in concurrent.futures.as_completed(f_map):
                                it = f_map[f]
                                try:
                                    results_dict = f.result()
                                    if results_dict:
                                        (it["dir"] / ".egad.done").touch()
                                        for src, mets in results_dict.items():
                                            res_dict = {
                                                'n_bps': it['n_bps'], 'iter': it['iter'],
                                                'n_samples': self._count_samples_from_file(it["expr"]), 'source': src
                                            }
                                            if self.do_auroc and mets.get('macro_auc') is not None: res_dict['auc'] = mets['macro_auc']
                                            if self.do_aupr and mets.get('macro_aupr') is not None: res_dict['aupr'] = mets['macro_aupr']
                                            final_results.append(res_dict)
                                except Exception as e:
                                    print(f"[Warn] EGAD task failed for n_bps {it['n_bps']} iter {it['iter']}: {e}")
                                p_egad.update(1)

                        # Save progress incrementally after each step
                        if final_results:
                            self._generate_plots(final_results, outdir=self.base_outdir)

            (self.base_outdir / ".saturation.done").touch()

        # Phase 4: Target Genes Sub-Workflows
        if target_specs:
            print(f"\n[Phase 4] Target Genes Workflows ({len(target_specs)} target files)...")
            r_script = get_egad_script_path() if (self.mapman_file or self.go_file) else None

            for tf, target_name in target_specs:
                target_genes = self._read_target_genes(tf)
                if not target_genes:
                    print(f"[Warn] Target file {tf} is empty or unreadable. Skipping.")
                    continue

                target_outdir = self.base_outdir / target_name
                target_outdir.mkdir(parents=True, exist_ok=True)

                print(f"\n[Target Genes: {target_name}] Filtering step-batch networks for {len(target_genes)} target genes...")

                for n_bps, iter_num, s_name in work_definitions:
                    iter_dir = self.base_outdir / s_name / f"iter{iter_num}"
                    orig_net = iter_dir / "network_saturation_edges.tsv"
                    filtered_net = iter_dir / f"{target_name}_edges.tsv"

                    if orig_net.exists():
                        if self.force or not filtered_net.exists() or force_wipe_level > 0:
                            self._filter_network(orig_net, target_genes, filtered_net)

                if self.mapman_file or self.go_file:
                    target_results = []
                    target_egad_queue = []

                    for n_bps, iter_num, s_name in work_definitions:
                        iter_dir = self.base_outdir / s_name / f"iter{iter_num}"
                        filtered_net = iter_dir / f"{target_name}_edges.tsv"
                        target_iter_dir = target_outdir / s_name / f"iter{iter_num}"
                        target_iter_dir.mkdir(parents=True, exist_ok=True)
                        out_tsv = target_iter_dir / "egad_results.tsv"
                        e_path = global_expr if n_bps == total_bps else iter_dir / "expression.tsv"

                        item = {
                            "n_bps": n_bps, "iter": iter_num, "dir": target_iter_dir,
                            "net": filtered_net, "out": out_tsv, "expr": e_path
                        }

                        if (target_iter_dir / ".egad.done").exists() and out_tsv.exists() and not self.force and force_wipe_level == 0:
                            try:
                                df_check = pd.read_csv(out_tsv, sep="\t")
                                if not df_check.empty:
                                    if "Annotation_Source" in df_check.columns:
                                        for src, group in df_check.groupby("Annotation_Source"):
                                            valid_group = group[group["Term"] != "None"]
                                            if valid_group.empty: continue
                                            res_dict = {
                                                'n_bps': item['n_bps'], 'iter': item['iter'],
                                                'n_samples': self._count_samples_from_file(item["expr"]),
                                                'source': src
                                            }
                                            if self.do_auroc and "AUC" in valid_group.columns:
                                                res_dict['auc'] = valid_group["AUC"].mean()
                                            if self.do_aupr and "AUPR" in valid_group.columns:
                                                res_dict['aupr'] = valid_group["AUPR"].mean()
                                            target_results.append(res_dict)
                                    else:
                                        res_dict = {
                                            'n_bps': item['n_bps'], 'iter': item['iter'],
                                            'n_samples': self._count_samples_from_file(item["expr"]),
                                            'source': "Legacy"
                                        }
                                        if self.do_auroc and "AUC" in df_check.columns: res_dict['auc'] = df_check["AUC"].mean()
                                        if self.do_aupr and "AUPR" in df_check.columns: res_dict['aupr'] = df_check["AUPR"].mean()
                                        target_results.append(res_dict)
                            except Exception as e:
                                print(f"[Warn] Failed reading {out_tsv}: {e}")
                                target_egad_queue.append(item)
                        else:
                            if filtered_net.exists():
                                target_egad_queue.append(item)

                    if target_results:
                        self._generate_plots(target_results, outdir=target_outdir)

                    if target_egad_queue:
                        with tqdm(total=len(work_definitions), initial=len(work_definitions) - len(target_egad_queue),
                                  desc=f"EGAD [{target_name}]") as p_egad:
                            step_target_groups = {}
                            for item in target_egad_queue:
                                step_target_groups.setdefault(item["n_bps"], []).append(item)

                            for n_bps in sorted(step_target_groups.keys()):
                                group = step_target_groups[n_bps]
                                with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                                    f_map = {
                                        executor.submit(
                                            run_egad_task, network_file=i["net"], mapman_file=self.mapman_file,
                                            go_file=self.go_file, out_file=i["out"], script_path=r_script,
                                            log_path=i["dir"] / "egad.log", curves_prefix=None,
                                            do_auroc=self.do_auroc, do_aupr=self.do_aupr
                                        ): i for i in group
                                    }

                                    for f in concurrent.futures.as_completed(f_map):
                                        it = f_map[f]
                                        try:
                                            results_dict = f.result()
                                            if results_dict:
                                                (it["dir"] / ".egad.done").touch()
                                                for src, mets in results_dict.items():
                                                    res_dict = {
                                                        'n_bps': it['n_bps'], 'iter': it['iter'],
                                                        'n_samples': self._count_samples_from_file(it["expr"]), 'source': src
                                                    }
                                                    if self.do_auroc and mets.get('macro_auc') is not None:
                                                        res_dict['auc'] = mets['macro_auc']
                                                    if self.do_aupr and mets.get('macro_aupr') is not None:
                                                        res_dict['aupr'] = mets['macro_aupr']
                                                    target_results.append(res_dict)
                                        except Exception as e:
                                            print(f"[Warn] Target EGAD task failed for [{target_name}] n_bps {it['n_bps']} iter {it['iter']}: {e}")
                                        p_egad.update(1)

                                # Save progress incrementally after each step for target
                                if target_results:
                                    self._generate_plots(target_results, outdir=target_outdir)

                    (target_outdir / ".saturation.done").touch()