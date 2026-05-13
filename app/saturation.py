import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, List, Tuple, TYPE_CHECKING, Dict, Any
import concurrent.futures
from tqdm import tqdm
import random
import os
import shutil

# --- MATPLOTLIB SILENCER & SETUP ---
os.environ["MPLCONFIGDIR"] = "/tmp"
import matplotlib.pyplot as plt

from .seidr import run_seidr_batch
from .egad import run_egad_task, get_egad_script_path

if TYPE_CHECKING:
    from .entities import Dataset, Config


class BatchOrchestrator:
    def __init__(self, dataset: "Dataset", config: "Config",
                 seed: Optional[int] = None,
                 workers: int = 4,
                 max_threads: Optional[int] = None,
                 mapman_file: Optional[Path] = None,
                 go_file: Optional[Path] = None,
                 force: bool = False,
                 num_steps: int = 10,
                 metrics: str = "both"):
        self.dataset = dataset
        self.config = config
        self.seed = seed
        self.workers = workers
        self.mapman_file = mapman_file
        self.go_file = go_file
        self.force = force
        self.num_steps = num_steps

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
        """Matrix is Sample x Gene. Seidr needs no headers/index."""
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
        if total_bps <= self.num_steps:
            return list(range(1, total_bps + 1))
        base = total_bps // self.num_steps
        remainder = total_bps % self.num_steps
        increments = [base] * self.num_steps
        for i in range(1, remainder + 1):
            increments[-i] += 1
        steps, current = [], 0
        for inc in increments:
            current += inc
            steps.append(current)
        return sorted(list(set(steps)))

    def _generate_plots(self, results: List[dict]):
        if not results:
            print("[Saturation] No results to plot.")
            return

        df = pd.DataFrame(results)
        df.to_csv(self.base_outdir / "saturation_results_raw.tsv", sep="\t", index=False)

        # 1. Aggregate samples independently (they don't care about MapMan vs GO)
        samp_agg = df[["n_bps", "iter", "n_samples"]].drop_duplicates().groupby("n_bps").agg(
            n_samples_mean=("n_samples", "mean"),
            n_samples_std=("n_samples", "std")
        ).reset_index().fillna(0)

        color_samp = '#d35400'
        colors = {
            "MapMan": {"auc": "royalblue", "aupr": "firebrick"},
            "GO": {"auc": "lightseagreen", "aupr": "orange"},
            "Legacy": {"auc": "#2c3e50", "aupr": "#8e44ad"}
        }

        # Plot Single Sample Accumulation
        try:
            plt.figure(figsize=(8, 6))
            plt.plot(samp_agg["n_bps"], samp_agg["n_samples_mean"], color=color_samp, linewidth=2, marker='s')
            plt.errorbar(samp_agg["n_bps"], samp_agg["n_samples_mean"], yerr=samp_agg["n_samples_std"], fmt='none',
                         ecolor='black', capsize=3)
            plt.title("Sample Accumulation vs BioProject Count")
            plt.xlabel("Number of BioProjects")
            plt.ylabel("Total Samples")
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.tight_layout()
            plt.savefig(self.base_outdir / "saturation_samples.pdf")
            plt.close()
        except Exception as e:
            print(f"[Plot Error - Samples] {e}")

        # 2. Aggregate Metrics grouping by BioProjects AND Annotation Source
        if "source" not in df.columns:
            df["source"] = "Legacy"

        agg_dict = {}
        if self.do_auroc and "auc" in df.columns: agg_dict["auc"] = ["mean", "std"]
        if self.do_aupr and "aupr" in df.columns: agg_dict["aupr"] = ["mean", "std"]

        if not agg_dict:
            print("[Saturation] No metric data to plot curves for.")
            return

        metric_agg = df.groupby(["n_bps", "source"]).agg(agg_dict).reset_index().fillna(0)
        metric_agg.columns = ['_'.join(col).strip('_') for col in metric_agg.columns.values]

        # Merge sample stats back into the flat summary
        summary_df = pd.merge(metric_agg, samp_agg, on="n_bps", how="left")
        summary_df.to_csv(self.base_outdir / "saturation_results_summary.tsv", sep="\t", index=False)

        sources = metric_agg["source"].unique()

        # --- AUROC PLOTS ---
        if self.do_auroc and "auc_mean" in metric_agg.columns:
            try:
                plt.figure(figsize=(8, 6))
                for src in sources:
                    src_data = metric_agg[metric_agg["source"] == src]
                    c = colors.get(src, colors["Legacy"])["auc"]
                    plt.plot(src_data["n_bps"], src_data["auc_mean"], color=c, linewidth=2, marker='o',
                             label=f"{src} AUROC")
                    plt.errorbar(src_data["n_bps"], src_data["auc_mean"], yerr=src_data["auc_std"], fmt='none',
                                 ecolor=c, capsize=3, alpha=0.6)

                plt.title("Performance vs BioProject Count (AUROC)")
                plt.xlabel("Number of BioProjects")
                plt.ylabel("Mean AUROC")
                plt.legend(loc="best")
                plt.grid(True, linestyle='--', alpha=0.6)
                plt.tight_layout()
                plt.savefig(self.base_outdir / "saturation_auroc.pdf")
                plt.close()

                # Combined Plot
                fig, ax1 = plt.subplots(figsize=(10, 6))
                ax1.set_xlabel('Number of BioProjects')
                ax1.set_ylabel('Mean AUROC')

                lines = []
                for src in sources:
                    src_data = metric_agg[metric_agg["source"] == src]
                    c = colors.get(src, colors["Legacy"])["auc"]
                    l = ax1.plot(src_data["n_bps"], src_data["auc_mean"], color=c, linewidth=2, marker='o',
                                 label=f"{src} AUROC")
                    ax1.errorbar(src_data["n_bps"], src_data["auc_mean"], yerr=src_data["auc_std"], fmt='none',
                                 ecolor=c, capsize=3, alpha=0.5)
                    lines.extend(l)

                ax2 = ax1.twinx()
                ax2.set_ylabel('Mean Samples', color=color_samp)
                l2 = ax2.plot(samp_agg["n_bps"], samp_agg["n_samples_mean"], color=color_samp, linewidth=2,
                              linestyle='--', marker='s', label="Samples")
                ax2.errorbar(samp_agg["n_bps"], samp_agg["n_samples_mean"], yerr=samp_agg["n_samples_std"], fmt='none',
                             ecolor='black', capsize=3, alpha=0.5)
                lines.extend(l2)

                plt.title("Saturation Analysis Summary (AUROC)")
                labels = [l.get_label() for l in lines]
                ax1.legend(lines, labels, loc='best')

                fig.tight_layout()
                plt.savefig(self.base_outdir / "saturation_combined_auroc.pdf")
                plt.close()
            except Exception as e:
                print(f"[Plot Error - AUROC] {e}")

        # --- AUPR PLOTS ---
        if self.do_aupr and "aupr_mean" in metric_agg.columns:
            try:
                plt.figure(figsize=(8, 6))
                for src in sources:
                    src_data = metric_agg[metric_agg["source"] == src]
                    c = colors.get(src, colors["Legacy"])["aupr"]
                    plt.plot(src_data["n_bps"], src_data["aupr_mean"], color=c, linewidth=2, marker='^',
                             label=f"{src} AUPR")
                    plt.errorbar(src_data["n_bps"], src_data["aupr_mean"], yerr=src_data["aupr_std"], fmt='none',
                                 ecolor=c, capsize=3, alpha=0.6)

                plt.title("Performance vs BioProject Count (AUPR)")
                plt.xlabel("Number of BioProjects")
                plt.ylabel("Mean AUPR")
                plt.legend(loc="best")
                plt.grid(True, linestyle='--', alpha=0.6)
                plt.tight_layout()
                plt.savefig(self.base_outdir / "saturation_aupr.pdf")
                plt.close()

                # Combined Plot
                fig, ax1 = plt.subplots(figsize=(10, 6))
                ax1.set_xlabel('Number of BioProjects')
                ax1.set_ylabel('Mean AUPR')

                lines = []
                for src in sources:
                    src_data = metric_agg[metric_agg["source"] == src]
                    c = colors.get(src, colors["Legacy"])["aupr"]
                    l = ax1.plot(src_data["n_bps"], src_data["aupr_mean"], color=c, linewidth=2, marker='^',
                                 label=f"{src} AUPR")
                    ax1.errorbar(src_data["n_bps"], src_data["aupr_mean"], yerr=src_data["aupr_std"], fmt='none',
                                 ecolor=c, capsize=3, alpha=0.5)
                    lines.extend(l)

                ax2 = ax1.twinx()
                ax2.set_ylabel('Mean Samples', color=color_samp)
                l2 = ax2.plot(samp_agg["n_bps"], samp_agg["n_samples_mean"], color=color_samp, linewidth=2,
                              linestyle='--', marker='s', label="Samples")
                ax2.errorbar(samp_agg["n_bps"], samp_agg["n_samples_mean"], yerr=samp_agg["n_samples_std"], fmt='none',
                             ecolor='black', capsize=3, alpha=0.5)
                lines.extend(l2)

                plt.title("Saturation Analysis Summary (AUPR)")
                labels = [l.get_label() for l in lines]
                ax1.legend(lines, labels, loc='best')

                fig.tight_layout()
                plt.savefig(self.base_outdir / "saturation_combined_aupr.pdf")
                plt.close()
            except Exception as e:
                print(f"[Plot Error - AUPR] {e}")

        print(f"\n[Saturation] Plots updated in {self.base_outdir}")

    def regenerate_plots_from_raw(self):
        raw_file = self.base_outdir / "saturation_results_raw.tsv"
        if not raw_file.exists():
            print(f"[Error] Cannot plot. Missing {raw_file}. Did you even run it completely once?")
            return
        try:
            print(f"[Saturation] Plot-only mode active. Reading {raw_file}...")
            df = pd.read_csv(raw_file, sep="\t")

            required_cols = ["n_bps", "n_samples"]
            if self.do_auroc: required_cols.append("auc")
            if self.do_aupr: required_cols.append("aupr")

            if not all(col in df.columns for col in required_cols):
                print(
                    f"[Error] The file {raw_file} is malformed. Missing required columns based on your metrics request.")
                return
            self._generate_plots(df.to_dict('records'))
        except Exception as e:
            print(f"[Error] Failed to read {raw_file}: {e}")

    def run(self):
        if self.force and self.base_outdir.exists():
            print(f"[Saturation] Force=True: Purging {self.base_outdir}...")
            for item in self.base_outdir.iterdir():
                try:
                    if item.is_dir():
                        shutil.rmtree(item)
                    else:
                        item.unlink()
                except Exception as e:
                    print(f"[Warn] Failed to delete {item}: {e}")
        else:
            print(f"[Saturation] Smart Resume: Reusing existing files.")

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
        for n, sname in step_name_map.items():
            print(f"  {sname}: {n} BPs")

        work_definitions = []
        for n_bps in bp_steps:
            iters = 1 if n_bps == total_bps else self.iterations
            for i in range(1, iters + 1):
                work_definitions.append((n_bps, i, step_name_map[n_bps]))

        seidr_queue, egad_queue = [], []
        global_genes, global_expr = self.deseq2_dir / "genes.txt", self.deseq2_dir / "expression.tsv"
        matrix = self._get_master_df()

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

                if n_bps < total_bps:
                    if not self.force and g_path.exists() and e_path.exists():
                        pass
                    else:
                        target_samples = n_bps * avg_samples_per_bp
                        current_selection = random.sample(bp_meta, n_bps)
                        current_n = sum(x['size'] for x in current_selection)
                        current_diff = abs(current_n - target_samples)
                        remaining_pool = [x for x in bp_meta if x not in current_selection]
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

                        selected_samples = []
                        for item in current_selection: selected_samples.extend(item['ids'])
                        final_s = [s for s in selected_samples if not matrix.empty and s in matrix.index]
                        if final_s:
                            sub = matrix.loc[final_s].copy().fillna(0)
                            sub = sub.loc[:, sub.var(axis=0) > 0.1]
                            if not sub.empty: self._write_seidr_files(sub, iter_dir)
                pbar.update(1)

        print(f"\n[Phase 2] Seidr Inference...")
        tasks = [i for i in seidr_queue if self.force or not (i[2] / ".saturation.done").exists()]
        if tasks:
            threads = max(1, self.total_thread_budget // self.workers)
            preset = getattr(self.config, "seidr_preset", "FAST")
            with tqdm(total=len(seidr_queue), initial=len(seidr_queue) - len(tasks), desc="Seidr Inf") as p_inf:
                with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                    f_map = {executor.submit(run_seidr_batch, self.config, i[0], i[1], i[2], preset, threads): i for i
                             in tasks}
                    for f in concurrent.futures.as_completed(f_map): p_inf.update(1)

        final_results = []
        if self.mapman_file or self.go_file:
            print(f"\n[Phase 3] EGAD Analysis...")
            r_script = get_egad_script_path()

            egad_tasks = []
            for item in egad_queue:
                # 1. FIRST, check if we already have the valid output
                if not self.force and item["out"].exists():
                    try:
                        df_check = pd.read_csv(item["out"], sep="\t")
                        if not df_check.empty:
                            if "Annotation_Source" in df_check.columns:
                                # Multi-annotation extraction
                                for src, group in df_check.groupby("Annotation_Source"):
                                    valid_group = group[group["Term"] != "None"]
                                    if valid_group.empty: continue

                                    res_dict = {
                                        'n_bps': item['n_bps'],
                                        'iter': item['iter'],
                                        'n_samples': self._count_samples_from_file(item["expr"]),
                                        'source': src
                                    }
                                    if self.do_auroc and "AUC" in valid_group.columns:
                                        res_dict['auc'] = valid_group["AUC"].mean()
                                    if self.do_aupr and "AUPR" in valid_group.columns:
                                        res_dict['aupr'] = valid_group["AUPR"].mean()
                                    final_results.append(res_dict)
                                continue
                            else:
                                # Legacy fallback
                                valid = True
                                auc_val, aupr_val = None, None
                                if self.do_auroc:
                                    if "AUC" not in df_check.columns:
                                        valid = False
                                    else:
                                        auc_val = df_check["AUC"].mean()
                                if self.do_aupr:
                                    if "AUPR" not in df_check.columns:
                                        valid = False
                                    else:
                                        aupr_val = df_check["AUPR"].mean()

                                if valid:
                                    res_dict = {
                                        'n_bps': item['n_bps'], 'iter': item['iter'],
                                        'n_samples': self._count_samples_from_file(item["expr"]),
                                        'source': "Legacy"
                                    }
                                    if self.do_auroc: res_dict['auc'] = auc_val
                                    if self.do_aupr: res_dict['aupr'] = aupr_val
                                    final_results.append(res_dict)
                                    continue
                        else:
                            print(f"[Warn] Output {item['out']} is missing requested metrics or empty. Re-queueing.")
                    except Exception as e:
                        print(f"[Warn] Corrupted EGAD output at {item['out']} ({e}). Re-queueing.")

                # 2. THEN, if we actually need to calculate it, check if the network exists
                if not item["net"].exists():
                    continue  # Can't evaluate a network that failed to generate

                # 3. If we made it here, the output is missing/bad AND the network exists. Queue it.
                egad_tasks.append(item)

            if egad_tasks:
                with tqdm(total=len(egad_queue), initial=len(egad_queue) - len(egad_tasks), desc="EGAD Eval") as p_egad:
                    with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                        f_map = {
                            executor.submit(
                                run_egad_task,
                                network_file=i["net"],
                                mapman_file=self.mapman_file,
                                go_file=self.go_file,
                                out_file=i["out"],
                                script_path=r_script,
                                log_path=i["dir"] / "egad.log",
                                curves_prefix=None,  # Explicitly no curves requested for saturation
                                do_auroc=self.do_auroc,
                                do_aupr=self.do_aupr
                            ): i for i in egad_tasks
                        }

                        for f in concurrent.futures.as_completed(f_map):
                            it = f_map[f]
                            results_dict = f.result()

                            if results_dict:
                                for src, mets in results_dict.items():
                                    res_dict = {
                                        'n_bps': it['n_bps'],
                                        'iter': it['iter'],
                                        'n_samples': self._count_samples_from_file(it["expr"]),
                                        'source': src
                                    }
                                    if self.do_auroc and mets.get('auc') is not None: res_dict['auc'] = mets['auc']
                                    if self.do_aupr and mets.get('aupr') is not None: res_dict['aupr'] = mets['aupr']
                                    final_results.append(res_dict)
                            p_egad.update(1)

            self._generate_plots(final_results)