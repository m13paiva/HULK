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
                 force: bool = False,
                 num_steps: int = 10):
        self.dataset = dataset
        self.config = config
        self.seed = seed
        self.workers = workers
        self.mapman_file = mapman_file
        self.force = force
        self.num_steps = num_steps

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

        agg = df.groupby("n_bps").agg(
            auc_mean=("auc", "mean"),
            auc_std=("auc", "std"),
            samples_mean=("n_samples", "mean"),
            samples_std=("n_samples", "std")
        ).reset_index().fillna(0)

        agg.to_csv(self.base_outdir / "saturation_results_summary.tsv", sep="\t", index=False)

        color_auc, color_err, color_samp = '#2c3e50', '#e74c3c', '#27ae60'

        try:
            plt.figure(figsize=(8, 6))
            plt.plot(agg["n_bps"], agg["auc_mean"], color=color_auc, linewidth=2, marker='o')
            plt.errorbar(agg["n_bps"], agg["auc_mean"], yerr=agg["auc_std"], fmt='none', ecolor=color_err, capsize=3)
            plt.title("Performance vs BioProject Count")
            plt.xlabel("Number of BioProjects")
            plt.ylabel("Mean AUROC")
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.tight_layout()
            plt.savefig(self.base_outdir / "saturation_auroc.pdf")
            plt.close()
        except Exception as e:
            print(f"[Plot Error 1] {e}")

        try:
            plt.figure(figsize=(8, 6))
            plt.plot(agg["n_bps"], agg["samples_mean"], color=color_samp, linewidth=2, marker='s')
            plt.errorbar(agg["n_bps"], agg["samples_mean"], yerr=agg["samples_std"], fmt='none', ecolor='black',
                         capsize=3)
            plt.title("Sample Accumulation vs BioProject Count")
            plt.xlabel("Number of BioProjects")
            plt.ylabel("Total Samples")
            plt.grid(True, linestyle='--', alpha=0.6)
            plt.tight_layout()
            plt.savefig(self.base_outdir / "saturation_samples.pdf")
            plt.close()
        except Exception as e:
            print(f"[Plot Error 2] {e}")

        try:
            fig, ax1 = plt.subplots(figsize=(10, 6))
            ax1.set_xlabel('Number of BioProjects')
            ax1.set_ylabel('Mean AUROC', color=color_auc)
            ax1.plot(agg["n_bps"], agg["auc_mean"], color=color_auc, linewidth=2, marker='o', label="AUROC")
            ax1.tick_params(axis='y', labelcolor=color_auc)
            ax2 = ax1.twinx()
            ax2.set_ylabel('Mean Samples', color=color_samp)
            ax2.plot(agg["n_bps"], agg["samples_mean"], color=color_samp, linewidth=2, linestyle='--', marker='s',
                     label="Samples")
            ax2.tick_params(axis='y', labelcolor=color_samp)
            plt.title("Saturation Analysis Summary")
            fig.tight_layout()
            plt.savefig(self.base_outdir / "saturation_combined.pdf")
            plt.close()
        except Exception as e:
            print(f"[Plot Error 3] {e}")
        print(f"\n[Saturation] Plots updated in {self.base_outdir}")

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
                    "out": iter_dir / "egad_auroc.tsv", "expr": e_path
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
        if self.mapman_file:
            print(f"\n[Phase 3] EGAD Analysis...")
            r_script = get_egad_script_path()
            egad_tasks = [item for item in egad_queue if
                          item["net"].exists() and (self.force or not item["out"].exists())]
            for item in egad_queue:
                if not self.force and item["out"].exists():
                    try:
                        auc = pd.read_csv(item["out"], sep="\t")["AUC"].mean()
                        final_results.append({'n_bps': item['n_bps'], 'iter': item['iter'], 'auc': auc,
                                              'n_samples': self._count_samples_from_file(item["expr"])})
                    except:
                        pass
            if egad_tasks:
                with tqdm(total=len(egad_queue), initial=len(egad_queue) - len(egad_tasks), desc="EGAD Eval") as p_egad:
                    with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                        f_map = {executor.submit(run_egad_task, i["net"], self.mapman_file, i["out"], r_script,
                                                 i["dir"] / "egad.log"): i for i in egad_tasks}
                        for f in concurrent.futures.as_completed(f_map):
                            it, res = f_map[f], f.result()
                            if res is not None:
                                final_results.append({'n_bps': it['n_bps'], 'iter': it['iter'], 'auc': res,
                                                      'n_samples': self._count_samples_from_file(it["expr"])})
                            p_egad.update(1)
            self._generate_plots(final_results)