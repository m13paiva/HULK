import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, List, Tuple, TYPE_CHECKING, Dict
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

        # FIXED: Directory is now strictly 'saturation'
        self.base_outdir = self.config.shared / "saturation"
        self.base_outdir.mkdir(parents=True, exist_ok=True)

        self.iterations = 10
        self.vst_path = self.config.shared / "deseq2" / "vst.tsv"
        self.master_df = None

    def _get_master_df(self) -> pd.DataFrame:
        if self.master_df is not None: return self.master_df
        if not self.vst_path.exists():
            print(f"[Error] VST matrix missing: {self.vst_path}")
            return pd.DataFrame()
        try:
            try:
                df = pd.read_csv(self.vst_path, sep="\t", index_col=0, engine="pyarrow")
            except:
                df = pd.read_csv(self.vst_path, sep="\t", index_col=0)
            self.master_df = df
            return self.master_df
        except Exception:
            return pd.DataFrame()

    def _write_seidr_files(self, df: pd.DataFrame, out_dir: Path):
        df = df.fillna(0)
        g_file = out_dir / "genes.txt"
        e_file = out_dir / "expression.tsv"

        with open(g_file, "w") as f: f.write("\n".join(df.index.astype(str)))
        df.T.to_csv(e_file, sep="\t", header=False, index=False, float_format='%.6f')
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

        steps = []
        current = 0
        for inc in increments:
            current += inc
            steps.append(current)

        if steps[-1] != total_bps: steps[-1] = total_bps
        return sorted(list(set(steps)))

    def _generate_plots(self, results: List[dict]):
        if not results: return

        df = pd.DataFrame(results)
        df.to_csv(self.base_outdir / "saturation_results_raw.tsv", sep="\t", index=False)

        agg = df.groupby("n_bps").agg(
            auc_mean=("auc", "mean"),
            auc_std=("auc", "std"),
            samples_mean=("n_samples", "mean"),
            samples_std=("n_samples", "std"),
            count=("auc", "count")
        ).reset_index()
        agg = agg.fillna(0)
        agg.to_csv(self.base_outdir / "saturation_results_summary.tsv", sep="\t", index=False)

        color_auc = '#2c3e50'
        color_err = '#e74c3c'
        color_samp = '#27ae60'

        try:
            plt.figure(figsize=(8, 6))
            plt.plot(agg["n_bps"], agg["auc_mean"], color=color_auc, linewidth=2, zorder=1)
            if (agg["auc_std"] > 0).any():
                plt.errorbar(agg["n_bps"], agg["auc_mean"], yerr=agg["auc_std"],
                             fmt='o', color=color_auc, ecolor=color_err,
                             elinewidth=2, capsize=5, linestyle='None', zorder=2)
            else:
                plt.scatter(agg["n_bps"], agg["auc_mean"], color=color_auc, zorder=2)

            plt.title("Network Reconstruction Performance vs BioProject Count")
            plt.xlabel("Number of BioProjects")
            plt.ylabel("Mean AUROC")
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.tight_layout()
            plt.savefig(self.base_outdir / "saturation_auroc.pdf")
            plt.close()
        except Exception as e:
            print(f"[Plot] Error 1: {e}")

        try:
            plt.figure(figsize=(8, 6))
            plt.plot(agg["n_bps"], agg["samples_mean"], color=color_samp, linewidth=2, zorder=1)
            if (agg["samples_std"] > 0).any():
                plt.errorbar(agg["n_bps"], agg["samples_mean"], yerr=agg["samples_std"],
                             fmt='o', color=color_samp, ecolor='black',
                             elinewidth=2, capsize=5, linestyle='None', zorder=2)
            else:
                plt.scatter(agg["n_bps"], agg["samples_mean"], color=color_samp, zorder=2)

            plt.title("Sample Accumulation vs BioProject Count")
            plt.xlabel("Number of BioProjects")
            plt.ylabel("Total Samples")
            plt.grid(True, linestyle='--', alpha=0.5)
            plt.tight_layout()
            plt.savefig(self.base_outdir / "saturation_samples.pdf")
            plt.close()
        except Exception as e:
            print(f"[Plot] Error 2: {e}")

        try:
            fig, ax1 = plt.subplots(figsize=(10, 6))
            ax1.set_xlabel('Number of BioProjects')
            ax1.set_ylabel('Mean AUROC', color=color_auc)
            ax1.plot(agg["n_bps"], agg["auc_mean"], color=color_auc, linewidth=2, label="AUROC")
            if (agg["auc_std"] > 0).any():
                ax1.errorbar(agg["n_bps"], agg["auc_mean"], yerr=agg["auc_std"],
                             fmt='o', color=color_auc, ecolor=color_err,
                             elinewidth=2, capsize=4, linestyle='None')
            ax1.tick_params(axis='y', labelcolor=color_auc)
            ax1.grid(True, linestyle='--', alpha=0.3)

            ax2 = ax1.twinx()
            ax2.set_ylabel('Mean Samples', color=color_samp)
            ax2.plot(agg["n_bps"], agg["samples_mean"], color=color_samp, linewidth=2, linestyle='--', label="Samples")
            if (agg["samples_std"] > 0).any():
                ax2.errorbar(agg["n_bps"], agg["samples_mean"], yerr=agg["samples_std"],
                             fmt='o', color=color_samp, ecolor='black',
                             elinewidth=1, capsize=3, linestyle='None', alpha=0.7)
            ax2.tick_params(axis='y', labelcolor=color_samp)

            plt.title("Saturation Analysis: Performance vs Data Size")
            plt.tight_layout()
            plt.savefig(self.base_outdir / "saturation_combined.pdf")
            plt.close()
        except Exception as e:
            print(f"[Plot] Error 3: {e}")

        print(f"\n[Saturation] Plots saved to {self.base_outdir}")

    def run(self):
        # --- PHASE 0: CLEANUP ---
        if self.force and self.base_outdir.exists():
            print(f"[Saturation] Force=True: Purging all contents from {self.base_outdir}...")
            for item in self.base_outdir.iterdir():
                try:
                    if item.is_dir():
                        shutil.rmtree(item)
                    else:
                        item.unlink()
                except Exception as e:
                    print(f"[Warn] Failed to delete {item}: {e}")

        # --- PHASE 1: COLLECT BPs ---
        bp_map: Dict[str, List[str]] = {}
        for bp in self.dataset.bioprojects:
            s_ids = [s.id for s in bp.samples]
            if s_ids: bp_map[bp.id] = s_ids

        bp_list = list(bp_map.keys())
        total_bps = len(bp_list)
        if total_bps == 0: return

        all_samples_count = sum(len(v) for v in bp_map.values())
        avg_samples_per_bp = all_samples_count / total_bps

        # --- PHASE 2: DEFINE STEPS ---
        bp_steps = self._calculate_steps(total_bps)
        step_name_map = {n: f"step{i + 1}" for i, n in enumerate(bp_steps)}
        print(f"\n[Saturation] BioProject Steps (Total {total_bps} BPs): {list(step_name_map.values())}")

        work_definitions = []
        for n_bps in bp_steps:
            iters = 1 if n_bps == total_bps else self.iterations
            s_name = step_name_map[n_bps]
            for i in range(1, iters + 1):
                work_definitions.append((n_bps, i, s_name))

        total_tasks = len(work_definitions)

        # --- PHASE 3: PRE-CHECK ---
        if not self.force and self.mapman_file:
            all_done = True
            res = []
            for n_bps, i, s_name in work_definitions:
                d = self.base_outdir / s_name / f"iter{i}"
                if not (d / "egad_auroc.tsv").exists():
                    all_done = False;
                    break

                e_file = (self.config.shared / "deseq2" / "expression.tsv"
                          if n_bps == total_bps else d / "expression.tsv")

                try:
                    auc = pd.read_csv(d / "egad_auroc.tsv", sep="\t")["AUC"].mean()
                    n_samp = self._count_samples_from_file(e_file)
                    res.append({'n_bps': n_bps, 'iter': i, 'auc': auc, 'n_samples': n_samp})
                except:
                    pass

            if all_done:
                print(f"[Saturation] All tasks completed. Regenerating plots...")
                self._generate_plots(res)
                return

        preset = getattr(self.config, "seidr_preset", "FAST")

        # --- PHASE 4: GENERATE INPUTS ---
        print(f"\n[Phase 1] Generating Batches...")
        seidr_queue = []
        egad_queue = []

        from collections import defaultdict
        step_groups = defaultdict(list)
        for w in work_definitions: step_groups[w[0]].append(w)

        pbar = tqdm(total=total_tasks, desc="Batch Gen", unit="batch")
        sorted_steps = sorted(step_groups.keys())

        global_genes = self.config.shared / "deseq2" / "genes.txt"
        global_expr = self.config.shared / "deseq2" / "expression.tsv"

        for n_bps in sorted_steps:
            group = step_groups[n_bps]
            s_name = step_name_map[n_bps]

            if n_bps == total_bps:
                for (curr_n_bps, iter_num, _) in group:
                    iter_dir = self.base_outdir / s_name / f"iter{iter_num}"
                    iter_dir.mkdir(parents=True, exist_ok=True)

                    net_file = iter_dir / "network_saturation_edges.tsv"
                    egad_file = iter_dir / "egad_auroc.tsv"

                    seidr_queue.append((global_genes, global_expr, iter_dir, s_name, iter_num))
                    egad_queue.append({
                        "n_bps": curr_n_bps, "iter": iter_num, "step": s_name,
                        "dir": iter_dir, "net": net_file, "out": egad_file,
                        "expr": global_expr
                    })
                    pbar.update(1)
                continue

            current_pool = list(bp_list)
            target_samples = n_bps * avg_samples_per_bp

            for (curr_n_bps, iter_num, _) in group:
                if len(current_pool) < curr_n_bps: current_pool = list(bp_list)

                best_combo = None;
                best_diff = float('inf')
                for _ in range(50):
                    cand = random.sample(current_pool, curr_n_bps)
                    diff = abs(sum(len(bp_map[b]) for b in cand) - target_samples)
                    if diff < best_diff:
                        best_diff = diff;
                        best_combo = cand
                        if diff < (target_samples * 0.05): break

                selected_bps = best_combo
                for b in selected_bps:
                    if b in current_pool: current_pool.remove(b)

                final_samples = []
                for b in selected_bps: final_samples.extend(bp_map[b])

                iter_dir = self.base_outdir / s_name / f"iter{iter_num}"
                iter_dir.mkdir(parents=True, exist_ok=True)

                net_file = iter_dir / "network_saturation_edges.tsv"
                g_file = iter_dir / "genes.txt"
                e_file = iter_dir / "expression.tsv"
                egad_file = iter_dir / "egad_auroc.tsv"

                seidr_queue.append((g_file, e_file, iter_dir, s_name, iter_num))
                egad_queue.append({
                    "n_bps": curr_n_bps, "iter": iter_num, "step": s_name,
                    "dir": iter_dir, "net": net_file, "out": egad_file,
                    "expr": e_file
                })

                if not self.force:
                    if net_file.exists() or (g_file.exists() and e_file.exists()):
                        pbar.update(1);
                        continue

                matrix = self._get_master_df()
                if matrix.empty: pbar.update(1); continue

                valid_samples = [s for s in final_samples if s in matrix.columns]
                if not valid_samples: pbar.update(1); continue

                sub = matrix[valid_samples].copy()
                sub.fillna(0, inplace=True)
                sub = sub[sub.var(axis=1) > 0.1]

                if sub.empty: pbar.update(1); continue
                self._write_seidr_files(sub, iter_dir)
                pbar.update(1)
        pbar.close()

        # --- PHASE 5: SEIDR ---
        print(f"\n[Phase 2] Seidr Inference...")
        tasks_to_run = []
        existing_nets = 0
        for item in seidr_queue:
            d = item[2]
            if not self.force and (d / ".saturation.done").exists():
                existing_nets += 1
            else:
                tasks_to_run.append(item)

        if tasks_to_run:
            threads_per_job = max(1, self.total_thread_budget // self.workers)
            with tqdm(total=len(seidr_queue), initial=existing_nets, desc="Seidr Inf", unit="net") as pbar:
                with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                    f_map = {
                        executor.submit(run_seidr_batch, cfg=self.config, genes_file=i[0], expression_file=i[1],
                                        outdir=i[2], preset=preset, threads=threads_per_job): i for i in tasks_to_run
                    }
                    for f in concurrent.futures.as_completed(f_map):
                        pbar.update(1)

        # --- PHASE 6: EGAD ---
        if self.mapman_file:
            print(f"\n[Phase 3] EGAD Analysis...")
            try:
                r_script = get_egad_script_path()
            except:
                return

            final_results = []
            egad_tasks = []
            existing_egad_count = 0

            for item in egad_queue:
                n_samp = self._count_samples_from_file(item["expr"])
                base_res = {'n_bps': item['n_bps'], 'iter': item['iter'], 'n_samples': n_samp}

                if not self.force and item["out"].exists():
                    try:
                        val = pd.read_csv(item["out"], sep="\t")["AUC"].mean()
                        base_res['auc'] = val
                        final_results.append(base_res)
                        existing_egad_count += 1
                        continue
                    except:
                        pass

                if item["net"].exists():
                    item['base_res'] = base_res
                    egad_tasks.append(item)

            if egad_tasks:
                with tqdm(total=len(egad_queue), initial=existing_egad_count, desc="EGAD Eval", unit="task") as pbar:
                    with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                        f_map = {
                            executor.submit(
                                run_egad_task, network_file=i["net"], mapman_file=self.mapman_file, out_file=i["out"],
                                script_path=r_script, log_path=i["dir"] / "egad.log"
                            ): i for i in egad_tasks
                        }
                        for f in concurrent.futures.as_completed(f_map):
                            item = f_map[f]
                            res = f.result()
                            if res is not None:
                                res_dict = item['base_res']
                                res_dict['auc'] = res
                                final_results.append(res_dict)
                            pbar.update(1)

            self._generate_plots(final_results)
        else:
            print("\n[Phase 3] Skipped (No MapMan).")