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
            # Using pyarrow for speed, falling back to standard if needed
            df = pd.read_csv(self.vst_path, sep="\t", index_col=0, engine="pyarrow")
            self.master_df = df
            return self.master_df
        except:
            df = pd.read_csv(self.vst_path, sep="\t", index_col=0)
            self.master_df = df
            return df

    def _write_seidr_files(self, df: pd.DataFrame, out_dir: Path):
        """Restored original logic: df is Sample x Gene."""
        df = df.fillna(0)
        g_file, e_file = out_dir / "genes.txt", out_dir / "expression.tsv"
        with open(g_file, "w") as f:
            f.write("\n".join(df.columns.astype(str)))
        # Write expression: No header, no index
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

    def run(self):
        # --- PHASE 0: CLEANUP ---
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

        # --- PHASE 1: PREP METADATA ---
        # Get sample IDs directly from the BioProject objects in the dataset
        bp_meta = []
        for bp in self.dataset.bioprojects:
            s_ids = bp.get_sample_ids()
            if s_ids:
                bp_meta.append({'obj': bp, 'size': len(s_ids), 'ids': s_ids})

        total_bps = len(bp_meta)
        if total_bps == 0:
            print("[Error] No BioProjects found in dataset.")
            return

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

        # --- PHASE 2: BATCH GENERATION ---
        seidr_queue, egad_queue = [], []
        global_genes = self.deseq2_dir / "genes.txt"
        global_expr = self.deseq2_dir / "expression.tsv"
        matrix = self._get_master_df()

        print(f"\n[Phase 1] Generating Batches...")
        print(f"\n[Phase 1] Generating Batches...")
        with tqdm(total=len(work_definitions), desc="Batch Gen", unit="batch") as pbar:
            for n_bps, iter_num, s_name in work_definitions:
                iter_dir = self.base_outdir / s_name / f"iter{iter_num}"
                iter_dir.mkdir(parents=True, exist_ok=True)

                # Setup Paths
                if n_bps == total_bps:
                    g_path, e_path = global_genes, global_expr
                else:
                    g_path, e_path = iter_dir / "genes.txt", iter_dir / "expression.tsv"

                seidr_queue.append((g_path, e_path, iter_dir, s_name, iter_num))
                egad_queue.append({
                    "n_bps": n_bps, "iter": iter_num, "dir": iter_dir,
                    "net": iter_dir / "network_saturation_edges.tsv",
                    "out": iter_dir / "egad_auroc.tsv", "expr": e_path
                })

                # Check if we need to do the heavy lifting
                if n_bps < total_bps:
                    # REMOVED the pbar.update(1) from here. We only need one at the bottom.
                    if not self.force and g_path.exists() and e_path.exists():
                        pass
                    else:
                        # --- GREEDY OPTIMIZATION FOR SAMPLES ---
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
                        for item in current_selection:
                            selected_samples.extend(item['ids'])

                        final_s = [s for s in selected_samples if not matrix.empty and s in matrix.index]
                        if final_s:
                            sub = matrix.loc[final_s].copy().fillna(0)
                            sub = sub.loc[:, sub.var(axis=0) > 0.1]
                            if not sub.empty:
                                self._write_seidr_files(sub, iter_dir)

                # The ONE AND ONLY update per loop iteration
                pbar.update(1)

        # --- PHASE 3: SEIDR INFERENCE ---
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

        # --- PHASE 4: EGAD EVALUATIONS ---
        if self.mapman_file:
            print(f"\n[Phase 3] EGAD Analysis...")
            r_script = get_egad_script_path()
            final_results = []
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