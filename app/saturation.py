import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, List, Tuple, TYPE_CHECKING
import concurrent.futures

from .seidr import run_seidr_batch

if TYPE_CHECKING:
    from .entities import Dataset, Config


class BatchOrchestrator:
    def __init__(self, dataset: "Dataset", config: "Config",
                 seed: Optional[int] = None,
                 workers: int = 4,
                 max_threads: Optional[int] = None):
        """
        Orchestrator: 'Knapsack' Mode + Parallel Execution.
        OPTIMIZED: Pre-transposes VST matrix to eliminate redundant Ops.
        """
        self.dataset = dataset
        self.config = config
        self.seed = seed
        self.workers = workers

        if max_threads is not None:
            self.total_thread_budget = max_threads
        else:
            self.total_thread_budget = self.config.max_threads

        if self.seed is not None:
            np.random.seed(self.seed)
            print(f"[Saturation] Random seed: {self.seed}")

        self.base_outdir = self.config.shared / "egad" / "data_size_analysis"
        self.base_outdir.mkdir(parents=True, exist_ok=True)

        self.iterations = 10
        self.vst_path = self.config.shared / "deseq2" / "vst.tsv"

        # LOAD AND TRANSPOSE IMMEDIATELY
        self.master_df = self._load_and_transpose_vst()

    def _load_and_transpose_vst(self) -> pd.DataFrame:
        if not self.vst_path.exists():
            print(f"[Error] VST matrix missing: {self.vst_path}")
            return pd.DataFrame()
        try:
            print(f">> Loading Master VST Matrix from {self.vst_path}...", flush=True)
            try:
                # Optimized read if pyarrow installed
                df = pd.read_csv(self.vst_path, sep="\t", index_col=0, engine="pyarrow")
            except:
                df = pd.read_csv(self.vst_path, sep="\t", index_col=0)

            # Transpose ONCE -> (Gene x Sample)
            # This allows fast column slicing in the loop
            df = df.T

            print(f"   [Info] Matrix ready. Shape: {df.shape[0]} genes x {df.shape[1]} samples.")
            return df
        except Exception as e:
            print(f"[Error] Failed to read VST: {e}")
            return pd.DataFrame()

    def _write_seidr_files(self, df: pd.DataFrame, out_dir: Path):
        # Input 'df' is (Gene x Sample), already clean and filtered.

        g_file = out_dir / "genes.txt"
        e_file = out_dir / "expression.tsv"

        # 1. Write Genes (Index of df)
        with open(g_file, "w") as f:
            f.write("\n".join(df.index.astype(str)))

        # 2. Transpose to (Sample x Gene) for Seidr tools
        #    Seidr expects samples as rows and genes as columns for correlation
        df_out = df.T

        # 3. Write Matrix (Values only, no headers/index)
        df_out.to_csv(e_file, sep="\t", header=False, index=False, float_format='%.6f')

        return g_file, e_file

    def run(self):
        if self.master_df.empty: return

        total_samples = len(self.dataset)
        steps_pct = [10, 20, 30, 40, 50, 60, 70, 80, 90]

        # Build Metadata
        # Master Matrix is (Gene x Sample), so we check columns for sample IDs
        valid_samples_in_matrix = set(self.master_df.columns)

        all_projects_meta = []
        for bp in self.dataset.bioprojects:
            s_ids = [s.id for s in bp.samples]
            valid_s_ids = [sid for sid in s_ids if sid in valid_samples_in_matrix]
            if valid_s_ids:
                all_projects_meta.append((bp, len(valid_s_ids), valid_s_ids))

        if not all_projects_meta:
            print("[Error] No valid BioProjects found.")
            return

        preset = getattr(self.config, "seidr_preset", "FAST")

        # --- PHASE 1: GENERATE FILES & QUEUE TASKS ---
        execution_queue = []
        print(f"\n[Phase 1] Generating inputs for {len(steps_pct)} batch steps x {self.iterations} iterations...",
              flush=True)

        for pct in steps_pct:
            n_target = int(total_samples * (pct / 100.0))
            if n_target == 0: continue

            dir_name = f"{pct}%_batch"
            available_pool = list(all_projects_meta)

            for i in range(1, self.iterations + 1):
                selected_samples = []
                current_n = 0
                used_in_this_run = set()

                while current_n < n_target:
                    remaining_space = n_target - current_n
                    candidates_whole = [
                        (bp, size, s_ids) for (bp, size, s_ids) in available_pool
                        if size <= remaining_space and bp not in used_in_this_run
                    ]

                    if candidates_whole:
                        weights = np.array([1.0 / (size ** 2) for (_, size, _) in candidates_whole])
                        weights /= weights.sum()
                        idx = np.random.choice(len(candidates_whole), p=weights)
                        chosen = candidates_whole[idx]
                        selected_samples.extend(chosen[2])
                        current_n += chosen[1]
                        used_in_this_run.add(chosen[0])
                        if chosen in available_pool: available_pool.remove(chosen)
                    else:
                        candidates_frag = [x for x in available_pool if x[0] not in used_in_this_run]
                        if not candidates_frag:
                            available_pool = list(all_projects_meta)
                            candidates_frag = [x for x in available_pool if x[0] not in used_in_this_run]
                            if not candidates_frag: break
                        idx = np.random.randint(len(candidates_frag))
                        chosen_bp, chosen_size, chosen_ids = candidates_frag[idx]
                        gap = n_target - current_n
                        take = min(gap, chosen_size)
                        subset = np.random.choice(chosen_ids, take, replace=False)
                        selected_samples.extend(subset)
                        current_n += take
                        used_in_this_run.add(chosen_bp)
                        item_to_remove = (chosen_bp, chosen_size, chosen_ids)
                        if item_to_remove in available_pool: available_pool.remove(item_to_remove)

                # --- FAST SLICING & FILTERING ---
                final_samples = [s for s in selected_samples if s in valid_samples_in_matrix]

                if final_samples:
                    iter_dir = self.base_outdir / dir_name / f"iter{i}"
                    iter_dir.mkdir(parents=True, exist_ok=True)

                    # 1. COLUMN SLICE (Subset Samples)
                    # Create explicit copy to handle NAs safely
                    subset_df = self.master_df[final_samples].copy()

                    # 2. FILL NA
                    subset_df.fillna(0, inplace=True)

                    # 3. VARIANCE FILTER
                    # Calculate variance across samples (axis=1 because subset_df is Gene x Sample)
                    variances = subset_df.var(axis=1)

                    # Drop genes with variance <= 0.1
                    subset_df = subset_df[variances > 0.1]

                    if subset_df.empty:
                        print(f"   [Warning] Batch {dir_name}/iter{i} dropped all genes due to low variance.")
                        continue

                    # 3. Write
                    g_path, e_path = self._write_seidr_files(subset_df, iter_dir)
                    execution_queue.append((g_path, e_path, iter_dir, dir_name, i))
                else:
                    print(f"   [Warning] Empty batch at {dir_name}/iter{i}")

        # --- PHASE 2: PARALLEL EXECUTION ---
        total_tasks = len(execution_queue)
        if total_tasks == 0:
            print("[Info] No tasks generated.")
            return

        threads_per_job = max(1, self.total_thread_budget // self.workers)

        print(f"\n[Phase 2] Executing {total_tasks} Seidr runs in parallel.")
        print(f"          Total Budget: {self.total_thread_budget} threads")
        print(f"          Workers: {self.workers} | Threads per Job: {threads_per_job}")
        print(f"          Preset: {preset}", flush=True)

        with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
            future_to_task = {}
            for (g, e, d, bname, iter_n) in execution_queue:
                fut = executor.submit(
                    run_seidr_batch,
                    cfg=self.config,
                    genes_file=g,
                    expression_file=e,
                    outdir=d,
                    preset=preset,
                    threads=threads_per_job
                )
                future_to_task[fut] = f"{bname}/iter{iter_n}"

            for future in concurrent.futures.as_completed(future_to_task):
                task_name = future_to_task[future]
                try:
                    future.result()
                except Exception as exc:
                    print(f"   [Error] {task_name} generated an exception: {exc}")