import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, List, Tuple, TYPE_CHECKING
import concurrent.futures

from .seidr import run_seidr_batch
from .egad import run_egad_task, aggregate_and_plot_egad, get_egad_script_path

if TYPE_CHECKING:
    from .entities import Dataset, Config


class BatchOrchestrator:
    def __init__(self, dataset: "Dataset", config: "Config",
                 seed: Optional[int] = None,
                 workers: int = 4,
                 max_threads: Optional[int] = None,
                 mapman_file: Optional[Path] = None,
                 force: bool = False):
        self.dataset = dataset
        self.config = config
        self.seed = seed
        self.workers = workers
        self.mapman_file = mapman_file
        self.force = force

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

        # LAZY LOAD STATE
        self.master_df = None

    def _get_master_df(self) -> pd.DataFrame:
        """
        Lazy loader: Only reads VST if we actually need to GENERATE a file.
        """
        if self.master_df is not None:
            return self.master_df

        if not self.vst_path.exists():
            print(f"[Error] VST matrix missing: {self.vst_path}")
            return pd.DataFrame()

        try:
            print(f">> Loading Master VST Matrix from {self.vst_path}...", flush=True)
            try:
                df = pd.read_csv(self.vst_path, sep="\t", index_col=0, engine="pyarrow")
            except:
                df = pd.read_csv(self.vst_path, sep="\t", index_col=0)

            # Transpose to (Gene x Sample) for fast column slicing
            self.master_df = df.T
            print(
                f"   [Info] Matrix ready. Shape: {self.master_df.shape[0]} genes x {self.master_df.shape[1]} samples.")
            return self.master_df
        except Exception as e:
            print(f"[Error] Failed to read VST: {e}")
            return pd.DataFrame()

    def _write_seidr_files(self, df: pd.DataFrame, out_dir: Path):
        df = df.fillna(0)
        g_file = out_dir / "genes.txt"
        e_file = out_dir / "expression.tsv"
        with open(g_file, "w") as f:
            f.write("\n".join(df.index.astype(str)))
        df.to_csv(e_file, sep="\t", header=False, index=False, float_format='%.6f')
        return g_file, e_file

    def run(self):
        total_samples = len(self.dataset)
        steps_pct = [10, 20, 30, 40, 50, 60, 70, 80, 90]

        # We assume all samples in dataset are valid candidates.
        # We DO NOT check against VST index here to avoid I/O.
        # If a sample is missing from VST, it will fail only during the WRITE step (if reached).
        all_projects_meta = []
        for bp in self.dataset.bioprojects:
            s_ids = [s.id for s in bp.samples]
            if s_ids:
                all_projects_meta.append((bp, len(s_ids), s_ids))

        if not all_projects_meta:
            print("[Error] No BioProjects found in dataset.")
            return

        preset = getattr(self.config, "seidr_preset", "FAST")

        # --- PHASE 1: GENERATE FILES & QUEUE ---
        seidr_queue = []
        egad_queue = []

        print(f"\n[Phase 1] Scanning {len(steps_pct)} batch steps x {self.iterations} iterations...", flush=True)

        for pct in steps_pct:
            n_target = int(total_samples * (pct / 100.0))
            if n_target == 0: continue

            dir_name = f"{pct}%_batch"
            # We must maintain the pool logic even if files exist,
            # so the random seed sequence remains deterministic across runs.
            available_pool = list(all_projects_meta)

            for i in range(1, self.iterations + 1):
                # --- START KNAPSACK LOGIC (CPU ONLY - FAST) ---
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
                # --- END KNAPSACK LOGIC ---

                # --- CHECK FILES ---
                iter_dir = self.base_outdir / dir_name / f"iter{i}"
                iter_dir.mkdir(parents=True, exist_ok=True)

                # Check directly for the final network.
                # If it exists, we don't care about VST or genes.txt or expression.tsv
                network_file = iter_dir / "network_saturation_edges.tsv"

                if not self.force and network_file.exists():
                    # FAST PATH: Network exists. Queue directly for EGAD.
                    egad_queue.append({
                        "net": network_file,
                        "out": iter_dir / "egad_auroc.tsv",
                        "dir": iter_dir,
                        "pct": pct,
                        "iter": i
                    })
                    # Skip Seidr queue. Skip VST load.
                    continue

                # --- MISSING FILES PATH ---
                # Network missing. Must queue Seidr.
                # Do we have input files?
                g_path = iter_dir / "genes.txt"
                e_path = iter_dir / "expression.tsv"

                if not (g_path.exists() and e_path.exists()):
                    # Inputs missing. NOW we must load VST.
                    matrix = self._get_master_df()
                    if matrix.empty: continue

                    # Filter samples that actually exist in the matrix
                    final_samples = [s for s in selected_samples if s in matrix.columns]

                    if not final_samples:
                        print(f"   [Warning] Empty batch at {dir_name}/iter{i}")
                        continue

                    subset_df = matrix[final_samples].copy()
                    subset_df.fillna(0, inplace=True)
                    variances = subset_df.var(axis=1)
                    subset_df = subset_df[variances > 0.1]

                    if subset_df.empty:
                        print(f"   [Warning] Batch {dir_name}/iter{i} dropped all genes.")
                        continue

                    self._write_seidr_files(subset_df, iter_dir)

                # Queue for Seidr
                seidr_queue.append((g_path, e_path, iter_dir, dir_name, i))

                # Queue for EGAD (will run after Seidr)
                egad_queue.append({
                    "net": network_file,
                    "out": iter_dir / "egad_auroc.tsv",
                    "dir": iter_dir,
                    "pct": pct,
                    "iter": i
                })

        # --- PHASE 2: SEIDR ---
        if seidr_queue:
            threads_per_job = max(1, self.total_thread_budget // self.workers)
            print(f"\n[Phase 2] Executing {len(seidr_queue)} Seidr runs in parallel.")
            print(f"          Total Budget: {self.total_thread_budget} threads")
            print(f"          Workers: {self.workers} | Threads per Job: {threads_per_job}")

            with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                future_to_task = {}
                for (g, e, d, bname, iter_n) in seidr_queue:
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
                    try:
                        future.result()
                    except Exception as exc:
                        print(f"   [Error] Seidr {future_to_task[future]}: {exc}")
        else:
            print("\n[Phase 2] All networks found. Skipping Seidr inference.")

        # --- PHASE 3: EGAD ---
        if self.mapman_file:
            print(f"\n[Phase 3] Running EGAD Analysis (Workers: {self.workers})...", flush=True)

            try:
                r_script = get_egad_script_path()
            except Exception as e:
                print(f"[Error] Failed to locate EGAD script: {e}")
                return

            egad_results = []

            with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                futures = {}
                for item in egad_queue:
                    # Double check network exists (if it failed in Phase 2)
                    if not item["net"].exists():
                        print(f"   [Skip] Network missing for {item['pct']}% iter {item['iter']}")
                        continue

                    if not self.force and item["out"].exists():
                        try:
                            val = pd.read_csv(item["out"], sep="\t")["AUC"].mean()
                            egad_results.append({'pct': item['pct'], 'iter': item['iter'], 'auc': val})
                        except:
                            pass
                        continue

                    f = executor.submit(
                        run_egad_task,
                        network_file=item["net"],
                        mapman_file=self.mapman_file,
                        out_file=item["out"],
                        script_path=r_script,
                        log_path=item["dir"] / "egad.log"
                    )
                    futures[f] = item

                for future in concurrent.futures.as_completed(futures):
                    item = futures[future]
                    res_auc = future.result()
                    if res_auc is not None:
                        egad_results.append({'pct': item['pct'], 'iter': item['iter'], 'auc': res_auc})
                    else:
                        print(f"   [Error] EGAD failed for {item['pct']}% iter {item['iter']}")

            aggregate_and_plot_egad(egad_results, self.base_outdir)
        else:
            print("\n[Phase 3] Skipped (No MapMan file provided).")