import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, List, Tuple, TYPE_CHECKING
import concurrent.futures
from tqdm import tqdm

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
            self.master_df = df.T
            return self.master_df
        except Exception:
            return pd.DataFrame()

    def _write_seidr_files(self, df: pd.DataFrame, out_dir: Path):
        df = df.fillna(0)
        g_file = out_dir / "genes.txt"
        e_file = out_dir / "expression.tsv"
        with open(g_file, "w") as f: f.write("\n".join(df.index.astype(str)))
        df.to_csv(e_file, sep="\t", header=False, index=False, float_format='%.6f')
        return g_file, e_file

    def run(self):
        total_samples = len(self.dataset)
        steps_pct = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

        all_projects_meta = []
        for bp in self.dataset.bioprojects:
            s_ids = [s.id for s in bp.samples]
            if s_ids: all_projects_meta.append((bp, len(s_ids), s_ids))

        if not all_projects_meta:
            print("[Error] No BioProjects found.")
            return

        preset = getattr(self.config, "seidr_preset", "FAST")

        # --- PREPARE WORK DEFINITIONS ---
        # We flatten the loops first so we can count total tasks
        work_definitions = []
        for pct in steps_pct:
            iters = 1 if pct == 100 else self.iterations
            for i in range(1, iters + 1):
                work_definitions.append((pct, i))

        total_tasks = len(work_definitions)

        # ==============================================================================
        # PHASE 1: GENERATION (SCAN -> EXECUTE)
        # ==============================================================================
        print(f"\n[Phase 1] Scanning for existing batch inputs...")

        # 1. Scan for existing inputs to set progress bar offset
        existing_batches = 0
        if not self.force:
            for pct, i in work_definitions:
                d = self.base_outdir / f"{pct}%_batch" / f"iter{i}"
                # If network exists, we assume inputs were good or processed.
                # If inputs exist, we are good.
                if (d / "network_saturation_edges.tsv").exists() or \
                        (d / ".saturation.done").exists() or \
                        ((d / "genes.txt").exists() and (d / "expression.tsv").exists()):
                    existing_batches += 1

        print(f"          Found {existing_batches}/{total_tasks} existing batches.")

        seidr_queue = []
        egad_queue = []

        # 2. Execute with Offset
        with tqdm(total=total_tasks, initial=existing_batches, desc="Batch Gen", unit="batch") as pbar:
            for pct, i in work_definitions:
                dir_name = f"{pct}%_batch"

                # --- A. RNG LOGIC (MUST RUN ALWAYS for Seed Consistency) ---
                selected_samples = []
                if pct == 100:
                    for bp in self.dataset.bioprojects: selected_samples.extend([s.id for s in bp.samples])
                else:
                    # Knapsack Logic
                    # We must reconstruct 'available_pool' fresh for each step or maintain it?
                    # The original logic reset available_pool per 'pct' loop, but maintained it across 'i'.
                    # We must mimic the exact nested loop structure for RNG state if we flattened it.
                    # Wait, flattened loop breaks the per-pct pool reset if we aren't careful.
                    # FIX: Re-initialize pool if 'i' == 1
                    if i == 1:
                        available_pool = list(all_projects_meta)

                    n_target = int(total_samples * (pct / 100.0))
                    if n_target == 0:
                        # Edge case: skip but counts as task?
                        # If n_target=0, previous logic continued.
                        # We should mark as done or skip logic.
                        pbar.update(1) if (pct, i) not in [] else None  # Simplify
                        continue

                    current_n = 0;
                    used_in_this_run = set()
                    while current_n < n_target:
                        remaining = n_target - current_n
                        cands_whole = [(bp, sz, sids) for (bp, sz, sids) in available_pool if
                                       sz <= remaining and bp not in used_in_this_run]
                        if cands_whole:
                            w = np.array([1.0 / (sz ** 2) for (_, sz, _) in cands_whole]);
                            w /= w.sum()
                            chosen = cands_whole[np.random.choice(len(cands_whole), p=w)]
                            selected_samples.extend(chosen[2]);
                            current_n += chosen[1];
                            used_in_this_run.add(chosen[0])
                            if chosen in available_pool: available_pool.remove(chosen)
                        else:
                            cands_frag = [x for x in available_pool if x[0] not in used_in_this_run]
                            if not cands_frag:
                                available_pool = list(all_projects_meta)
                                cands_frag = [x for x in available_pool if x[0] not in used_in_this_run]
                                if not cands_frag: break
                            chosen = cands_frag[np.random.randint(len(cands_frag))]
                            take = min(n_target - current_n, chosen[1])
                            selected_samples.extend(np.random.choice(chosen[2], take, replace=False))
                            current_n += take;
                            used_in_this_run.add(chosen[0])
                            if chosen in available_pool: available_pool.remove(chosen)
                # --- END RNG ---

                # --- B. DEFINE PATHS & QUEUE ---
                iter_dir = self.base_outdir / dir_name / f"iter{i}"
                iter_dir.mkdir(parents=True, exist_ok=True)

                net_file = iter_dir / "network_saturation_edges.tsv"
                done_marker = iter_dir / ".saturation.done"
                egad_file = iter_dir / "egad_auroc.tsv"
                g_path = iter_dir / "genes.txt"
                e_path = iter_dir / "expression.tsv"

                # Add to queues regardless of existence (Phase 2/3 need them)
                seidr_queue.append((g_path, e_path, iter_dir, dir_name, i))
                egad_queue.append({"net": net_file, "out": egad_file, "dir": iter_dir, "pct": pct, "iter": i})

                # --- C. CHECK EXISTENCE ---
                if not self.force:
                    # If network or inputs exist, we count this as "Existing" (offset covered it)
                    if net_file.exists() or done_marker.exists() or (g_path.exists() and e_path.exists()):
                        continue

                # --- D. GENERATE (If missing) ---
                matrix = self._get_master_df()
                if matrix.empty:
                    pbar.update(1)
                    continue

                final_s = [s for s in selected_samples if s in matrix.columns]
                if not final_s:
                    pbar.update(1)
                    continue

                sub = matrix[final_s].copy();
                sub.fillna(0, inplace=True)
                sub = sub[sub.var(axis=1) > 0.1]

                if sub.empty:
                    pbar.update(1)
                    continue

                self._write_seidr_files(sub, iter_dir)
                pbar.update(1)

        # ==============================================================================
        # PHASE 2: SEIDR (SCAN -> EXECUTE)
        # ==============================================================================
        print(f"\n[Phase 2] Scanning for completed networks...")

        # 1. Scan
        existing_nets = 0
        tasks_to_run = []
        for item in seidr_queue:
            d = item[2]
            if not self.force and (d / ".saturation.done").exists():
                existing_nets += 1
            else:
                tasks_to_run.append(item)

        print(f"          Found {existing_nets}/{len(seidr_queue)} completed networks.")

        # 2. Execute
        if tasks_to_run:
            threads_per_job = max(1, self.total_thread_budget // self.workers)

            with tqdm(total=len(seidr_queue), initial=existing_nets, desc="Seidr Inf", unit="net") as pbar:
                with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                    f_map = {
                        executor.submit(
                            run_seidr_batch,
                            cfg=self.config,
                            genes_file=i[0],
                            expression_file=i[1],
                            outdir=i[2],
                            preset=preset,
                            threads=threads_per_job
                        ): i for i in tasks_to_run
                    }

                    for f in concurrent.futures.as_completed(f_map):
                        try:
                            f.result()
                        except:
                            pass
                        pbar.update(1)
        else:
            print("          All Seidr tasks complete.")

        # ==============================================================================
        # PHASE 3: EGAD (SCAN -> EXECUTE)
        # ==============================================================================
        if self.mapman_file:
            print(f"\n[Phase 3] Scanning for EGAD evaluations...")
            try:
                r_script = get_egad_script_path()
            except:
                return

            # 1. Scan
            existing_egad = 0
            egad_tasks = []
            egad_results = []

            for item in egad_queue:
                # If cached, load result immediately
                if not self.force and item["out"].exists():
                    try:
                        val = pd.read_csv(item["out"], sep="\t")["AUC"].mean()
                        egad_results.append({'pct': item['pct'], 'iter': item['iter'], 'auc': val})
                        existing_egad += 1
                        continue
                    except:
                        pass  # If corrupt, re-run

                # Only add to queue if network exists
                if item["net"].exists():
                    egad_tasks.append(item)
                else:
                    # Network failed/missing, can't run EGAD
                    # We count it as "processed" (skipped) for the bar?
                    # Actually, usually better to ignore or log warning.
                    pass

            print(f"          Found {existing_egad}/{len(egad_queue)} existing evaluations.")

            # 2. Execute
            if egad_tasks:
                with tqdm(total=len(egad_queue), initial=existing_egad, desc="EGAD Eval", unit="task") as pbar:
                    with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                        f_map = {
                            executor.submit(
                                run_egad_task,
                                network_file=i["net"],
                                mapman_file=self.mapman_file,
                                out_file=i["out"],
                                script_path=r_script,
                                log_path=i["dir"] / "egad.log"
                            ): i for i in egad_tasks
                        }

                        for f in concurrent.futures.as_completed(f_map):
                            item = f_map[f]
                            res = f.result()
                            if res is not None:
                                egad_results.append({'pct': item['pct'], 'iter': item['iter'], 'auc': res})
                            pbar.update(1)
            else:
                print("          All EGAD tasks complete (or networks missing).")

            aggregate_and_plot_egad(egad_results, self.base_outdir)
        else:
            print("\n[Phase 3] Skipped (No MapMan).")