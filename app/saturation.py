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
            return pd.read_csv(self.vst_path, sep="\t", index_col=0)

    def _write_seidr_files(self, df: pd.DataFrame, out_dir: Path):
        df = df.fillna(0)
        g_file, e_file = out_dir / "genes.txt", out_dir / "expression.tsv"
        with open(g_file, "w") as f:
            f.write("\n".join(df.columns.astype(str)))
        df.to_csv(e_file, sep="\t", header=False, index=False, float_format='%.6f')
        return g_file, e_file

    def run(self):
        total_samples = len(self.dataset)
        steps_pct = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

        all_projects_meta = []
        for bp in self.dataset.bioprojects:
            s_ids = [s.id for s in bp.samples]
            if s_ids: all_projects_meta.append((bp, len(s_ids), s_ids))

        preset = getattr(self.config, "seidr_preset", "FAST")
        work_definitions = []
        for pct in steps_pct:
            iters = 1 if pct == 100 else self.iterations
            for i in range(1, iters + 1):
                work_definitions.append((pct, i))

        print(f"\n[Phase 1] Scanning for existing batch inputs...")
        existing_batches = 0
        if not self.force:
            for pct, i in work_definitions:
                d = self.base_outdir / f"{pct}%_batch" / f"iter{i}"
                if pct == 100:
                    # For 100%, we only check if the network exists
                    if (d / ".saturation.done").exists() or (d / "network_saturation_edges.tsv").exists():
                        existing_batches += 1
                elif (d / ".saturation.done").exists() or (
                        (d / "genes.txt").exists() and (d / "expression.tsv").exists()):
                    existing_batches += 1

        seidr_queue, egad_queue = [], []

        with tqdm(total=len(work_definitions), initial=existing_batches, desc="Batch Gen", unit="batch") as pbar:
            available_pool = list(all_projects_meta)
            last_pct = None

            for pct, i in work_definitions:
                if pct != last_pct:
                    available_pool, last_pct = list(all_projects_meta), pct

                iter_dir = self.base_outdir / f"{pct}%_batch" / f"iter{i}"
                iter_dir.mkdir(parents=True, exist_ok=True)

                # --- QUEUE LOGIC ---
                if pct == 100:
                    # POINT DIRECTLY TO DESEQ2 FILES
                    g_path = self.deseq2_dir / "genes.txt"
                    e_path = self.deseq2_dir / "expression.tsv"
                else:
                    g_path, e_path = iter_dir / "genes.txt", iter_dir / "expression.tsv"

                seidr_queue.append((g_path, e_path, iter_dir, pct, i))
                egad_queue.append({"net": iter_dir / "network_saturation_edges.tsv", "out": iter_dir / "egad_auroc.tsv",
                                   "dir": iter_dir, "pct": pct, "iter": i})

                # Skip generation if files exist OR if it's the 100% batch (no gen needed)
                if pct == 100 or (not self.force and (g_path.exists() and e_path.exists())):
                    pbar.update(1)
                    continue

                # --- STANDARD SAMPLING LOGIC (Only for < 100%) ---
                selected_samples = []
                n_target = int(total_samples * (pct / 100.0))
                current_n = 0
                used_in_this_run = set()
                while current_n < n_target:
                    rem = n_target - current_n
                    cands = [(bp, sz, sids) for (bp, sz, sids) in available_pool if
                             sz <= rem and bp not in used_in_this_run]
                    if cands:
                        w = np.array([1.0 / (sz ** 2) for (_, sz, _) in cands])
                        chosen = cands[np.random.choice(len(cands), p=w / w.sum())]
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

                matrix = self._get_master_df()
                final_s = [s for s in selected_samples if not matrix.empty and s in matrix.index]
                if final_s:
                    sub = matrix.loc[final_s].copy().fillna(0)
                    sub = sub.loc[:, sub.var(axis=0) > 0.1]
                    if not sub.empty: self._write_seidr_files(sub, iter_dir)
                pbar.update(1)

        print(f"\n[Phase 2] Executing Seidr inference...")
        tasks = [it for it in seidr_queue if self.force or not (it[2] / ".saturation.done").exists()]
        if tasks:
            tpj = max(1, self.total_thread_budget // self.workers)
            with tqdm(total=len(seidr_queue), initial=len(seidr_queue) - len(tasks), desc="Seidr Inf",
                      unit="net") as pbar:
                with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                    f_map = {executor.submit(run_seidr_batch, self.config, i[0], i[1], i[2], preset, tpj): i for i in
                             tasks}
                    for f in concurrent.futures.as_completed(f_map): pbar.update(1)

        if self.mapman_file:
            print(f"\n[Phase 3] Running EGAD evaluations...")
            r_script = get_egad_script_path()
            egad_results, egad_tasks = [], []
            for item in egad_queue:
                if not self.force and item["out"].exists():
                    try:
                        val = pd.read_csv(item["out"], sep="\t")["AUC"].mean()
                        egad_results.append({'pct': item['pct'], 'iter': item['iter'], 'auc': val})
                        continue
                    except:
                        pass
                if item["net"].exists(): egad_tasks.append(item)

            if egad_tasks:
                with tqdm(total=len(egad_queue), initial=len(egad_queue) - len(egad_tasks), desc="EGAD Eval",
                          unit="task") as pbar:
                    with concurrent.futures.ThreadPoolExecutor(max_workers=self.workers) as executor:
                        f_map = {executor.submit(run_egad_task, i["net"], self.mapman_file, i["out"], r_script,
                                                 i["dir"] / "egad.log"): i for i in egad_tasks}
                        for f in concurrent.futures.as_completed(f_map):
                            item, res = f_map[f], f.result()
                            if res is not None: egad_results.append(
                                {'pct': item['pct'], 'iter': item['iter'], 'auc': res})
                            pbar.update(1)
            aggregate_and_plot_egad(egad_results, self.base_outdir)