import os
import sys
import subprocess
import pandas as pd
from pathlib import Path
from typing import Optional


def get_egad_script_path() -> Path:
    """
    Locates the egad.R script.
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
        mapman_file: Path,
        out_file: Path,
        script_path: Path,
        log_path: Path
) -> Optional[float]:
    """
    Wraps the R script execution.
    """
    cmd = [
        "Rscript", str(script_path),
        "--network", str(network_file),
        "--mapman", str(mapman_file),
        "--output", str(out_file)
    ]

    with open(log_path, "a") as log:
        log.write(f"\n{'=' * 40}\n[EXEC] {' '.join(cmd)}\n")

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            text=True
        )
        with open(log_path, "a") as log:
            log.write(result.stderr)
            log.write("[SUCCESS]\n")

    except subprocess.CalledProcessError as e:
        with open(log_path, "a") as log:
            log.write(f"[ERROR] Output:\n{e.stderr}\n")
            log.write(f"[FAILED] Exit Code: {e.returncode}\n")
        return None

    if out_file.exists():
        try:
            df = pd.read_csv(out_file, sep="\t")
            if df.empty: return None
            return df["AUC"].mean()
        except Exception:
            return None
    return None