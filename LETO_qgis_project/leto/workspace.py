"""Run-folder and geodatabase management."""
from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path

from . import gis_backend as gis

from .stage5_database import safe_filename


def ensure_working_gdb(path: str) -> None:
    if str(path).lower().endswith('.gdb'):
        raise ValueError('For QGIS, working_gdb must be a workspace directory, not a .gdb.')
    if gis.Exists(path):
        return
    parent = os.path.dirname(path)
    os.makedirs(parent, exist_ok=True)
    gis.management.create_workspace(parent, os.path.basename(path))


def create_run_paths(runs_root: str, run_name: str) -> dict[str, str]:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    root = Path(runs_root)
    root.mkdir(parents=True, exist_ok=True)
    base = root / f"{timestamp}_{safe_filename(run_name)}"
    run = base
    counter = 2
    while run.exists():
        run = Path(f"{base}_{counter}")
        counter += 1

    paths = {
        "run": run,
        "intermediate": run / "02_Intermediate",
        "rasters": run / "02_Intermediate" / "Segmentation_Rasters",
        "segmentation": run / "03_Segmentation",
        "crosswalk": run / "04_Crosswalk",
        "fvs_inputs": run / "05_FVS_Inputs",
        "databases": run / "06_FVS_Databases",
        "qa": run / "07_QA",
        "logs": run / "08_Logs",
    }
    for path in paths.values():
        Path(path).mkdir(parents=True, exist_ok=True)
    return {key: str(value) for key, value in paths.items()}

