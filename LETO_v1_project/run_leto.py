r"""Run one LETO region from a Python configuration file.

ArcGIS Pro notebook example:
    %run r"D:\LETO\run_leto.py" r"D:\LETO\configs\maine_test.py"

Command prompt example using the ArcGIS Pro Python environment:
    python run_leto.py configs\maine_test.py
"""
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_FOLDER = Path(__file__).resolve().parent
# ArcGIS notebooks keep imported modules alive between %run calls. A prior LETO
# checkout can therefore remain in sys.modules even when a newer run_leto.py is
# executed from another folder. Put this checkout first and evict every cached
# LETO module before importing the pipeline.
project_folder_text = str(PROJECT_FOLDER)
sys.path[:] = [
    entry for entry in sys.path
    if str(Path(entry or ".").resolve()).casefold()
    != project_folder_text.casefold()
]
sys.path.insert(0, project_folder_text)
for module_name in list(sys.modules):
    if module_name == "leto" or module_name.startswith("leto."):
        del sys.modules[module_name]

import leto
from leto.config import load_config
from leto.pipeline import run_pipeline


def main(config_path: str | None = None):
    loaded_package = Path(leto.__file__).resolve().parent
    expected_package = PROJECT_FOLDER / "leto"
    if loaded_package != expected_package:
        raise RuntimeError(
            "LETO imported from a different checkout.\n"
            f"Expected: {expected_package}\n"
            f"Loaded: {loaded_package}"
        )
    print(f"LETO source folder: {PROJECT_FOLDER}", flush=True)

    if config_path is None:
        if len(sys.argv) < 2:
            config_path = str(PROJECT_FOLDER / "configs" / "maine_test.py")
        else:
            config_path = sys.argv[1]

    config = load_config(config_path)
    return run_pipeline(config)


if __name__ == "__main__":
    main()

