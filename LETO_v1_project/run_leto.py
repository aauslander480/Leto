r"""Run one LETO region from a Python configuration file.

ArcGIS Pro notebook example:
    %run r"C:\path\to\run_leto.py" r"C:\path\to\maine_test.py"

Command prompt example using the ArcGIS Pro Python environment:
    python run_leto.py configs\maine_test.py
"""
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_FOLDER = Path(__file__).resolve().parent
if str(PROJECT_FOLDER) not in sys.path:
    sys.path.insert(0, str(PROJECT_FOLDER))

from leto.config import load_config
from leto.pipeline import run_pipeline


def main(config_path: str | None = None):
    if config_path is None:
        if len(sys.argv) < 2:
            config_path = str(PROJECT_FOLDER / "configs" / "maine_test.py")
        else:
            config_path = sys.argv[1]

    config = load_config(config_path)
    return run_pipeline(config)


if __name__ == "__main__":
    main()
