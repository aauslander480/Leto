"""Run LETO using QGIS's Python: python3 run_leto.py configs/maine_sapelo2.py."""
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_FOLDER = Path(__file__).resolve().parent
if str(PROJECT_FOLDER) not in sys.path:
    sys.path.insert(0, str(PROJECT_FOLDER))

from leto.config import load_config


def main(config_path: str | None = None):
    if config_path is None:
        if len(sys.argv) < 2:
            raise SystemExit('Provide a regional configuration or --check-runtime.')
        else:
            config_path = sys.argv[1]

    from leto import gis_backend
    gis_backend.initialize()
    if config_path == '--check-runtime':
        import numpy, pandas, scipy, xlrd
        print(f'QGIS {gis_backend.q.Qgis.QGIS_VERSION}; {gis_backend.gdal.VersionInfo("--version")}')
        print(f'NumPy {numpy.__version__}; pandas {pandas.__version__}; SciPy {scipy.__version__}')
        return
    from leto.pipeline import run_pipeline
    config = load_config(config_path)
    return run_pipeline(config)


if __name__ == "__main__":
    main()

