# LETO

LETO is a spatial forest landscape modeling workflow that integrates TreeMap,
Forest Inventory and Analysis (FIA) data, ArcGIS Pro, and the Forest Vegetation
Simulator (FVS).

The actively developed implementation lives in `LETO_v1_project/`. The older
workflow scripts remain in `scripts/` for reference while the v1 package is used
for new development.

An independent, experimental PyQGIS/GDAL port is in `LETO_qgis_project/`.
It includes a Sapelo2 configuration and Slurm template. The original ArcPy
implementation is unchanged. Read the port README for runtime requirements,
scientific-validation status, and unsupported geodesic-near configurations.


## Repository layout

```text
LETO_v1_project/          Active Python package and tests
LETO_v1_project/docs/     v1 workflow and validation notes
LETO_v1_project/configs/  Example run configuration
scripts/                  Legacy/reference workflow scripts
workflow/                 Original workflow notes
```

## Developer quick start

From the repository root:

```bash
cd LETO_v1_project
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest
```

ArcGIS/ArcPy and FVS are required for end-to-end spatial workflow runs, but the
unit tests are designed to run without ArcGIS.

## Running the v1 pipeline in ArcGIS Pro

Edit paths and region-specific settings in `LETO_v1_project/configs/maine_test.py`,
then run the package entry script from an ArcGIS Pro Python environment:

```python
%run D:\path\to\LETO_v1_project\run_leto.py D:\path\to\LETO_v1_project\configs\maine_test.py
```

Large input datasets are not distributed with this repository. Users must
download TreeMap, FIA, ownership, and FVS reference data from their official
sources.

