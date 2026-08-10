# LETO

LETO is a spatial forest landscape modeling workflow that integrates TreeMap,
Forest Inventory and Analysis (FIA) data, ArcGIS Pro, and the Forest Vegetation
Simulator (FVS).

The actively developed implementation lives in `LETO_v1_project/`. The older
workflow scripts remain in `scripts/` for reference while the v1 package is used
for new development.

## Current development baseline

The v1 package reorganizes the verified Maine LETO pipeline into a maintainable,
one-region software structure without intentionally changing the scientific
algorithms. Before adding new algorithms, the package should reproduce these
verified baseline metrics:

- 429,977 valid raster cells
- 12,168 parent segments
- 19,662 final management-unit polygon records
- 19,244 direct FVS stands
- 418 nearest-imputed stands
- 6 grouped FVS databases

See `LETO_v1_project/README.md` and `LETO_v1_project/docs/` for workflow and
baseline-validation details.

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
