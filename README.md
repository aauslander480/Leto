# LETO

LETO prepares forest landscape units and weighted Forest Vegetation Simulator
(FVS) inventories from TreeMap and Forest Inventory and Analysis (FIA) data.
It uses ArcGIS Pro/ArcPy for spatial processing and writes stand/tree CSVs,
QA reports, and FVS-ready SQLite databases.

## What it does

1. Builds aligned rasters for forest type, age, basal area, quadratic mean
   diameter, and trees per acre.
2. Segments the landscape and creates management units with configured
   ownership, riparian, and wetland rules.
3. Calculates management-unit-to-FIA-plot donor weights.
4. Builds weighted tree lists, with optional nearest-stand imputation.
5. Writes databases grouped by management/ownership routing, combined, or both.

LETO prepares FVS inputs; the active pipeline does **not** run FVS simulations
or join simulation results back to maps.

## Before you start

- Full spatial runs require Windows and a licensed ArcGIS Pro Python
  environment with ArcPy, NumPy, pandas, and SciPy.
- Input datasets, an FIA-to-FVS species crosswalk workbook, and a compatible
  blank FVS SQLite database must be supplied separately.
- Configure LETO for your study area, source datasets, and FVS variant.
  The included regional example contains local paths that must be replaced.
- This research workflow requires user-supplied data and regional validation;
  no sample-data demo is bundled.

## Get started

Clone the repository using a GitHub account that has access:

```powershell
git clone https://github.com/aauslander480/Leto.git
cd Leto
git switch main
cd LETO_v1_project
Copy-Item configs\maine_test.py configs\my_region.py
```

Edit `configs/my_region.py`, then use an **ArcGIS Pro Python Command Prompt**
with your prepared environment:

```bat
python run_leto.py configs\my_region.py
```

Follow the [setup and input-data guide](LETO_v1_project/docs/getting_started.md)
before this command. It covers environment checks, required inputs, notebook
execution, configuration, and the first-run review.

## Documentation

| Guide | Contents |
|---|---|
| [Getting started](LETO_v1_project/docs/getting_started.md) | Environment, inputs, configuration, and running |
| [Workflow and outputs](LETO_v1_project/docs/workflow.md) | Five stages, identifiers, output locations, and QA |
| [Validation](LETO_v1_project/docs/validation.md) | Regional acceptance checks and reproducibility |
| [Advanced behavior](LETO_v1_project/README.md) | Proximity grouping, spatial contract, and imputation diagnostics |
| [Troubleshooting](LETO_v1_project/docs/troubleshooting.md) | Common setup, data, and run failures |

## Repository layout

```text
LETO_v1_project/          Active implementation
  run_leto.py            One-region pipeline entry point
  configs/               Regional configuration example
  leto/                  Configuration, stages, orchestration, and QA
  tests/                 Unit tests that do not require ArcGIS
  docs/                  User guides and validation notes
scripts/                 Legacy/reference scripts; separate from the v1 entry point
AGENTS.md                Scientific integrity and development rules
```

## Development checks without ArcGIS

The package declares Python 3.9 or newer. These commands install the lightweight
package and test dependencies; they do not install ArcPy or the full GIS runtime.

Windows PowerShell, from the repository root:

```powershell
cd LETO_v1_project
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[dev]"
.\.venv\Scripts\python.exe -m pytest
```

macOS/Linux, for the unit tests:

```bash
cd LETO_v1_project
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[dev]'
.venv/bin/python -m pytest
```

Passing unit tests does not validate an end-to-end spatial run. Read
[AGENTS.md](AGENTS.md) before contributing; preserve reproducible regional configurations and document
any intentional changes to scientific assumptions.

When reporting a problem, include the failing stage, traceback, environment
versions, relevant configuration with private paths removed, and available QA
summaries. Do not attach large source datasets.
