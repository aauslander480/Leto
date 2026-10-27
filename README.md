# LETO

LETO creates forest management units and weighted Forest Vegetation Simulator
(FVS) inventories from TreeMap and Forest Inventory and Analysis (FIA) data.
It produces GIS layers, stand and tree CSVs, SQLite databases, and quality
assurance (QA) reports. FVS simulations are run separately.

## Requirements

LETO is intended for analysts familiar with GIS data, Python configuration
files, and forest inventory assumptions.

- Windows and a licensed ArcGIS Pro Python environment.
- A study boundary and TreeMap raster with donor-plot lookup attributes.
- FIA TREE and COND exports covering the landscape's donor plots.
- A species crosswalk workbook and compatible blank FVS SQLite database.
- Ownership, hydrography, and wetland data for enabled options.

See [Input data and downloads](LETO_v1_project/docs/input_data.md) for sources
and compatibility requirements. Input datasets are not included.

## Quick start

Clone the repository, or use GitHub's **Code → Download ZIP** and extract it:

```powershell
git clone https://github.com/aauslander480/Leto.git
cd Leto
git switch main
cd LETO_v1_project
Copy-Item configs\maine_test.py configs\my_region.py
```

Follow [Getting started](LETO_v1_project/docs/getting_started.md) to prepare
your environment and data. Edit `configs/my_region.py` with your paths,
regional metadata, and modeling settings.

From `LETO_v1_project` in an ArcGIS Pro Python Command Prompt:

```bat
python run_leto.py configs\my_region.py
```

Pass your configuration explicitly. Start with a small study area and review
QA before processing the full landscape. Notebook instructions are in the
setup guide.

## Results

The console prints a timestamped results folder containing configuration,
summaries, inventories, databases, and QA reports. Final spatial layers are
stored in your working geodatabase.

Open `ManagementUnits_Final` in ArcGIS Pro and review `Run_Summary.txt`
and the QA reports. Check acreage, donor coverage, species losses, and
imputation before using the inventories.

## Documentation

| Guide | Contents |
|---|---|
| [Getting started](LETO_v1_project/docs/getting_started.md) | Environment setup and first run |
| [Input data](LETO_v1_project/docs/input_data.md) | Downloads, formats, and compatibility |
| [Configuration](LETO_v1_project/docs/configuration.md) | Paths and analysis settings |
| [Workflow and outputs](LETO_v1_project/docs/workflow.md) | Processing stages, files, and identifiers |
| [Validation](LETO_v1_project/docs/validation.md) | Checks before using results |
| [Troubleshooting](LETO_v1_project/docs/troubleshooting.md) | Setup and processing errors |
| [Advanced behavior](LETO_v1_project/docs/advanced_behavior.md) | Grouping and imputation details |

## Help and contributions

For unresolved problems, [open an issue](https://github.com/aauslander480/Leto/issues)
with the command, expected result, traceback, environment versions, and relevant
QA summaries. Remove private paths and data.

See [Contributing](CONTRIBUTING.md) for development setup and tests.
