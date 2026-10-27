# LETO

Turn TreeMap and Forest Inventory and Analysis (FIA) data into mapped forest
management units and weighted Forest Vegetation Simulator (FVS) inventories
for your study area.

LETO brings landscape segmentation, donor-plot weighting, inventory preparation,
and quality checks into one configurable workflow. Your outputs include GIS
feature classes, stand and tree CSV files, and FVS-ready SQLite databases.

## Is LETO a fit for your project?

LETO is intended for forest analysts and researchers who work with ArcGIS Pro,
TreeMap/FIA data, and FVS. You should be comfortable preparing GIS datasets,
editing a Python configuration file, and reviewing inventory assumptions.

You will need:

- Windows with a licensed ArcGIS Pro Python environment.
- A study boundary and TreeMap raster with its donor-plot lookup attributes.
- FIA TREE and COND exports covering the donor plots in your landscape.
- A species translation workbook and compatible blank FVS SQLite database.
- Ownership, hydrography, and wetland data for the options you enable.

Input data are supplied separately; see the
[download directory](LETO_v1_project/docs/input_data.md) for sources and file choices. LETO prepares inventories; running FVS
simulations is a separate step. There is no graphical setup wizard or bundled
sample-data run.

## Start with your own study area

1. **Get the source.** Clone the repository or use GitHub's **Code → Download
   ZIP** and extract the complete archive.
2. **Prepare your environment and inputs.** Follow the
   [getting-started guide](LETO_v1_project/docs/getting_started.md).
3. **Create your configuration.** Copy the supplied regional example, replace
   paths, and review the [configuration choices](LETO_v1_project/docs/configuration.md).
4. **Run a small study area first.** From an ArcGIS Pro Python Command Prompt,
   run the command below.
5. **Review your results.** Check spatial coverage, species translation, donor
   weights, and imputation before using the FVS databases.

To clone with Git:

```powershell
git clone https://github.com/aauslander480/Leto.git
cd Leto
git switch main
```

From the downloaded repository's `LETO_v1_project` folder, in PowerShell:

```powershell
Copy-Item configs\maine_test.py configs\my_region.py
```

The source file is a regional example. Replace its paths and review all
modeling settings for your project before running:

```bat
python run_leto.py configs\my_region.py
```

You can also [run from an ArcGIS Pro notebook](LETO_v1_project/docs/getting_started.md).
Always pass your configuration explicitly.

## Find what you need

| Your task | Guide |
|---|---|
| Find and download input datasets | [Input data and downloads](LETO_v1_project/docs/input_data.md) |
| Install prerequisites and prepare input data | [Getting started](LETO_v1_project/docs/getting_started.md) |
| Choose settings for your landscape | [Configuration](LETO_v1_project/docs/configuration.md) |
| Locate maps, inventories, databases, and QA reports | [Workflow and outputs](LETO_v1_project/docs/workflow.md) |
| Decide whether results are suitable for your analysis | [Validation](LETO_v1_project/docs/validation.md) |
| Resolve an error or unexpected result | [Troubleshooting](LETO_v1_project/docs/troubleshooting.md) |
| Understand grouping and imputation in detail | [Advanced behavior](LETO_v1_project/docs/advanced_behavior.md) |

## What to expect after a run

The console prints a timestamped results folder containing your configuration,
run summaries, inventory CSVs, databases, and QA reports. Final spatial layers
are written to your configured working geodatabase. Open
`ManagementUnits_Final` in ArcGIS Pro to inspect the management units and
review `Run_Summary.txt` alongside the QA reports.

Each region needs its own data and settings. A successful run still requires
review of acreage, donor coverage, species losses, and imputed inventories.

## Questions and problems

Check the [troubleshooting guide](LETO_v1_project/docs/troubleshooting.md) first.
If you still need help, [open an issue](https://github.com/aauslander480/Leto/issues)
with what you tried, what you expected, the traceback, environment versions,
and relevant QA summaries. Remove private paths and data before posting.

For code changes and local test setup, see [Contributing](CONTRIBUTING.md).
