# Getting started

[Overview](../../README.md) · [Input data](input_data.md) ·
[Configuration](configuration.md) · [Troubleshooting](troubleshooting.md)

## 1. Prepare the environment

Clone or download the repository and open `LETO_v1_project`, which contains
`run_leto.py` and `configs`.

Use a licensed ArcGIS Pro installation on Windows. Select or clone a Pro
Python environment through Pro's environment manager, then use its Python
Command Prompt or notebook kernel. Check the interpreter and required imports:

```bat
python -c "import sys, arcpy, numpy, pandas, scipy; print(sys.executable); print(arcpy.GetInstallInfo()); print(numpy.__version__, pandas.__version__, scipy.__version__)"
```

LETO requires Python >=3.9. Its package metadata declares pandas; spatial runs
also require ArcPy, NumPy, and SciPy. Reading `.xls` workbooks may require
`xlrd`; `.xlsx` normally uses `openpyxl`. Manage missing packages in
the selected Pro environment.

An ordinary Python virtual environment cannot supply ArcPy through LETO's pip
dependencies. Running `run_leto.py` from this folder does not require an
editable package install.

The repository has no locked GIS environment or specified tested Pro version.
Record your environment versions with the run.

## 2. Prepare input data

Use [Input data and downloads](input_data.md) to obtain the datasets, species
workbook, and FVS template. Check the template's schema compatibility before
running. That guide also covers donor-state coverage and optional layers.

Required COND fields are `PLT_CN`, `CONDID`, `COND_STATUS_CD`,
`CONDPROP_UNADJ`, `STDAGE`, `FORTYPCD`, and `BALIVE`. Retain
the original TREE export columns.

The TreeMap lookup requires `VALUE` and `PLT_CN` (or `PLTCN` /
`PLOT_CN`). Available `STDAGE`, `FORTYPCD`, `BALIVE`, `QMD`,
and `TPA_LIVE` / `TPA` attributes provide fallbacks for missing FIA
attributes. Keep the full lookup raster if clipping removes these attributes.

## 3. Configure your study area

From `LETO_v1_project`, in PowerShell:

```powershell
Copy-Item configs\maine_test.py configs\my_region.py
```

Edit the copy's `CONFIG` dictionary using the
[configuration guide](configuration.md). Replace all paths, set regional/FVS
metadata, and review the analysis settings. Keep all nested sections; the
loader supplies only some defaults.

Use a dedicated working geodatabase. Reusing it can overwrite spatial
datasets even though each run creates a new timestamped results folder.
Archive spatial outputs before reuse and run regions sequentially.

## 4. Run LETO

In an ArcGIS Pro Python Command Prompt, from `LETO_v1_project`:

```bat
python run_leto.py configs\my_region.py
```

In an ArcGIS Pro notebook, replace these paths with your own:

```python
%run "D:/Leto/LETO_v1_project/run_leto.py" "D:/Leto/LETO_v1_project/configs/my_region.py"
```

Always pass your configuration. With no argument, the entry point selects
`configs/maine_test.py`, including placeholder paths.

The entry point runs all five stages and has no stage-selection or resume
flag. Start with a small boundary to assess coverage and resource needs;
minimum RAM and run time are not specified.

## 5. Review results

Open the run folder printed in the console. Start with `Run_Summary.txt`,
`Run_Summary.json`, `Stage_Timings.csv`, and
`FVS_Database_Group_Summary.csv`.

Use [Workflow and outputs](workflow.md) to locate spatial layers and QA
reports, then follow [Validation](validation.md). Review acreage, donor
weights, species losses, imputation, and units without usable trees before
using the databases in FVS.
