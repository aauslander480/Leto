# Getting started

[Repository overview](../../README.md) · [Workflow and outputs](workflow.md) ·
[Troubleshooting](troubleshooting.md)

## 1. Get LETO and prepare the environment

Clone the repository or download it with GitHub's **Code → Download ZIP** and
extract the full archive. Open the `LETO_v1_project` folder; this contains
`run_leto.py`, `configs`, and the user guides.

Use a licensed ArcGIS Pro installation on Windows for a full run. Select or
clone an ArcGIS Pro Python environment through Pro's environment manager and
use its Python Command Prompt or notebook kernel. A standard Python virtual
environment cannot supply ArcPy through LETO's pip dependencies.

From the selected environment, check imports:

```bat
python -c "import sys, arcpy, numpy, pandas, scipy; print(sys.executable); print(arcpy.GetInstallInfo()); print(numpy.__version__, pandas.__version__, scipy.__version__)"
```

The package metadata requires Python >=3.9 and declares only pandas as a runtime
dependency. ArcPy, NumPy, and SciPy are also needed by the spatial stages.
Reading the supplied legacy `.xls` species workbook may require the
`xlrd` Excel reader in that same environment; `.xlsx` normally uses
`openpyxl`. Manage missing packages in the selected Pro environment.

The repository does not pin a tested ArcGIS Pro version or provide a complete
locked GIS environment. Record your Pro and dependency versions with your run.

Running `run_leto.py` from the project folder does not require an editable
package install. Prepare the Pro environment above, then continue with your data.

## 2. Gather the inputs

Use [Input data and downloads](input_data.md) for official source links,
which products to download, and FVS template compatibility notes.

Large inputs are not included. Obtain the relevant products from their
publishers (US Forest Service for TreeMap, FIA, and FVS; the ownership raster's
publisher; USGS for NHD; and US Fish and Wildlife Service for NWI). The example
uses TreeMap2022 and a specific ownership classification: a different product
or release is not automatically equivalent.

| Input / configuration key | Expected content |
|---|---|
| `study_area` | ArcGIS-readable study boundary feature class |
| `treemap_raster` | TreeMap raster covering the study area |
| `treemap_lookup_raster` | Full lookup raster retaining the raster attribute table; defaults to `treemap_raster` |
| `fia_cond_files` | FIA COND CSVs with `PLT_CN`, `CONDID`, `COND_STATUS_CD`, `CONDPROP_UNADJ`, `STDAGE`, `FORTYPCD`, `BALIVE` |
| `fia_tree_files` | FIA TREE CSVs used for raster summaries and tree inventories; retain the original export columns |
| `species_crosswalk`, `species_crosswalk_sheet` | Workbook and sheet translating FIA species to the chosen FVS codes |
| `blank_fvs_database` | Compatible blank SQLite database with `FVS_StandInit`, `FVS_TreeInit`, and `FVS_PlotInit` tables |
| `ownership_source_raster` | Ownership raster matching the configured code dictionary when ownership is enabled |
| `riparian.nhd_flowline`, `nhd_area`, `nhd_waterbody` | NHD feature classes for enabled riparian processing, with classification codes matching the buffer rules |
| `riparian.nwi_wetland_inputs` | Wetland feature inputs if wetland classification is enabled |

The TreeMap lookup table must retain `VALUE` and `PLT_CN` (recognized plot
field alternatives are `PLTCN` and `PLOT_CN`). Available `STDAGE`,
`FORTYPCD`, `BALIVE`, `QMD`, and `TPA_LIVE`/`TPA` attributes
provide fallbacks when FIA-derived attributes are missing. A clipped raster
without its original lookup attributes is not a replacement for that table.

Include FIA files covering the **TreeMap donor plots**, which can extend beyond
the study area's state. Select source files based on donor coverage rather
than administrative boundaries alone. Preserve control numbers as text; opening and resaving
them through a spreadsheet can lose identifier precision.

The example crosswalk expects sheet `EasternSpeciesTranslator`, an FIA column
matching `FIA CODE`, `FIA_CODE`, or `FIACODE`, and an FVS column matching
`SN_Mapped_To` or `SN MAPPED TO`. Check the actual translations against your
chosen variant; the column name alone does not establish compatibility.

## 3. Create a regional configuration

From `LETO_v1_project`, copy `configs/maine_test.py` to
`configs/my_region.py`. This file is a complete regional example, not a
universal set of modeling assumptions. Use the [configuration guide](configuration.md)
to review the copy's `CONFIG` dictionary:

1. Replace every input path, `working_gdb`, and `runs_root_folder`.
   Use absolute paths and Python raw strings such as `r"D:\LETO_Data\study.gdb\Boundary"`.
2. Set `region_id`, `run_name`, FVS variant, state abbreviation, numeric state
   code, and inventory year for your region.
3. Review ownership codes, segmentation assumptions, riparian buffers,
   wetland inputs, species translation, donor filtering, and imputation.
4. Select `database.mode`: `grouped`, `combined`, or `both`.
5. Keep all nested sections from the example. The loader supplies some defaults
   but does not fill every setting consumed by the stages.

Configuration files execute Python when loaded. Use configurations you trust.
A configuration-only check, run from the project folder, is:

```bat
python -c "from leto.config import load_config; c = load_config('configs/my_region.py'); print(c.region_id, c.run_name)"
```

This checks loading and selected configuration rules; it does **not** check all
input paths, fields, GIS licenses, or database compatibility.

Choose a separate working geodatabase for each independent run. LETO creates
timestamped output folders, but spatial intermediate and final datasets use
fixed names inside `working_gdb` and may be overwritten. Copy any spatial
outputs you need to retain before reusing that geodatabase. Run regions
sequentially with separate workspaces.

Select small-stand handling for your objectives and document the choice.
Validate each regional setup using the [validation guide](validation.md),
including acreage impacts, donor coverage, and any proximity grouping.
Expected counts depend on your study area, inputs, and settings.

## 4. Run one region

In an ArcGIS Pro Python Command Prompt, change to the cloned
`LETO_v1_project` folder:

```bat
python run_leto.py configs\my_region.py
```

Always pass the intended configuration explicitly. With no argument the entry
point selects `configs/maine_test.py`, including its original local paths.

In an ArcGIS Pro notebook, use quoted absolute paths (replace the example):

```python
%run "D:/Leto/LETO_v1_project/run_leto.py" "D:/Leto/LETO_v1_project/configs/my_region.py"
```

The pipeline runs all five stages. There is no stage-selection or resume flag
in this entry point. Start with a small study boundary to assess data coverage
and resource needs; the repository does not specify minimum RAM or run time.

## 5. Review the result

The console prints the new run folder. Check `Run_Summary.json`,
`Run_Summary.txt`, `Stage_Timings.csv`, and
`FVS_Database_Group_Summary.csv`, then inspect the QA files and final feature
classes described in [Workflow and outputs](workflow.md).

A completed run is not sufficient scientific validation. Review acreage,
donor normalization, species translation losses, imputation extent, and
stands without usable trees before using the databases. FVS execution is a
separate step in your FVS workflow, using the matching variant and validated
sampling assumptions.
