# Troubleshooting

[Setup guide](getting_started.md) · [Workflow and outputs](workflow.md)

| Symptom | What to check |
|---|---|
| `No module named arcpy` | Use the ArcGIS Pro Python environment or notebook kernel. An ordinary venv and `pip install -e .` do not install ArcPy. |
| Missing NumPy, SciPy, pandas, or Excel reader | Check the interpreter printed by `sys.executable`; manage the dependency in that same Pro environment. |
| Configuration not found / missing `CONFIG` | Pass the configuration path explicitly. The file must define a dictionary named `CONFIG`. |
| Missing nested key / `KeyError` | Start from a complete copy of `configs/maine_test.py`; loading does not populate every nested setting. |
| An example drive path appears | Replace every path in your regional configuration, including nested ownership, NHD, and optional wetland paths. |
| TreeMap lookup lacks `PLT_CN` or VALUE codes | Point `treemap_lookup_raster` at the original raster with its attribute table. Preserve the mapping from raster VALUE to donor plot. |
| Many donors absent from FIA files | Confirm the donor states and source release, not just the study state. Preserve control-number text precision. |
| Species workbook fails / translations missing | Check file format, Excel reader, sheet and column names, and the translation codes for the chosen FVS variant. Review missing-species QA. |
| Raster alignment failure | Review extent, dimensions, coordinate system, cell size, and snap alignment. Categorical rasters require nearest-neighbor resampling. |
| Geodatabase locks or overwritten spatial outputs | Close layers holding locks, use a dedicated working geodatabase, and avoid concurrent runs against it. Archive prior spatial results. |
| No runnable stands for nearest imputation | Check FIA coverage, live-tree records, retained donors, and species translation. Nearest imputation needs at least one valid donor stand. |
| SQLite missing table / column errors | Supply a compatible blank FVS SQLite template with the expected inventory schema; a new empty SQLite file is insufficient. |
| Unexpected counts or acreage | Compare against your study footprint and a previously accepted run with the same input versions and settings. Review segmentation, donor coverage, and imputation QA using the [validation guide](validation.md). |

A failed run may leave partial output files and intermediate geodatabase
datasets. Root-level run summaries are written at the end, so their absence
can indicate an incomplete run. Retain the traceback and partial QA for
diagnosis; the entry point does not resume from a selected stage.

For an issue report, include the command, failing stage, traceback, environment
versions, expected result, and relevant configuration or QA rows. Remove
private paths and source data.
