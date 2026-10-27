# Workflow and outputs

[Getting started](getting_started.md) · [Validation](validation.md)

## Five stages

| Stage | Processing | Main results |
|---|---|---|
| 1 | Derive segmentation attributes from FIA; fall back to TreeMap attributes and then nearest assigned cells | Aligned age, forest type, basal area, QMD, and TPA rasters; imputation QA |
| 2 | Segment the landscape, apply configured small-stand and management rules, finalize spatial identities | Component polygons and dissolved management units |
| 3 | Crosswalk TreeMap donors to parent segments, then aggregate by contributing component acreage | Management-unit stand list and normalized donor weights |
| 4 | Filter and renormalize donors, translate species, weight FIA trees, optionally impute missing tree lists | FVS stand/tree inventories, missing-species and imputation diagnostics |
| 5 | Copy the blank FVS database and populate inventory tables | Grouped and/or combined SQLite databases |

## Spatial identifiers

- `SEG_ID` identifies original biological segments in the spatial outputs;
  `ParentSEG_ID` records parent provenance.
- `BIO_MU_ID` records biological grouping before routing compatibility splits.
- `MU_ID` identifies the final management unit and FVS routing unit.
- `MU_PART` identifies component pieces deterministically.
- `STAND_ID` links FVS stand and tree tables.

`ManagementSegments_Final` retains singlepart components and provenance.
`ManagementUnits_Final` dissolves them to one potentially multipart feature
per final `MU_ID`. Their feature counts can differ.

Some inventory/crosswalk code retains the legacy column name `SEG_ID` for the
final unit key. Interpret it in that table's context; do not join an inventory
table to original parent segments on that name alone.

## Output locations

Each run creates a timestamped folder beneath `runs_root_folder`:

```text
<timestamp>_<run_name>/
  Run_Settings.json
  Run_Summary.json
  Run_Summary.txt
  Stage_Timings.csv
  FVS_Database_Group_Summary.csv
  02_Intermediate/
    Segmentation_Rasters/
  03_Segmentation/
  04_Crosswalk/
  05_FVS_Inputs/
  06_FVS_Databases/
  07_QA/
  08_Logs/
```

Capture console or notebook output for tracebacks; `08_Logs` may contain no
persistent log. Final feature classes are stored in the working geodatabase.
Timestamped results folders do not preserve previous versions of those layers.

| Location | Files / datasets to inspect |
|---|---|
| Working geodatabase | `ManagementSegments_Final`, `ManagementUnits_Final`, spatial intermediates |
| `04_Crosswalk` | `Segment_Stand_List.csv`, `Segment_PLT_CN_Weights.csv`, `Segment_PLT_CN_Weights_Filtered.csv` |
| `05_FVS_Inputs` | `FVS_StandInit.csv`, `FVS_TreeInit.csv` |
| `06_FVS_Databases` | `FVS_Inventory_<group>.db` |
| `07_QA` | Missing-species, missing-tree, imputation, management-unit, and proximity diagnostics |

Stage 1 writes `QA_STAGE1_Raster_Imputation_Summary.csv`,
`QA_STAGE1_Landscape_Imputation_Summary.csv`, and optional provenance rasters
with its raster outputs. Stage 4 writes `QA_Imputation_Diagnostics.csv`,
`QA_Imputation_Summary.csv`, `QA_Missing_FVS_Species.csv`,
`QA_Missing_FVS_Species_By_Stand.csv`, and
`Segment_FVS_Stands_No_Live_Trees.csv` in `07_QA`.

`Run_Settings.json` records the resolved configuration and input paths, not
copies of the source datasets. Archive the source versions and environment
information separately for reproducibility.

## Inventory assumptions

Donor filtering uses `inventory.minimum_plot_weight`, retains at least the
largest donor, and renormalizes weights. Tree height falls back from `HT` to
`ACTUALHT` when available.

Units without usable trees may be excluded from the runnable inventory.
Review `Segment_FVS_Stands_No_Live_Trees.csv` even when processing completes.

The example uses weighted trees per acre for `TREE_COUNT`, with basal area
factor -1, plot size 1, breakpoint DBH 0, one plot, and zero nonstockable plots.
Keep sampling design consistent with your FVS template.

In `grouped` mode, group stand totals should match unique final stands.
In `both` mode, the combined database repeats the grouped inventory.

Follow [Validation](validation.md) to review coverage, weights, species,
imputation, and database integrity before using the results.
