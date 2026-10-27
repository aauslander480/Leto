# Configure LETO for your project

[Getting started](getting_started.md) · [Workflow and outputs](workflow.md)

You configure LETO by editing a Python file containing a dictionary named
`CONFIG`. You do not need to edit the processing modules for a normal run.

Copy the complete `configs/maine_test.py` example to `configs/my_region.py`.
Its filename identifies the original example; its paths, state settings, and
management rules must be reviewed for your study area. Keep all nested sections,
because the loader does not supply every field consumed by the pipeline.

## Set your data and workspace

| Settings | What you choose |
|---|---|
| `region_id`, `run_name` | Labels identifying your region and analysis |
| `study_area` | Boundary defining the landscape you want to process |
| `treemap_raster`, `treemap_lookup_raster` | Spatial source and original raster attribute lookup; retain donor identifiers |
| `fia_tree_files`, `fia_cond_files` | Source files covering your TreeMap donor plots |
| `working_gdb` | Dedicated file geodatabase for spatial processing and final layers |
| `runs_root_folder` | Parent folder for timestamped results |
| `species_crosswalk`, `species_crosswalk_sheet` | Species translation workbook and sheet |
| `blank_fvs_database` | Compatible blank FVS SQLite template |

Use absolute paths, for example `r"D:\LETO_Data\study.gdb\Boundary"`.
Replace nested paths too, including ownership, NHD, and wetland sources.
Avoid concurrent runs writing to the same working geodatabase.

## Review the analysis choices

| Section | Decisions for your study |
|---|---|
| FVS metadata | Set variant, state abbreviation, numeric state code, inventory year, and stand group |
| `raster_model` | Review forest/live status codes, stored-value multipliers, NoData settings, and provenance rasters |
| `segmentation` | Choose target seed acreage, size limits, variable weights, neighborhood settings, convergence, and similar-stand merging |
| `ownership` and ownership flags | Match class codes to your source and decide whether ownership constrains boundaries |
| `small_stands` | Decide how to handle units below your minimum size |
| `riparian` | Choose buffer rules, hydrography inputs, wetland classification, and management codes |
| `inventory` | Review donor threshold, live-tree status, species columns, and missing-inventory imputation |
| `database` | Choose database grouping and sampling design |
| `proximity_chain_qa` | Choose diagnostics and tolerances for proximity grouping |

These are scientific choices, not just installation settings. Record why they
fit your analysis rather than assuming the example values apply everywhere.

### Small stands

The available actions are `force_merge`, `retain`, `delete`, and
`proximity_merge`. Set the action in
`small_stands["remaining_small_stand_action"]`; that section takes precedence
over the legacy action key in `segmentation`.

Proximity grouping can combine eligible nearby components under a shared
management-unit ID without filling gaps in the geometry. Review distance,
ownership, forest type, maximum combined acreage, and unmatched-unit settings.
See [advanced behavior](advanced_behavior.md) for the identity and QA details.
If you choose deletion, quantify its acreage impact.

### Missing inventories

With `inventory.impute_missing_stands` enabled, units without usable direct
tree lists can receive a nearest runnable stand's inventory. Review the
imputed acreage and donor-distance reports to decide whether that is
appropriate for your use. Disabling imputation does not guarantee every
spatial unit will appear in the final runnable inventory.

### Database organization

- `grouped`: separate databases by `DB_GROUP`.
- `combined`: one database containing the runnable inventory.
- `both`: both forms; they contain overlapping inventories.

Keep sampling design consistent with weighted `TREE_COUNT` values and your
FVS template. Validate species translations against the selected FVS variant.

## Check and run your configuration

From `LETO_v1_project`:

```bat
python -c "from leto.config import load_config; c = load_config('configs/my_region.py'); print(c.region_id, c.run_name)"
python run_leto.py configs\my_region.py
```

The first command checks loading and selected rules. It does not verify all
data paths, fields, or GIS requirements. The second runs all stages.
Review [validation](validation.md) before using the resulting inventories.
