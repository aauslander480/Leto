# Validating a regional run

[Getting started](getting_started.md) · [Workflow and outputs](workflow.md)

Check results against your study area and analysis objectives. Expected counts
depend on your inputs and configuration.

## Establish an accepted reference run

1. Start with a manageable study boundary and confirm input coverage,
   coordinate systems, donor states, species translations, and FVS variant.
2. Save the resolved configuration, source product versions, ArcGIS Pro/Python
   environment details, console output, run summaries, and QA reports.
3. Archive final spatial datasets from the working geodatabase as well as the
   timestamped result folder.
4. Review the checks below and document acceptable tolerances and any
   exclusions. Record why the run is suitable for the intended analysis.

## Acceptance checks

| Area | Review |
|---|---|
| Spatial coverage | Confirm the expected raster footprint and acreage through segmentation and management-unit creation. Investigate unexplained losses. |
| Management units | Distinguish component polygon counts from dissolved unit counts; confirm ownership and management routing. |
| Donors | Require weights to sum to 1 per final unit within numerical tolerance, including after filtering and acreage-weighted grouping. |
| Species and tree lists | Review untranslated species, weighted TPA losses, positive tree counts, and height fallback results. |
| Imputation | Review directly supported and imputed acreage, failure reasons, donor distances, and units still without usable trees. |
| Inventory integrity | Check unique stand IDs, unique tree IDs within stands, and that every tree references a valid stand. |
| Databases | Check template schema, sampling design, row counts, and grouping. Do not double-count the combined inventory in `both` mode. |

For proximity grouping, additionally review preserved component geometry and
provenance, combined acreage, compatibility splits, multipart units, and chain
distance diagnostics. For deliberate deletion, quantify removed acreage.

## Compare changes

Rerun an accepted configuration when changing code or dependencies. Compare
acreage, donor weights, stand/tree counts, imputation, and database contents;
row ordering alone is not a scientific difference. Investigate changes outside
documented tolerances.

When changing inputs or modeling assumptions, document the intended effect and
establish a new reference rather than expecting identical outputs.

## Historical regression example

The [Maine regression reference](baseline_validation.md) records the original
test setup. Its counts apply only to that configuration.
