# Input data and downloads

[Getting started](getting_started.md) · [Configuration](configuration.md)

Use the sources below to assemble inputs for your study area. Agency landing
pages are preferred over version-specific download filenames, so you can
choose coverage and read metadata before downloading.

Sources reviewed October 27, 2026. The links identify publishers and products;
large archives and the FVS installer were not downloaded or tested with LETO.

## Required inputs

| Input | Where to get it | What to select / configure |
|---|---|---|
| Study boundary | Your project's boundary, or [Census TIGER/Line downloads](https://www.census.gov/cgi-bin/geo/shapefiles/) for state/county boundaries | Select your geography, extract the shapefile, and prepare your study boundary in ArcGIS Pro. Set `study_area`. A custom project boundary does not need to come from Census. |
| TreeMap raster and lookup | [USFS TreeMap downloads and data dictionary](https://data.fs.usda.gov/geodata/rastergateway/treemap/), [TreeMap product page](https://research.fs.usda.gov/firelab/products/dataandtools/treemap-tree-level-model-united-states-forests), and [TreeMap2022 archive record](https://www.fs.usda.gov/rds/archive/catalog/RDS-2025-0032) | The supplied configuration uses TreeMap2022_CONUS. Obtain the main TreeMap dataset with its raster attribute table and donor-plot mapping. Set `treemap_raster` and, where needed, `treemap_lookup_raster`. |
| FIA TREE CSVs | [FIA DataMart](https://apps.fs.usda.gov/fia/datamart/datamart.html); [agency overview and user-guide links](https://research.fs.usda.gov/products/dataandtools/fia-datamart) | Choose raw CSV tables for the donor states and retain `<STATE>_TREE.csv`. Set `fia_tree_files`. |
| FIA COND CSVs | [FIA DataMart](https://apps.fs.usda.gov/fia/datamart/datamart.html) | Retain `<STATE>_COND.csv` for the donor states. Set `fia_cond_files`. FVS-ready SQLite downloads are a different product from the raw CSVs used by these stages. |
| FIA-to-FVS species crosswalk | [Official FVS user guides: Variant Overview Supplemental Information](https://www.fs.usda.gov/managing-land/forest-management/fvs/documents/guides); [FVS_SpeciesCrosswalk.xls workbook](https://www.fs.usda.gov/fmsc/ftp/fvs/docs/overviews/FVS_SpeciesCrosswalk.xls) | Download the Excel workbook, not the explanatory PDF. Set `species_crosswalk`, sheet, and configured FIA/FVS columns for your variant. |
| Blank FVS SQLite template | [Official FVS Complete Package](https://www.fs.usda.gov/fvs/software/complete.php); [FVS software landing page](https://www.fs.usda.gov/managing-land/forest-management/fvs/software); [Database Extension guide](https://www.fs.usda.gov/sites/default/files/forest-management/fvs-dbs-user-guide.pdf) | Obtain a blank SQLite template through the FVS distribution or FVS support, then check its schema against LETO's requirements below. Set `blank_fvs_database`. |

### TreeMap: keep the donor lookup

Individual attribute rasters such as basal area or age do not replace the main
TreeMap donor raster. LETO needs raster VALUE codes mapped to `PLT_CN`.
Keep the downloaded raster's associated attribute/metadata files together.
If clipping removes lookup attributes, retain the full original lookup raster.
The TreeMap gateway describes the packaged files and available products.

### FIA: cover the donors, not just the boundary

Download TREE and COND records for all states needed to cover the TreeMap donor
plots. Plot identifiers can refer to records outside your study state's
administrative boundary. Preserve `PLT_CN` and other control numbers as text.
Keep the source download dates and releases with your project; the newest FIA
export is not guaranteed to retain every donor record used by an older TreeMap
release.

### Species: select the actual variant mapping

The official workbook contains mappings for multiple eastern and western
variants. LETO's supplied example uses `EasternSpeciesTranslator` and
`SN_Mapped_To` column aliases; those names are not a universal choice.
Inspect the workbook and select the column appropriate to your FVS variant
through `inventory.species_crosswalk_fvs_columns`, along with the correct
sheet. Review untranslated-species QA after running.

### FVS template: compatibility must be checked

The example filename `FVS_Data.db` is a local path, not a verified standalone
download name. No standalone download of that exact example file was verified.
The Complete Package is the official acquisition route; if a suitable template
is not present, use the support entry on the
[FVS home page](https://www.fs.usda.gov/managing-land/forest-management/fvs).

LETO expects `FVS_StandInit`, `FVS_TreeInit`, and `FVS_PlotInit` tables
and appends the columns produced by `leto/stage5_database.py`. Check all column
names and types before using a template. The current FVS Database Extension
guide notes that some control-number columns are absent from its standard
blank template, while LETO writes control-number fields. Therefore an
unmodified current template must not be assumed to be a drop-in replacement
for the original project's database. Do not create an empty `.db` file or
rename an Access database to satisfy this input.

## Inputs required by enabled options

| Input | Where to get it | What to select / configure |
|---|---|---|
| Ownership raster | [USFS NWOS geospatial products](https://research.fs.usda.gov/programs/nwos) → Forest Ownership in the Conterminous United States, 2025 Geospatial Dataset; [RDS-2025-0045 archive record](https://www.fs.usda.gov/rds/archive/catalog/RDS-2025-0045) | The supplied example uses `US_forest_ownership.tif` from this product. Set `ownership_source_raster` and verify the code dictionary. Earlier ownership products have different categories. |
| NHD hydrography for riparian processing | [USGS hydrography downloads](https://www.usgs.gov/national-hydrography/access-national-hydrography-products); [NHDPlus HR product details](https://www.usgs.gov/national-hydrography/nhdplus-high-resolution) | Download vector file-geodatabase coverage for your area. Locate flowline, area, and waterbody feature classes and check their classification codes. Set `riparian.nhd_flowline`, `nhd_area`, and `nhd_waterbody`. |
| NWI wetlands for wetland classification | [USFWS NWI wetlands data and download access](https://www.fws.gov/program/national-wetlands-inventory/wetlands-data) | Follow Download Data or Wetlands Mapper for state/watershed polygon downloads. Set `riparian.nwi_wetland_inputs` when `classify_wetlands` is enabled. |

The supplied configuration enables ownership and riparian processing and
disables wetland classification. Review your own flags before deciding which
downloads to omit.

USGS retains NHD/NHDPlus downloads, but NHD maintenance ended in 2023 as the
agency moved to 3DHP. LETO's current rules target NHD-style inputs; a 3DHP
download should not be assumed to have an interchangeable schema.

## Before running

- Extract local datasets and point the configuration at files or feature
  classes, not ZIP archives or web map display layers.
- Confirm coverage, fields, code dictionaries, and coordinate systems.
- Record source URLs, product versions, download dates, and any preparation.
- Retain downloaded metadata and follow each source's attribution guidance.

Some USFS archive requests returned access errors during automated review.
Use the linked TreeMap or NWOS product page to reach the archive in a browser
if the direct catalog link is blocked. The crosswalk workbook is linked by the
current official FVS guide page; its binary content was not inspected here.
