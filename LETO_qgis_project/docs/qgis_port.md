# PyQGIS/GDAL translation contract

Source revision: `d2b92ee8809d54946339e4d4648244816bcd8212`.
The original `LETO_v1_project` remains unchanged. The copy retains the five
stages and scientific functions; `gis_backend.py` implements the spatial calls
with QGIS and GDAL. No Esri Python package is imported by the copy.

## Preserved methodology

The NumPy/SciPy segmentation, raster attribute derivation, FIA-first and
TreeMap-fallback order, independent nearest-cell fills, neighborhood rules,
seed sizes, objective weights, thresholds, convergence, stand-merging rules,
donor filtering/normalization, area-weighted MU aggregation, species mapping,
inventory sampling design, and SQLite generation are retained.
The Maine configuration is copied unchanged; the Sapelo2 configuration changes
paths/workspace identity only. AST regression checks compare non-GIS function
bodies against the sibling ArcPy project. Existing scientific tests are retained.

## Spatial substitutions

| ArcPy operation | QGIS/GDAL implementation | Validation requirement |
| --- | --- | --- |
| Raster clipping | GDAL rasterization with cell-center mask and outward source-grid snapping | Compare exact valid-cell masks at the boundary |
| Ownership alignment | GDAL Warp to explicit reference extent/width/height, nearest neighbor | Exact shape, affine transform, CRS and categorical cell values |
| NumPy raster conversion | GDAL GeoTIFF with north-up affine transform and explicit NoData | Exact array round-trip, including row orientation |
| RasterToPolygon NO_SIMPLIFY | GDAL Polygonize, four-connected, valid-data mask | Identical cell footprints, holes, areas, and parent SEG_ID |
| File-geodatabase workspace | Directory of individual GeoPackages, GeoTIFFs and CSV tables | Logical `workspace/name` resolves to `name.gpkg`, `.tif`, or `.csv` |
| Cursors and fields | QgsVectorLayer providers, transactions and explicit field types | Identifiers, NULLs and field values preserved; provider FIDs may differ |
| Clip/intersect/erase/dissolve | QGIS native processing with invalid geometry set to abort | Acreage and provenance; no automatic repair, simplification or gap filling |
| PLANAR buffers | QGIS buffer in model CRS, foot distances converted with 0.3048 | Curved boundaries use 64 segments per quadrant, not Esri's curve approximation |
| GEODESIC acreage | QgsDistanceArea on the input CRS ellipsoid | Compare acreage against ArcGIS within the original configured tolerance |
| PLANAR nearest donor | QgsSpatialIndex plus exact polygon geometry distances | Boundary distance, not centroid distance; preserve donor-selection rules |

The model grid must be projected in meters, consistent with the original
segmentation acreage calculation. Source hydrography is explicitly reprojected
to the model CRS for overlays/buffers. Geographic buffers are never silently
interpreted as degrees. Source file geodatabases are read via OGR, not written.
Copy TreeMap's `.tif.vat.dbf` and `.aux.xml` alongside the TIFF: its original
PLT_CN and biological fallback attributes are required.

## Explicit limitations

This is an experimental port, not certified equivalent to the ArcPy baseline.
GEODESIC nearest-boundary search (including optional `proximity_merge`) is
rejected before processing. There is no demonstrated equivalent implementation
yet; substituting PLANAR would change methodology. The Maine baseline uses
`force_merge` and PLANAR inventory imputation, so neither setting is changed.

Equal-distance nearest donors are exported to `<near_table>.ties.csv` and the
run stops. ArcGIS may choose a different tied feature; choosing the smallest
QGIS FID would silently change inventory composition. These cases require a
baseline donor crosswalk and an explicitly validated resolution strategy.

Buffer tessellation, overlay precision and feature ordering differ between GIS
engines. They may affect riparian slivers, acres, derived IDs and counts. Do not
loosen the original acreage tolerance or change scientific parameters to force
a pass. Output `Backend_Settings.json` records GIS versions, buffer resolution,
spatial conventions, limitations and validation status alongside Run_Settings.

## Validation before production

1. Run `python3 run_leto.py --check-runtime` and `python3 -m pytest -q` with
   QGIS's Python. A skipped `test_qgis_backend.py` means GIS tests did NOT run.
2. Run exactly the same Maine inputs and unchanged baseline scientific settings
   through both backends. Record dependency versions and source revision.
3. Compare all five stage-1 arrays, cell masks and provenance rasters exactly;
   compare parent label rasters cell-for-cell before polygon processing.
4. Match polygon components by parent SEG_ID, routing attributes and spatial
   footprint before comparing final MU IDs. Provider FIDs are not stable across
   engines. Check riparian/wetland classes, symmetric differences and acreage.
5. Compare donor PLT_CN provenance/weights, directly inventoried and imputed
   stands, species QA, weighted TPA, tree counts and all FVS database contents
   under a verified stand-ID mapping. Counts alone do not prove equivalence.
6. Investigate differences against `baseline_validation.md`; do not treat the
   new engine as validated until these checks pass on the actual source data.

Tests include genuine PyQGIS/GDAL raster, polygon, overlay, buffer, cursor,
geodesic-area and nearest-donor fixtures, plus the existing scientific tests.
The GitHub workflow installs QGIS on Linux to execute these tests. Full Maine
validation needs external datasets not present in this repository.

API references: [QGIS standalone Python](https://docs.qgis.org/3.40/en/docs/pyqgis_developer_cookbook/intro.html),
[QgsDistanceArea](https://qgis.org/pyqgis/master/core/QgsDistanceArea.html),
[GDAL Polygonize](https://gdal.org/en/stable/programs/gdal_polygonize.html).
