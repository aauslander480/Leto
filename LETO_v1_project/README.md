# LETO v1.0 Baseline Project

This project reorganizes the verified unified LETO pipeline into a maintainable one-region software structure **without intentionally changing the scientific algorithms**.

## First objective

Run the same Maine test area and compare it with the verified baseline:

- 429,977 valid raster cells
- 12,168 parent segments
- 19,662 final management-unit polygon records
- 19,244 direct FVS stands
- 418 nearest-imputed stands
- 6 grouped FVS databases

Small numerical or ordering differences should be investigated before adding new algorithms.

## Folder structure

```text
LETO_v1_project/
├── run_leto.py
├── configs/
│   └── maine_test.py
├── leto/
│   ├── config.py
│   ├── pipeline.py
│   ├── workspace.py
│   ├── logging_utils.py
│   └── stage1...stage5 modules
├── tests/
└── docs/
```

## User documentation

Start with the [repository README](../README.md) and
[getting-started guide](docs/getting_started.md) for environment setup,
required inputs, and regional configuration. See [workflow and outputs](docs/workflow.md)
for results and QA, and [troubleshooting](docs/troubleshooting.md) for common failures.

## Running in an ArcGIS Pro notebook

Clone the repository or extract the complete source archive. Prepare the inputs
and copy/edit a regional configuration as described in the setup guide. Then run:

```python
%run "D:/Leto/LETO_v1_project/run_leto.py" "D:/Leto/LETO_v1_project/configs/my_region.py"
```

Edit paths and region-specific settings in your copied configuration; preserve
`configs/maine_test.py` as the baseline reference.

## Running another area

Copy `configs/maine_test.py`, rename it, and change its source paths and FVS metadata. The LETO package should remain unchanged.

## Small-stand behavior

The verified baseline continues to use
`small_stands["remaining_small_stand_action"] = "force_merge"`.

An optional `"proximity_merge"` mode groups eligible residual small segments
within the configured distance under one `MU_ID` without changing geometry or
renumbering their original `SEG_ID` values. Its default candidate rules require
matching ownership and dominant forest type, enforce the configured maximum
combined acreage, and use a 60-meter search distance. Donor weights for grouped
segments are combined by contributing acreage.


## Proximity-chain QA and final management-unit identity

Proximity grouping has two identities: `BIO_MU_ID` is the biological grouping
created by the unchanged proximity candidate algorithm, while final `MU_ID`
is the routing-safe stand identity after riparian, wetland, ownership, management
class, management type, and database-group compatibility is enforced.
`SEG_ID` is immutable provenance and `MU_PART` is deterministic by original
`SEG_ID` and source object ID.

Chain QA is calculated groupwise for each final MU. `MAX_PAIR_DIST_M` and
`MU_SPAN_M` are the maximum shortest boundary-to-boundary separation among
unique component pairs; `MAX_CENTROID_DIST_M` is reported separately.
`CHAIN_FLAG` is true when span exceeds the configured proximity search
distance, and `CHAIN_RATIO` is span divided by that distance. The groupwise
implementation intentionally avoids a global all-pairs distance table.


### Stage 2 and Stage 3 spatial contract

Stage 2 now writes both `ManagementSegments_Final` and
`ManagementUnits_Final`. The former remains singlepart and retains immutable
`SEG_ID`, `ParentSEG_ID`, final `MU_ID`, `BIO_MU_ID`, deterministic
`MU_PART`, acreage, routing classes, and biological provenance. The latter is
a multipart dissolve with exactly one feature per final `MU_ID`; biological
continuous fields are acreage-weighted and dominant forest type is selected by
acreage with a deterministic tie break.

Stage 3 reads `ManagementUnits_Final` for its one-row-per-FVS-stand list and
`ManagementSegments_Final` for donor provenance. Parent donor weights are
multiplied by contributing segment-piece acreage, divided by final MU acreage,
and repeated `PLT_CN` donors are consolidated before exact normalization.


## Diagnostic-only imputation QA

Stage 1 preserves the existing FIA-first, TreeMap-RAT-fallback, and independent
nearest-cell fill order. It now writes
`QA_STAGE1_Raster_Imputation_Summary.csv`,
`QA_STAGE1_Landscape_Imputation_Summary.csv`, five optional
`QA_<VARIABLE>_Source.tif` rasters, and
`QA_Any_Nearest_Imputation.tif`. Nearest distances are measured from each
filled cell to the exact donor indices already selected by the existing SciPy
algorithm; calculating distance does not participate in donor selection.

Stage 4 preserves the configured donor threshold and geographic
`GenerateNearTable` selection. It now writes
`QA_Imputation_Diagnostics.csv` with one row per imputed MU and
`QA_Imputation_Summary.csv` with overall, failure-reason, distance-statistic,
and distance-bin rows. A donor absent from all loaded FIA TREE files is reported
as `MISSING_FIA_DONOR`; a present donor with no configured live-tree records is
reported separately as `NO_LIVE_TREES`. Tree-data validity mirrors the
existing usable-tree filter rather than introducing new exclusions.
