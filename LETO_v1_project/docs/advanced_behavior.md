# Advanced behavior

[User guide](../README.md) · [Configuration](configuration.md) · [Validation](validation.md)

## Small-stand behavior

See [Configuration](configuration.md) for small-stand actions and settings.

The optional `"proximity_merge"` mode groups eligible residual small segments
within the configured distance under one `MU_ID` without changing geometry or
renumbering their original `SEG_ID` values. Its default candidate rules require
matching ownership and dominant forest type, enforce the configured maximum
combined acreage, and use a 60-meter search distance. Donor weights for grouped
segments are combined by contributing acreage.

## Proximity-chain QA and final management-unit identity

Proximity grouping has two identities: `BIO_MU_ID` is the biological grouping
created by proximity grouping, while final `MU_ID`
is the stand identity after riparian, wetland, ownership, management
class, management type, and database-group compatibility is enforced.
`SEG_ID` is immutable provenance and `MU_PART` is deterministic by original
`SEG_ID` and source object ID.

Chain QA is calculated groupwise for each final MU. `MAX_PAIR_DIST_M` and
`MU_SPAN_M` are the maximum shortest boundary-to-boundary separation among
unique component pairs; `MAX_CENTROID_DIST_M` is reported separately.
`CHAIN_FLAG` is true when span exceeds the configured proximity search
distance, and `CHAIN_RATIO` is span divided by that distance. Distances are calculated within each group.

## Spatial outputs and donor weights

Stage 2 writes both `ManagementSegments_Final` and
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

## Imputation diagnostics

Stage 1 uses FIA attributes first, then TreeMap raster attributes, then an
independent nearest-cell fill for each variable. It writes
`QA_STAGE1_Raster_Imputation_Summary.csv`,
`QA_STAGE1_Landscape_Imputation_Summary.csv`, five optional
`QA_<VARIABLE>_Source.tif` rasters, and
`QA_Any_Nearest_Imputation.tif`. Nearest distances are measured from each
filled cell to the exact donor indices selected by SciPy; calculating distance does not participate in donor selection.

Stage 4 applies the configured donor threshold and uses `GenerateNearTable`
for nearest-stand selection. It writes
`QA_Imputation_Diagnostics.csv` with one row per imputed MU and
`QA_Imputation_Summary.csv` with overall, failure-reason, distance-statistic,
and distance-bin rows. A donor absent from all loaded FIA TREE files is reported
as `MISSING_FIA_DONOR`; a present donor with no configured live-tree records is
reported separately as `NO_LIVE_TREES`. Diagnostics use the same tree-validity rules as inventory preparation.
