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

## Running in an ArcGIS Pro notebook

Extract the entire ZIP. Then run:

```python
%run r"D:\path\to\LETO_v1_project\run_leto.py" r"D:\path\to\LETO_v1_project\configs\maine_test.py"
```

Edit paths and region-specific settings only in `configs/maine_test.py`.

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
