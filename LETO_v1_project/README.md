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

## Current small-stand behavior

The baseline still uses `remaining_small_stand_action = "force_merge"`. The planned 60-meter proximity merge will be introduced only after this reorganized project reproduces the verified baseline.
