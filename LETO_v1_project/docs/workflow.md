# LETO workflow

1. Build aligned TreeMap/FIA segmentation rasters.
2. Segment the landscape and create management units.
3. Build management-unit-to-PLT_CN donor weights.
4. Build weighted FVS stand and tree inventories.
5. Create FVS-ready SQLite databases.

Stages 3–5 pass pandas DataFrames internally. CSV files remain as audit outputs.
