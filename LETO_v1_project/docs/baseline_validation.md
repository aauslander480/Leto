# Maine baseline validation

The first reorganized-project run should be compared against the verified unified run.

| Metric | Verified value |
|---|---:|
| Valid segmentation cells | 429,977 |
| Parent segments | 12,168 |
| Final management-unit polygon records | 19,662 |
| Direct FVS stands | 19,244 |
| Nearest-imputed stands | 418 |
| Grouped FVS databases | 6 |

Exact equality is expected because the stage algorithms were copied unchanged. Differences in output row ordering are harmless; differences in counts, acreage, weights, or database contents should be investigated before adding proximity merging.
