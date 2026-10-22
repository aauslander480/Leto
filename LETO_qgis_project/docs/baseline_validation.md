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

Exact equality is expected when the configured small-stand action remains
`force_merge`. Differences in output row ordering are harmless; differences in
counts, acreage, weights, or database contents must be investigated.

The optional `proximity_merge` action is outside the verified baseline. Its QA
must show no acreage loss, normalized donor weights, no FVS stands without
trees, and valid grouped database outputs.

