# Historical Maine regression reference

For general use, follow the [regional validation guide](validation.md).
This page records the original Maine test setup for historical regression
comparisons. These values apply only to its source data, study boundary, and
scientific settings; reproducing this example is not required for other areas.

| Metric | Verified value |
|---|---:|
| Valid segmentation cells | 429,977 |
| Parent segments | 12,168 |
| Final management-unit polygon records | 19,662 |
| Direct FVS stands | 19,244 |
| Nearest-imputed stands | 418 |
| Grouped FVS databases | 6 |

For the original setup with unchanged inputs and settings, exact equality is
expected when the configured small-stand action remains
`force_merge`. Differences in output row ordering are harmless; differences in
counts, acreage, weights, or database contents must be investigated.

The optional `proximity_merge` action is outside the verified baseline. Its QA
must show no acreage loss, normalized donor weights, no FVS stands without
trees, and valid grouped database outputs.
