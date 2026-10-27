"""Donor-weight filtering without ArcGIS dependencies."""
from __future__ import annotations

import pandas as pd


def filter_weights(weights: pd.DataFrame, minimum_weight: float) -> pd.DataFrame:
    """Filter donor weights and renormalize within each SEG_ID.

    Requires SEG_ID, STAND_ID, PLT_CN, and WEIGHT columns. Invalid or
    nonpositive weights are discarded. If the threshold removes every donor
    for a segment, its largest donor is retained. Returns a sorted copy with
    WEIGHT_ORIGINAL and normalized WEIGHT; missing columns raise KeyError.
    """
    weights = weights.copy()
    weights["WEIGHT_ORIGINAL"] = pd.to_numeric(weights["WEIGHT"], errors="coerce")
    weights = weights.dropna(
        subset=["SEG_ID", "STAND_ID", "PLT_CN", "WEIGHT_ORIGINAL"]
    )
    weights = weights[weights["WEIGHT_ORIGINAL"] > 0].copy()
    retained = weights[weights["WEIGHT_ORIGINAL"] >= minimum_weight].copy()
    lost = set(weights["SEG_ID"]) - set(retained["SEG_ID"])
    if lost:
        largest = (
            weights[weights["SEG_ID"].isin(lost)]
            .sort_values(["SEG_ID", "WEIGHT_ORIGINAL"], ascending=[True, False])
            .drop_duplicates("SEG_ID")
        )
        retained = pd.concat([retained, largest], ignore_index=True)
    retained["WEIGHT"] = retained["WEIGHT_ORIGINAL"] / retained.groupby(
        "SEG_ID"
    )["WEIGHT_ORIGINAL"].transform("sum")
    return retained.sort_values(
        ["SEG_ID", "WEIGHT"], ascending=[True, False]
    ).reset_index(drop=True)
