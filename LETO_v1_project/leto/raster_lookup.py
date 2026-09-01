"""Vectorized TreeMap VALUE-to-attribute assignment helpers."""
from __future__ import annotations

import numpy as np
import pandas as pd


REASON_NAMES = (
    "MISSING_FIA_PLT_CN",
    "FIA_ATTRIBUTE_MISSING",
    "NO_LIVE_TREES",
    "TREEMAP_RAT_ATTRIBUTE_MISSING",
)


def assign_records_to_raster(
    source_array: np.ndarray,
    value_records: dict[int, dict[str, object]],
    field: str,
    *,
    source_nodata: int | float,
    decimal_multiplier: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, int], int]:
    """Map TreeMap VALUE records to cells without per-code raster scans.

    The returned values and provenance are equivalent to assigning each
    ``source_array == VALUE`` mask independently, but the source raster is
    indexed once with a sorted lookup table. Missing-value reason counts are
    also accumulated without retaining four full-size boolean QA rasters.
    """
    source = np.asarray(source_array)
    if source.ndim != 2:
        raise ValueError("source_array must be two-dimensional.")
    if not value_records:
        raise ValueError("value_records cannot be empty.")

    source_valid = source != source_nodata
    valid_values = source[source_valid].astype(np.int64, copy=False)
    codes = np.asarray(sorted(int(code) for code in value_records), dtype=np.int64)
    positions = np.searchsorted(codes, valid_values)
    matched = positions < len(codes)
    if np.any(matched):
        matched[matched] &= codes[positions[matched]] == valid_values[matched]
    if not np.all(matched):
        missing_codes = np.unique(valid_values[~matched])[:20].tolist()
        raise RuntimeError(
            "TreeMap raster contains VALUE codes absent from the lookup: "
            f"{missing_codes}"
        )

    records = [value_records[int(code)] for code in codes]
    field_values = np.asarray(
        [
            np.nan if pd.isna(record.get(field))
            else float(record[field]) * float(decimal_multiplier)
            for record in records
        ],
        dtype=np.float64,
    )
    source_codes = np.asarray(
        [int(record.get(f"_SOURCE_{field}", 1)) for record in records],
        dtype=np.uint8,
    )

    working = np.full(source.shape, np.nan, dtype=np.float64)
    provenance = np.zeros(source.shape, dtype=np.uint8)
    mapped_values = field_values[positions]
    working[source_valid] = mapped_values
    direct = np.isfinite(mapped_values)
    valid_provenance = np.zeros(valid_values.shape, dtype=np.uint8)
    valid_provenance[direct] = source_codes[positions[direct]]
    provenance[source_valid] = valid_provenance

    missing_records = ~np.isfinite(field_values)
    reason_by_code = {
        "MISSING_FIA_PLT_CN": np.asarray([
            missing_records[index]
            and not bool(record.get("_FIA_PLOT_FOUND", False))
            for index, record in enumerate(records)
        ]),
        "FIA_ATTRIBUTE_MISSING": np.asarray([
            missing_records[index]
            and bool(record.get("_FIA_PLOT_FOUND", False))
            and not bool(record.get(f"_FIA_ATTRIBUTE_AVAILABLE_{field}", False))
            for index, record in enumerate(records)
        ]),
        "NO_LIVE_TREES": np.asarray([
            missing_records[index]
            and field in {"QMD", "TPA"}
            and bool(record.get("_FIA_PLOT_FOUND", False))
            and (
                pd.isna(record.get("LIVE_TREE_RECORDS"))
                or float(record.get("LIVE_TREE_RECORDS") or 0) <= 0
            )
            for index, record in enumerate(records)
        ]),
        "TREEMAP_RAT_ATTRIBUTE_MISSING": np.asarray([
            missing_records[index]
            and not bool(record.get(f"_RAT_ATTRIBUTE_AVAILABLE_{field}", False))
            for index, record in enumerate(records)
        ]),
    }
    reason_counts = {
        reason: int(np.count_nonzero(flags[positions]))
        for reason, flags in reason_by_code.items()
    }
    matched_value_codes = int(np.count_nonzero(np.isfinite(field_values)))
    return working, provenance, source_valid, reason_counts, matched_value_codes

