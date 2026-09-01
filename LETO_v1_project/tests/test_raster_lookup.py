import numpy as np
import pytest

from leto.raster_lookup import assign_records_to_raster


def records():
    return {
        10: {
            "QMD": 1.2, "_SOURCE_QMD": 1, "_FIA_PLOT_FOUND": True,
            "_FIA_ATTRIBUTE_AVAILABLE_QMD": True,
            "_RAT_ATTRIBUTE_AVAILABLE_QMD": False,
            "LIVE_TREE_RECORDS": 3,
        },
        20: {
            "QMD": 2.5, "_SOURCE_QMD": 2, "_FIA_PLOT_FOUND": False,
            "_FIA_ATTRIBUTE_AVAILABLE_QMD": False,
            "_RAT_ATTRIBUTE_AVAILABLE_QMD": True,
            "LIVE_TREE_RECORDS": np.nan,
        },
        30: {
            "QMD": np.nan, "_SOURCE_QMD": 1, "_FIA_PLOT_FOUND": True,
            "_FIA_ATTRIBUTE_AVAILABLE_QMD": False,
            "_RAT_ATTRIBUTE_AVAILABLE_QMD": False,
            "LIVE_TREE_RECORDS": 0,
        },
    }


def test_vectorized_assignment_matches_values_provenance_and_reasons():
    source = np.array([[10, 20, -9999], [30, 10, 30]], dtype=np.int32)
    values, provenance, valid, reasons, matched = assign_records_to_raster(
        source, records(), "QMD", source_nodata=-9999, decimal_multiplier=10
    )
    assert values[0, 0] == pytest.approx(12)
    assert values[0, 1] == pytest.approx(25)
    assert np.isnan(values[1, 0])
    assert provenance.tolist() == [[1, 2, 0], [0, 1, 0]]
    assert valid.tolist() == [[True, True, False], [True, True, True]]
    assert reasons == {
        "MISSING_FIA_PLT_CN": 0,
        "FIA_ATTRIBUTE_MISSING": 2,
        "NO_LIVE_TREES": 2,
        "TREEMAP_RAT_ATTRIBUTE_MISSING": 2,
    }
    assert matched == 2


def test_unknown_source_value_fails_early():
    with pytest.raises(RuntimeError, match="absent from the lookup"):
        assign_records_to_raster(
            np.array([[10, 99]]), records(), "QMD", source_nodata=-9999
        )

