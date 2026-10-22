import pandas as pd
import pytest

from leto.management_units import (
    assign_final_mu_ids,
    summarize_chain_run,
    summarize_mu_chain,
    validate_final_units,
)


COMPATIBILITY = ["MGMT_CLASS", "MGMT_TYPE", "DB_GROUP", "OWN_CODE", "RIPARIAN", "WETLAND"]


def _components():
    return pd.DataFrame({
        "SEG_ID": [10, 20, 30],
        "BIO_MU_ID": [10, 10, 30],
        "MU_ID": [10, 10, 30],
        "Acres": [2.0, 3.0, 8.0],
        "MGMT_CLASS": ["forest", "forest", "forest"],
        "MGMT_TYPE": ["upland", "upland", "upland"],
        "DB_GROUP": ["A", "A", "A"],
        "OWN_CODE": [1, 1, 1],
        "RIPARIAN": [0, 0, 0],
        "WETLAND": [0, 0, 0],
        "SOURCE_OID": [9, 3, 2],
    })


def test_final_mu_identity_preserves_bio_group_and_segment_ids():
    result = assign_final_mu_ids(
        _components(), compatibility_fields=COMPATIBILITY, proximity_enabled=True
    )
    assert result["SEG_ID"].tolist() == [10, 20, 30]
    assert result.loc[result["SEG_ID"].isin([10, 20]), "MU_ID"].nunique() == 1
    assert result.set_index("SEG_ID").loc[10, "MU_PART"] == 1
    assert result.set_index("SEG_ID").loc[20, "MU_PART"] == 2


def test_final_mu_splits_incompatible_management_classes():
    components = _components()
    components.loc[components["SEG_ID"] == 20, ["MGMT_CLASS", "RIPARIAN"]] = ["riparian", 1]
    result = assign_final_mu_ids(
        components, compatibility_fields=COMPATIBILITY, proximity_enabled=True
    )
    assert result.loc[result["BIO_MU_ID"] == 10, "MU_ID"].nunique() == 2
    assert result.groupby("MU_ID")["MGMT_CLASS"].nunique().max() == 1


def test_disabled_path_preserves_existing_mu_ids():
    components = _components()
    components["BIO_MU_ID"] = components["SEG_ID"]
    components["MU_ID"] = [101, 102, 103]
    result = assign_final_mu_ids(
        components, compatibility_fields=COMPATIBILITY, proximity_enabled=False
    )
    assert result["MU_ID"].tolist() == [101, 102, 103]


def test_chain_span_flags_transitive_group():
    components = _components().iloc[:3].copy()
    components["MU_ID"] = 10
    pairs = [
        {"boundary_distance_m": 30.0, "centroid_distance_m": 80.0},
        {"boundary_distance_m": 60.0, "centroid_distance_m": 100.0},
        {"boundary_distance_m": 120.0, "centroid_distance_m": 160.0},
    ]
    qa = summarize_mu_chain(
        components, pairs, maximum_distance_meters=60.0, warning_chain_ratio=1.5
    )
    assert qa["MU_SPAN_M"] == 120.0
    assert qa["CHAIN_FLAG"] is True
    assert qa["CHAIN_RATIO"] == 2.0
    assert qa["SOURCE_SEG_COUNT"] == 3
    run = summarize_chain_run(pd.DataFrame([qa]))
    assert run["chain_flagged_units"] == 1


def test_dissolved_unit_validation_accepts_one_row_for_many_components():
    components = _components().iloc[:2]
    units = pd.DataFrame({"MU_ID": [10], "Acres": [5.0]})
    validate_final_units(
        units, components, compatibility_fields=COMPATIBILITY
    )


def test_dissolved_unit_validation_rejects_duplicate_mu_rows():
    components = _components().iloc[:2]
    units = pd.DataFrame({"MU_ID": [10, 10], "Acres": [2.0, 3.0]})
    with pytest.raises(RuntimeError, match="one row per MU_ID"):
        validate_final_units(units, components, compatibility_fields=COMPATIBILITY)

