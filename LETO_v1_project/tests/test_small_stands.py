import pandas as pd

from leto.small_stands import (
    aggregate_management_unit_weights,
    plan_proximity_merges,
)


def test_proximity_merge_groups_without_changing_segment_ids():
    units = pd.DataFrame(
        {
            "SEG_ID": [10, 20, 30],
            "Acres": [3.0, 8.0, 20.0],
            "OWN_CODE": [4, 4, 4],
            "FORTYPCD_DOM": [101, 101, 101],
            "STDAGE_MEAN": [42.0, 43.0, 80.0],
        }
    )
    candidates = pd.DataFrame(
        {
            "SOURCE_SEG_ID": [10, 10],
            "TARGET_SEG_ID": [20, 30],
            "DISTANCE_METERS": [55.0, 10.0],
        }
    )
    mapping, merges, qa = plan_proximity_merges(
        units,
        candidates,
        minimum_acres=5.0,
        maximum_distance_meters=60.0,
        maximum_combined_acres=300.0,
        require_same_ownership=True,
        require_same_forest_type=True,
        retain_when_no_candidate=True,
    )

    assert set(mapping["SEG_ID"]) == {10, 20, 30}
    assert mapping.set_index("SEG_ID").loc[10, "MU_ID"] == 10
    assert mapping.set_index("SEG_ID").loc[20, "MU_ID"] == 10
    assert merges.iloc[0]["TARGET_SEG_ID"] == 20
    assert qa["proximity_merges"] == 1
    assert qa["total_acres_before"] == qa["total_acres_after"] == 31.0


def test_proximity_merge_respects_ownership_and_forest_type():
    units = pd.DataFrame(
        {
            "SEG_ID": [1, 2, 3],
            "Acres": [2.0, 10.0, 10.0],
            "OWN_CODE": [1, 2, 1],
            "FORTYPCD_DOM": [100, 100, 200],
            "STDAGE_MEAN": [20.0, 20.0, 20.0],
        }
    )
    candidates = pd.DataFrame(
        {
            "SOURCE_SEG_ID": [1, 1],
            "TARGET_SEG_ID": [2, 3],
            "DISTANCE_METERS": [5.0, 5.0],
        }
    )
    mapping, merges, qa = plan_proximity_merges(
        units,
        candidates,
        minimum_acres=5.0,
        maximum_distance_meters=60.0,
        maximum_combined_acres=300.0,
        require_same_ownership=True,
        require_same_forest_type=True,
        retain_when_no_candidate=True,
    )

    assert merges.empty
    assert mapping.set_index("SEG_ID").loc[1, "MU_ID"] == 1
    assert qa["small_stands_remaining"] == 1


def test_area_weighted_donor_aggregation_matches_example():
    parent_weights = pd.DataFrame(
        {
            "SEG_ID": [1, 1, 2],
            "PLT_CN": ["Plot1", "Plot2", "Plot3"],
            "WEIGHT": [0.70, 0.30, 1.00],
        }
    )
    components = pd.DataFrame(
        {
            "SEG_ID": [1, 2],
            "MU_ID": [1, 1],
            "Acres": [10.0, 5.0],
        }
    )

    result = aggregate_management_unit_weights(parent_weights, components)
    weights = result.set_index("PLT_CN")["WEIGHT"]

    assert abs(weights["Plot1"] - 0.4666666666666667) < 1e-12
    assert abs(weights["Plot2"] - 0.2) < 1e-12
    assert abs(weights["Plot3"] - 0.3333333333333333) < 1e-12
    assert abs(result["WEIGHT"].sum() - 1.0) < 1e-12
    assert result["DONOR_ACRES"].sum() == 15.0
