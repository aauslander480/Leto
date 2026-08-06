from pathlib import Path

from leto.config import load_config


def test_example_config_loads():
    path = Path(__file__).parents[1] / "configs" / "maine_test.py"
    config = load_config(path)
    assert config.region_id == "ME"
    assert config.fvs_variant == "NE"
    assert config.segmentation["remaining_small_stand_action"] == "force_merge"
    assert config.raster_model["stored_value_multipliers"]["BALIVE"] == 10.0
    assert config.ownership["nonforest_and_water_codes"] == (1, 2)
    assert config.riparian["flowline_buffer_rules"][46006] == (
        "Perennial Stream",
        45.0,
    )
    assert config.riparian["management_codes"]["riparian"] == 1
    assert config.database["basal_area_factor"] == -1
    assert config.small_stands["remaining_small_stand_action"] == "force_merge"
    assert config.small_stands["maximum_distance_meters"] == 60.0
    assert config.small_stands["maximum_combined_acres"] == 300.0
