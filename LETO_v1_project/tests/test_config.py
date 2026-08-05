from pathlib import Path

from leto.config import load_config


def test_example_config_loads():
    path = Path(__file__).parents[1] / "configs" / "maine_test.py"
    config = load_config(path)
    assert config.region_id == "ME"
    assert config.fvs_variant == "NE"
    assert config.segmentation["remaining_small_stand_action"] == "force_merge"
