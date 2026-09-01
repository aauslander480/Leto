"""Configuration loading and validation for LETO."""
from __future__ import annotations

import importlib.util
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from types import ModuleType
from typing import Any


@dataclass
class RunConfig:
    region_id: str
    run_name: str
    working_gdb: str
    runs_root_folder: str
    study_area: str
    treemap_raster: str
    treemap_lookup_raster: str
    fia_tree_files: list[str]
    fia_cond_files: list[str]
    species_crosswalk: str
    species_crosswalk_sheet: str
    blank_fvs_database: str
    fvs_variant: str
    fvs_state: str
    fvs_state_code: int
    inventory_year: int
    stand_group_name: str = "All_Stands"
    use_ownership_attributes: bool = True
    use_ownership_hard_boundaries: bool = True
    ownership_source_raster: str | None = None
    collapse_nonforest_and_water: bool = True
    raster_model: dict[str, Any] = field(default_factory=dict)
    ownership: dict[str, Any] = field(default_factory=dict)
    segmentation: dict[str, Any] = field(default_factory=dict)
    small_stands: dict[str, Any] = field(default_factory=dict)
    proximity_chain_qa: dict[str, Any] = field(default_factory=dict)
    riparian: dict[str, Any] = field(default_factory=dict)
    inventory: dict[str, Any] = field(default_factory=dict)
    database: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def write_json(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_dict(), indent=2),
            encoding="utf-8",
        )


def _load_module(path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location("leto_user_config", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load LETO configuration: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_config(path: str | Path) -> RunConfig:
    config_path = Path(path).expanduser().resolve()
    if not config_path.exists():
        raise FileNotFoundError(f"LETO configuration was not found:\n{config_path}")

    module = _load_module(config_path)
    raw = getattr(module, "CONFIG", None)
    if not isinstance(raw, dict):
        raise ValueError(
            f"Configuration must define a dictionary named CONFIG:\n{config_path}"
        )

    required = {
        "region_id", "run_name", "working_gdb", "runs_root_folder",
        "study_area", "treemap_raster", "fia_tree_files", "fia_cond_files",
        "species_crosswalk", "blank_fvs_database", "fvs_variant",
        "fvs_state", "fvs_state_code", "inventory_year",
    }
    missing = sorted(required - set(raw))
    if missing:
        raise ValueError(f"Configuration is missing required keys: {missing}")

    raw = dict(raw)
    raw.setdefault("treemap_lookup_raster", raw["treemap_raster"])
    raw.setdefault("species_crosswalk_sheet", "EasternSpeciesTranslator")
    raw.setdefault("stand_group_name", "All_Stands")
    raw.setdefault("use_ownership_attributes", True)
    raw.setdefault("use_ownership_hard_boundaries", True)
    raw.setdefault("ownership_source_raster", None)
    raw.setdefault("collapse_nonforest_and_water", True)
    raw.setdefault("raster_model", {})
    raw.setdefault("ownership", {})
    ownership_defaults = {
        "clean_treemap_nonforest": False,
        "nonforest_code": next(
            (
                int(code)
                for code, name in raw["ownership"].get("types", {}).items()
                if str(name).casefold() == "nonforest"
            ),
            None,
        ),
        "reassignment_excluded_codes": tuple(
            raw["ownership"].get("nonforest_and_water_codes", ())
        ) + (raw["ownership"].get("unknown_code"),),
    }
    raw["ownership"] = {
        **ownership_defaults,
        **raw["ownership"],
    }
    raw.setdefault("segmentation", {})
    small_stand_defaults = {
        "remaining_small_stand_action": raw["segmentation"].get(
            "remaining_small_stand_action", "retain"
        ),
        "maximum_distance_meters": 60.0,
        "require_same_ownership": True,
        "require_same_forest_type": True,
        "maximum_combined_acres": raw["segmentation"].get(
            "maximum_stand_acres", 300.0
        ),
        "retain_when_no_candidate": True,
    }
    supplied_small_stands = dict(raw.get("small_stands", {}))
    raw["small_stands"] = {
        **small_stand_defaults,
        **supplied_small_stands,
    }
    chain_qa_defaults = {
        "enabled": True,
        "warning_chain_ratio": 2.0,
        "write_component_pair_csv": False,
        "acreage_tolerance_acres": 0.01,
    }
    raw["proximity_chain_qa"] = {
        **chain_qa_defaults,
        **dict(raw.get("proximity_chain_qa", {})),
    }
    raw.setdefault("riparian", {})
    raw.setdefault("inventory", {})
    raw.setdefault("database", {})

    ownership_types = {
        int(code): str(name)
        for code, name in raw["ownership"].get("types", {}).items()
    }
    raw["ownership"]["types"] = ownership_types
    nonforest_code = raw["ownership"].get("nonforest_code")
    if raw["ownership"]["clean_treemap_nonforest"]:
        if nonforest_code is None or int(nonforest_code) not in ownership_types:
            raise ValueError(
                "ownership.nonforest_code must identify a configured ownership type "
                "when clean_treemap_nonforest is enabled."
            )
        nonforest_code = int(nonforest_code)
        if ownership_types[nonforest_code].casefold() != "nonforest":
            raise ValueError(
                "ownership.nonforest_code must map to the existing NonForest type."
            )
        raw["ownership"]["nonforest_code"] = nonforest_code
    excluded = tuple(
        sorted(
            {
                int(code)
                for code in raw["ownership"]["reassignment_excluded_codes"]
                if code is not None
            }
        )
    )
    unknown_code = int(raw["ownership"]["unknown_code"])
    required_exclusions = {
        unknown_code,
        *(int(code) for code in raw["ownership"]["nonforest_and_water_codes"]),
    }
    if not required_exclusions.issubset(excluded):
        raise ValueError(
            "ownership.reassignment_excluded_codes must include Unknown, "
            "NonForest, and Water codes."
        )
    raw["ownership"]["reassignment_excluded_codes"] = excluded

    action = str(
        raw["small_stands"]["remaining_small_stand_action"]
    ).lower().strip()
    allowed_actions = {"force_merge", "retain", "delete", "proximity_merge"}
    if action not in allowed_actions:
        raise ValueError(
            "small_stands.remaining_small_stand_action must be one of "
            f"{sorted(allowed_actions)}."
        )
    raw["small_stands"]["remaining_small_stand_action"] = action
    for key in ("maximum_distance_meters", "maximum_combined_acres"):
        try:
            value = float(raw["small_stands"][key])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(
                f"small_stands.{key} must be a positive number."
            ) from exc
        if value <= 0:
            raise ValueError(
                f"small_stands.{key} must be a positive number."
            )
        raw["small_stands"][key] = value
    minimum_acres = float(
        raw["segmentation"].get("minimum_stand_acres", 0)
    )
    if (
        action == "proximity_merge"
        and raw["small_stands"]["maximum_combined_acres"] < minimum_acres
    ):
        raise ValueError(
            "small_stands.maximum_combined_acres must be at least "
            "segmentation.minimum_stand_acres."
        )

    try:
        warning_ratio = float(raw["proximity_chain_qa"]["warning_chain_ratio"])
        acreage_tolerance = float(raw["proximity_chain_qa"]["acreage_tolerance_acres"])
    except (TypeError, ValueError) as exc:
        raise ValueError("proximity_chain_qa ratios and tolerances must be numeric.") from exc
    if warning_ratio <= 0 or acreage_tolerance < 0:
        raise ValueError(
            "proximity_chain_qa.warning_chain_ratio must be positive and "
            "acreage_tolerance_acres must be non-negative."
        )
    raw["proximity_chain_qa"]["warning_chain_ratio"] = warning_ratio
    raw["proximity_chain_qa"]["acreage_tolerance_acres"] = acreage_tolerance

    return RunConfig(**raw)

