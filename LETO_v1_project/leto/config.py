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
    segmentation: dict[str, Any] = field(default_factory=dict)
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
    raw.setdefault("segmentation", {})
    raw.setdefault("riparian", {})
    raw.setdefault("inventory", {})
    raw.setdefault("database", {})

    return RunConfig(**raw)
