"""One-region LETO orchestration using the verified baseline algorithms."""
from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Any

import arcpy
import numpy as np
import pandas as pd

from .config import RunConfig
from .logging_utils import StageTiming, message, timed_stage
from .workspace import create_run_paths, ensure_working_gdb
from . import stage1_rasters as s1
from . import stage2_segmentation as s2
from . import stage3_crosswalk as s3
from . import stage4_inventory as s4
from . import stage5_database as s5


def _option(options: dict[str, Any], name: str, default: Any) -> Any:
    return options.get(name, default)


def _ensure_management_units(segment_polygons: str, output_fc: str) -> str:
    """Create final management-unit fields when riparian splitting is disabled."""
    if arcpy.Exists(output_fc):
        arcpy.management.Delete(output_fc)
    arcpy.management.CopyFeatures(segment_polygons, output_fc)
    fields = {field.name.upper() for field in arcpy.ListFields(output_fc)}
    additions = [
        ("MU_ID", "LONG", None), ("FINAL_ID", "LONG", None),
        ("ParentID", "TEXT", 40), ("Acres", "DOUBLE", None),
        ("RIPARIAN", "SHORT", None), ("WETLAND", "SHORT", None),
        ("MGMT_CLASS", "SHORT", None), ("MGMT_TYPE", "TEXT", 20),
        ("DB_GROUP", "TEXT", 30),
    ]
    for name, field_type, length in additions:
        if name.upper() not in fields:
            kwargs = {"field_length": length} if length else {}
            arcpy.management.AddField(output_fc, name, field_type, **kwargs)

    own_type = next(
        (f.name for f in arcpy.ListFields(output_fc) if f.name.upper() == "OWN_TYPE"),
        None,
    )
    cursor_fields = [
        "OID@", "SEG_ID", "MU_ID", "FINAL_ID", "ParentID", "Acres", "SHAPE@",
        "RIPARIAN", "WETLAND", "MGMT_CLASS", "MGMT_TYPE", "DB_GROUP",
    ]
    if own_type:
        cursor_fields.append(own_type)

    with arcpy.da.UpdateCursor(output_fc, cursor_fields) as cursor:
        for row in cursor:
            oid = int(row[0])
            seg_id = int(row[1])
            if row[2] is None:
                row[2] = oid
            row[3] = oid
            row[4] = f"SEG_{seg_id}"
            row[5] = float(row[6].getArea("GEODESIC", "ACRES"))
            row[7] = 0
            row[8] = 0
            row[9] = 0
            row[10] = "Upland"
            owner = row[12] if own_type else "Unknown"
            row[11] = owner or "Unknown"
            cursor.updateRow(row)
    return output_fc


def _filter_weights(weights: pd.DataFrame, minimum_weight: float) -> pd.DataFrame:
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


def _prepare_stands_dataframe(stands: pd.DataFrame, config: RunConfig) -> pd.DataFrame:
    stands = stands.copy()
    stands["SEG_ID"] = pd.to_numeric(stands["SEG_ID"], errors="raise").astype("int64")
    if "STAND_ID" not in stands.columns:
        stands["STAND_ID"] = "MU_" + stands["SEG_ID"].astype(str)
    stands["STAND_ID"] = stands["STAND_ID"].astype(str)
    stands["Stand_ID"] = stands["SEG_ID"].astype(str)
    stands["MU_ID"] = stands["SEG_ID"].astype(str)
    stands["Acres"] = pd.to_numeric(stands["Acres"], errors="coerce")
    if "AGE" not in stands.columns and "STDAGE_MEAN" in stands.columns:
        stands["AGE"] = np.floor(
            pd.to_numeric(stands["STDAGE_MEAN"], errors="coerce") + 0.5
        )
    stands["AGE"] = pd.to_numeric(stands["AGE"], errors="coerce").round().astype("Int64")
    stands["INV_YEAR"] = config.inventory_year
    stands["VARIANT"] = config.fvs_variant
    stands["STATE"] = config.fvs_state
    return stands


def _write_dataframe_outputs(
    paths: dict[str, str],
    stands: pd.DataFrame,
    weights: pd.DataFrame,
    final_stands: pd.DataFrame,
    trees: pd.DataFrame,
    missing: pd.DataFrame,
) -> None:
    weights.to_csv(Path(paths["crosswalk"]) / "Segment_PLT_CN_Weights_Filtered.csv", index=False)
    stands.to_csv(Path(paths["crosswalk"]) / "Segment_Stand_List.csv", index=False)
    final_stands.to_csv(Path(paths["fvs_inputs"]) / "FVS_StandInit.csv", index=False)
    trees.to_csv(Path(paths["fvs_inputs"]) / "FVS_TreeInit.csv", index=False)
    missing.to_csv(Path(paths["qa"]) / "Segment_FVS_Stands_No_Live_Trees.csv", index=False)


def _build_databases(
    config: RunConfig,
    paths: dict[str, str],
    final_stands: pd.DataFrame,
    tree_list: pd.DataFrame,
) -> pd.DataFrame:
    if not Path(config.blank_fvs_database).exists():
        raise FileNotFoundError(
            f"Blank FVS database not found:\n{config.blank_fvs_database}"
        )

    options = config.database
    mode = options["mode"]
    stands = final_stands.copy()
    if "DB_GROUP" not in stands.columns:
        stands["DB_GROUP"] = "All_Stands"
    stands["DB_GROUP"] = (
        stands["DB_GROUP"].fillna("Unknown").astype(str).str.strip().replace("", "Unknown")
    )

    groups: list[tuple[str, pd.DataFrame]] = []
    if mode in {"grouped", "both"}:
        groups.extend((name, frame.copy()) for name, frame in stands.groupby("DB_GROUP", sort=True))
    if mode in {"combined", "both"}:
        groups.append(("All_Stands", stands.copy()))
    if not groups:
        raise ValueError("database.mode must be 'grouped', 'combined', or 'both'.")

    s5.blank_db = config.blank_fvs_database
    s5.fvs_variant = config.fvs_variant
    s5.fvs_state_code = config.fvs_state_code
    s5.default_inv_year = config.inventory_year
    s5.stand_group_name = config.stand_group_name
    s5.basal_area_factor = options["basal_area_factor"]
    s5.inv_plot_size = options["inventory_plot_size"]
    s5.brk_dbh = options["breakpoint_dbh"]
    s5.num_plots = options["number_of_plots"]
    s5.nonstk_plots = options["nonstockable_plots"]

    summaries = []
    for group_name, group_stands in groups:
        ids = set(group_stands["STAND_ID"].astype(str))
        group_trees = tree_list[tree_list["STAND_ID"].astype(str).isin(ids)].copy()
        stand_out = s5.prepare_stands(group_stands)
        tree_out = s5.prepare_trees(group_trees)
        output_db = Path(paths["databases"]) / f"FVS_Inventory_{s5.safe_filename(group_name)}.db"
        s5.write_database(str(output_db), stand_out, tree_out)
        acres = pd.to_numeric(group_stands.get("Acres"), errors="coerce").sum()
        summaries.append({
            "DB_GROUP": group_name,
            "STANDS": len(stand_out),
            "ACRES": acres,
            "TREE_RECORDS": len(tree_out),
            "OUTPUT_DB": str(output_db),
        })
        message(f"Created database: {output_db}")

    summary = pd.DataFrame(summaries)
    summary.to_csv(Path(paths["run"]) / "FVS_Database_Group_Summary.csv", index=False)
    return summary


def _configure_stage2(config: RunConfig, paths: dict[str, str], rasters: dict[str, str]) -> None:
    seg = config.segmentation
    small = config.small_stands
    rip = config.riparian

    s2.age_raster = rasters["age_raster"]
    s2.forest_type_raster = rasters["forest_type_raster"]
    s2.basal_area_raster = rasters["basal_area_raster"]
    s2.qmd_raster = rasters["qmd_raster"]
    s2.tpa_raster = rasters["tpa_raster"]
    s2.output_folder = paths["segmentation"]
    s2.output_gdb = config.working_gdb
    s2.segment_raster = os.path.join(paths["segmentation"], "TreeMap_Segments_Initial.tif")
    s2.segment_polygons = os.path.join(config.working_gdb, "TreeMap_Segments_Initial")
    s2.segment_summary_csv = os.path.join(paths["segmentation"], "TreeMap_Segment_Summary.csv")
    s2.management_units_final = os.path.join(config.working_gdb, "ManagementUnits_Final")
    s2.riparian_split_stands = s2.management_units_final
    s2.riparian_buffer = os.path.join(config.working_gdb, "Riparian_Buffer")
    s2.ownership_aligned_raster = os.path.join(paths["segmentation"], "Ownership_Classes_Aligned.tif")
    s2.wetland_footprint = os.path.join(config.working_gdb, "NWI_Wetlands_Dissolved")
    s2.ownership_projected_temp = os.path.join(config.working_gdb, "_ownership_projected")

    s2.INITIAL_SEED_ACRES = seg["initial_seed_acres"]
    s2.MIN_STAND_ACRES = seg["minimum_stand_acres"]
    s2.MAX_STAND_ACRES = seg["maximum_stand_acres"]
    s2.MAX_ITERATIONS = seg["maximum_iterations"]
    s2.CONVERGENCE_THRESHOLD = seg["convergence_threshold"]
    s2.MINIMUM_SCORE_IMPROVEMENT = seg["minimum_score_improvement"]
    s2.USE_EIGHT_NEIGHBORS = seg["use_eight_neighbors"]
    s2.SHARED_EDGE_BONUS = seg["shared_edge_bonus"]
    s2.Z_CLIP = seg["standardization_clip"]
    s2.READ_NODATA = seg["read_nodata"]
    s2.REPORT_TIMING = seg["report_timing"]
    s2.WEIGHTS = dict(seg["variable_weights"])
    s2.MERGE_SIMILAR_STANDS = seg["merge_similar_stands"]
    s2.SIMILAR_MERGE_MAX_AGE_DIFFERENCE = seg["similar_merge_max_age_difference"]
    s2.SIMILAR_MERGE_REQUIRE_SAME_FORTYPCD = seg["similar_merge_require_same_fortypcd"]
    s2.SIMILAR_MERGE_MIN_SHARED_EDGES = seg["similar_merge_min_shared_edges"]
    s2.SIMILAR_MERGE_MAX_PASSES = seg["similar_merge_max_passes"]
    s2.REMAINING_SMALL_STAND_ACTION = small["remaining_small_stand_action"]
    s2.PROXIMITY_MAXIMUM_DISTANCE_METERS = small["maximum_distance_meters"]
    s2.PROXIMITY_REQUIRE_SAME_OWNERSHIP = small["require_same_ownership"]
    s2.PROXIMITY_REQUIRE_SAME_FOREST_TYPE = small["require_same_forest_type"]
    s2.PROXIMITY_MAXIMUM_COMBINED_ACRES = small["maximum_combined_acres"]
    s2.PROXIMITY_RETAIN_WHEN_NO_CANDIDATE = small["retain_when_no_candidate"]
    s2.PROXIMITY_QA_CSV = os.path.join(
        paths["qa"], "Proximity_Merge_Pairs.csv"
    )

    s2.USE_OWNERSHIP_HARD_BOUNDARIES = (
        config.use_ownership_attributes and config.use_ownership_hard_boundaries
    )
    s2.OWNERSHIP_SOURCE_RASTER = config.ownership_source_raster
    s2.COLLAPSE_NONFOREST_AND_WATER = config.collapse_nonforest_and_water
    ownership = config.ownership
    s2.OWNERSHIP_TYPES = dict(ownership["types"])
    s2.UNKNOWN_OWNERSHIP_CODE = ownership["unknown_code"]
    s2.UNKNOWN_OWNERSHIP_TYPE = ownership["unknown_type"]
    s2.NONFOREST_AND_WATER_CODES = tuple(ownership["nonforest_and_water_codes"])

    s2.CREATE_RIPARIAN_STANDS = rip["create_riparian_stands"]
    s2.CLASSIFY_WETLANDS = rip["classify_wetlands"]
    s2.NWI_WETLAND_INPUTS = list(rip["nwi_wetland_inputs"])
    s2.WETLAND_MAJORITY_THRESHOLD = rip["wetland_majority_threshold"]
    s2.FLOWLINE_BUFFER_RULES = dict(rip["flowline_buffer_rules"])
    s2.WATERBODY_BUFFER_RULES = dict(rip["waterbody_buffer_rules"])
    s2.AREA_BUFFER_RULES = dict(rip["area_buffer_rules"])
    management_codes = rip["management_codes"]
    s2.UPLAND_MANAGEMENT_CODE = management_codes["upland"]
    s2.RIPARIAN_MANAGEMENT_CODE = management_codes["riparian"]
    s2.WETLAND_MANAGEMENT_CODE = management_codes["wetland"]
    s2.nhd_flowline = rip["nhd_flowline"]
    s2.nhd_area = rip["nhd_area"]
    s2.nhd_waterbody = rip["nhd_waterbody"]


def run_pipeline(config: RunConfig) -> dict[str, Any]:
    """Run one independent LETO region from raw inputs to FVS databases."""
    started = datetime.now()
    timings: list[StageTiming] = []
    ensure_working_gdb(config.working_gdb)
    paths = create_run_paths(config.runs_root_folder, config.run_name)
    message(f"LETO run folder: {paths['run']}")
    config.write_json(Path(paths["run"]) / "Run_Settings.json")

    with timed_stage("Stage 1 - Build segmentation rasters", timings):
        raster_outputs = s1.run_stage({
            "treemap_raster": config.treemap_raster,
            "treemap_lookup_raster": config.treemap_lookup_raster,
            "fia_tree_files": config.fia_tree_files,
            "fia_cond_files": config.fia_cond_files,
            "output_folder": paths["rasters"],
            "working_gdb": config.working_gdb,
            "clip_boundary": config.study_area,
            **config.raster_model,
        })

    with timed_stage("Stage 2 - Segment landscape and build management units", timings):
        _configure_stage2(config, paths, raster_outputs)
        s2.main()
        create_riparian = config.riparian["create_riparian_stands"]
        management_units = (
            s2.management_units_final
            if create_riparian
            else _ensure_management_units(s2.segment_polygons, s2.management_units_final)
        )

    with timed_stage("Stage 3 - Build donor crosswalk", timings):
        raster_info = s3.validate_raster_alignment(s2.segment_raster, raster_outputs["source_raster"])
        acres_per_cell = s3.cell_area_acres(raster_info)
        required_values = s3.read_required_treemap_values(raster_outputs["source_raster"])
        value_to_plt = s3.read_treemap_value_lookup(config.treemap_lookup_raster, required_values)
        parent_weights = s3.build_weights_table(
            s2.segment_raster, raster_outputs["source_raster"], value_to_plt, acres_per_cell
        )
        s3.inventory_year = config.inventory_year
        s3.variant = config.fvs_variant
        s3.state = config.fvs_state
        stands_df = s3.read_stand_polygons(management_units)
        s3.validate_segment_ids(parent_weights, stands_df)
        weights_df = s3.expand_parent_weights_to_management_units(parent_weights, stands_df)
        weights_df.to_csv(Path(paths["crosswalk"]) / "Segment_PLT_CN_Weights.csv", index=False)

    with timed_stage("Stage 4 - Build weighted FVS inventory", timings):
        stands_df = _prepare_stands_dataframe(stands_df, config)
        minimum_weight = float(config.inventory["minimum_plot_weight"])
        weights_filtered = _filter_weights(weights_df, minimum_weight)
        s4.fia_tree_files = list(config.fia_tree_files)
        s4.species_crosswalk_xls = config.species_crosswalk
        s4.species_crosswalk_sheet = config.species_crosswalk_sheet
        s4.MIN_PLT_WEIGHT = minimum_weight
        s4.IMPUTE_MISSING_STANDS = bool(config.inventory["impute_missing_stands"])
        s4.LIVE_TREE_STATUS_CODE = config.inventory["live_tree_status_code"]
        s4.FIA_TREE_CSV_CHUNKSIZE = config.inventory["fia_tree_csv_chunksize"]
        s4.SPECIES_CROSSWALK_FIA_COLUMNS = tuple(
            config.inventory["species_crosswalk_fia_columns"]
        )
        s4.SPECIES_CROSSWALK_FVS_COLUMNS = tuple(
            config.inventory["species_crosswalk_fvs_columns"]
        )
        s4.NEAREST_SEARCH_RADIUS = config.inventory["nearest_search_radius"]
        s4.NEAREST_METHOD = config.inventory["nearest_method"]
        s4.stand_polygons = management_units
        s4.near_table = os.path.join(config.working_gdb, "missing_segment_to_nearest_runnable")
        s4.qa_missing_species_csv = os.path.join(paths["qa"], "QA_Missing_FVS_Species.csv")
        s4.qa_missing_species_by_stand_csv = os.path.join(paths["qa"], "QA_Missing_FVS_Species_By_Stand.csv")
        species_lookup = s4.read_species_lookup()
        tree_list, qa_metrics = s4.build_tree_list(weights_filtered, species_lookup)
        tree_list, still_missing = s4.impute_missing_tree_lists(stands_df, tree_list)
        final_stands = s4.build_final_stands(stands_df, weights_filtered, tree_list)
        _write_dataframe_outputs(
            paths, stands_df, weights_filtered, final_stands, tree_list, still_missing
        )

    with timed_stage("Stage 5 - Create FVS databases", timings):
        db_summary = _build_databases(config, paths, final_stands, tree_list)

    finished = datetime.now()
    timing_rows = [{"STAGE": item.name, "SECONDS": item.seconds} for item in timings]
    pd.DataFrame(timing_rows).to_csv(Path(paths["run"]) / "Stage_Timings.csv", index=False)

    summary = {
        "leto_version": "1.0.0-baseline",
        "region_id": config.region_id,
        "run_name": config.run_name,
        "started": started.isoformat(),
        "finished": finished.isoformat(),
        "elapsed_seconds": (finished - started).total_seconds(),
        "final_fvs_stands": int(len(final_stands)),
        "final_fvs_tree_records": int(len(tree_list)),
        "stands_still_missing_trees": int(len(still_missing)),
        "databases_created": int(len(db_summary)),
        "run_folder": paths["run"],
        "stage_timings": timing_rows,
        "qa_metrics": {
            **qa_metrics,
            "proximity_merge": dict(s2.PROXIMITY_MERGE_QA),
        },
    }
    Path(paths["run"], "Run_Summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    Path(paths["run"], "Run_Summary.txt").write_text(
        "\n".join([
            "LETO v1.0 BASELINE RUN SUMMARY",
            "=" * 72,
            f"Region: {config.region_id}",
            f"Run name: {config.run_name}",
            f"Started: {started}",
            f"Finished: {finished}",
            f"Elapsed: {finished - started}",
            "",
            f"Final FVS stands: {len(final_stands):,}",
            f"Final FVS tree records: {len(tree_list):,}",
            f"Stands still missing trees: {len(still_missing):,}",
            f"Databases created: {len(db_summary):,}",
            "",
            "STAGE TIMINGS",
            *[f"{item.name}: {item.seconds:,.1f} seconds" for item in timings],
            "",
            f"Run folder: {paths['run']}",
        ]),
        encoding="utf-8",
    )

    message("=" * 72)
    message("LETO v1.0 BASELINE PIPELINE COMPLETE")
    message("=" * 72)
    message(f"Run folder: {paths['run']}")

    return {
        "config": config,
        "paths": paths,
        "rasters": raster_outputs,
        "management_units": management_units,
        "stands": stands_df,
        "weights": weights_filtered,
        "fvs_stands": final_stands,
        "fvs_trees": tree_list,
        "missing_stands": still_missing,
        "database_summary": db_summary,
        "summary": summary,
    }
