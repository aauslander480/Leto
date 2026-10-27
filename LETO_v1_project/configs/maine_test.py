"""Verified Maine test-area configuration.

Replace all C:\\path\\to placeholders with your local paths before running.
The LETO source files should not need path edits.
"""

CONFIG = {
    "region_id": "ME",
    "run_name": "Maine_LETO_Baseline",
    "working_gdb": r"C:\path\to\working.gdb",
    "runs_root_folder": r"C:\path\to\LETO_Runs",

    "study_area": r"C:\path\to\study_area.gdb\StudyBoundary",
    "treemap_raster": r"C:\path\to\TreeMap2022_CONUS.tif",

    "fia_tree_files": [
        r"C:\path\to\FIA\MA_TREE.csv",
        r"C:\path\to\FIA\ME_TREE.csv",
        r"C:\path\to\FIA\NH_TREE.csv",
        r"C:\path\to\FIA\NY_TREE.csv",
    ],
    "fia_cond_files": [
        r"C:\path\to\FIA\MA_COND.csv",
        r"C:\path\to\FIA\ME_COND.csv",
        r"C:\path\to\FIA\NH_COND.csv",
        r"C:\path\to\FIA\NY_COND.csv",
    ],

    "species_crosswalk": r"C:\path\to\FVS_SpeciesCrosswalk.xls",
    "species_crosswalk_sheet": "EasternSpeciesTranslator",
    "blank_fvs_database": r"C:\path\to\FVS_Data.db",
    "fvs_variant": "NE",
    "fvs_state": "ME",
    "fvs_state_code": 23,
    "inventory_year": 2022,
    "stand_group_name": "All_Stands",

    "use_ownership_attributes": True,
    "use_ownership_hard_boundaries": True,
    "ownership_source_raster": r"C:\path\to\US_forest_ownership.tif",
    "collapse_nonforest_and_water": True,

    # Regional scientific assumptions used to derive segmentation rasters.
    "raster_model": {
        "forest_condition_status_code": 1,
        "live_tree_status_code": 1,
        "source_nodata": -999999,
        "output_nodata": -9999,
        "coverage_check_nodata": -999999999,
        "write_qa_provenance_rasters": True,
        "stored_value_multipliers": {
            "STDAGE": 1.0,
            "FORTYPCD": 1.0,
            "BALIVE": 10.0,
            "QMD": 10.0,
            "TPA": 10.0,
        },
    },

    "ownership": {
        "types": {
            0: "Unknown",
            1: "NonForest",
            2: "Water",
            3: "Family",
            4: "Corporate",
            5: "Tribal",
            6: "Federal",
            7: "State",
            8: "Local",
        },
        "unknown_code": 0,
        "unknown_type": "Unknown",
        "nonforest_and_water_codes": (1, 2),
    },

    "segmentation": {
        "initial_seed_acres": 100.0,
        "minimum_stand_acres": 5.0,
        "maximum_stand_acres": 300.0,
        "maximum_iterations": 40,
        "convergence_threshold": 0.01,
        "minimum_score_improvement": 0.01,
        "use_eight_neighbors": False,
        "shared_edge_bonus": 0.2,
        "standardization_clip": 4.0,
        "read_nodata": -999999,
        "report_timing": True,
        "variable_weights": {
            "FORTYPCD": 0.30,
            "STDAGE": 0.25,
            "BALIVE": 0.20,
            "QMD": 0.15,
            "TPA": 0.10,
        },
        "merge_similar_stands": True,
        "similar_merge_max_age_difference": 10.0,
        "similar_merge_require_same_fortypcd": True,
        "similar_merge_min_shared_edges": 2,
        "similar_merge_max_passes": 50,
        # Legacy key retained for existing configurations.
        "remaining_small_stand_action": "force_merge",
    },

    "small_stands": {
        # Keep force_merge for exact Maine baseline reproduction. Change only
        # this value to "proximity_merge" to enable the optional grouping.
        "remaining_small_stand_action": "force_merge",
        "maximum_distance_meters": 60.0,
        "require_same_ownership": True,
        "require_same_forest_type": True,
        "maximum_combined_acres": 300.0,
        "retain_when_no_candidate": True,
    },

    "proximity_chain_qa": {
        "enabled": True,
        "warning_chain_ratio": 2.0,
        "write_component_pair_csv": False,
        "acreage_tolerance_acres": 0.01,
    },

    "riparian": {
        "create_riparian_stands": True,
        "classify_wetlands": False,
        "nwi_wetland_inputs": [],
        "wetland_majority_threshold": 50.0,
        "flowline_buffer_rules": {
            46000: ("Unclassified Stream", 0.0),
            46003: ("Intermittent Stream", 0.0),
            46006: ("Perennial Stream", 45.0),
            46007: ("Ephemeral Stream", 0.0),
        },
        "waterbody_buffer_rules": {
            39000: ("Unclassified Lake/Pond", 0.0),
            39001: ("Intermittent Lake/Pond", 0.0),
            39005: ("Intermittent Lake/Pond", 0.0),
            39006: ("Intermittent Lake/Pond", 0.0),
            39004: ("Perennial Lake/Pond", 45.0),
            39009: ("Perennial Lake/Pond", 45.0),
            39010: ("Perennial Lake/Pond", 45.0),
            39011: ("Perennial Lake/Pond", 45.0),
            39012: ("Perennial Lake/Pond", 45.0),
        },
        "area_buffer_rules": {
            46000: ("Unclassified Stream/River Area", 0.0),
            46003: ("Intermittent Stream/River Area", 0.0),
            46006: ("Perennial Stream/River Area", 45.0),
            46007: ("Ephemeral Stream/River Area", 0.0),
            53700: ("Area of Complex Channels", 45.0),
        },
        "management_codes": {
            "upland": 0,
            "riparian": 1,
            "wetland": 2,
        },
        "nhd_flowline": r"C:\path\to\hydrography.gdb\NHDPlus_H_National\NetworkNHDFlowline",
        "nhd_area": r"C:\path\to\hydrography.gdb\NHDPlus_H_National\NHDArea",
        "nhd_waterbody": r"C:\path\to\hydrography.gdb\NHDPlus_H_National\NHDWaterbody",
    },

    "inventory": {
        "minimum_plot_weight": 0.05,
        "impute_missing_stands": True,
        "live_tree_status_code": 1,
        "fia_tree_csv_chunksize": 250_000,
        "species_crosswalk_fia_columns": ("FIA CODE", "FIA_CODE", "FIACODE"),
        "species_crosswalk_fvs_columns": ("SN_Mapped_To", "SN MAPPED TO"),
        "nearest_search_radius": "",
        "nearest_method": "PLANAR",
    },

    "database": {
        "mode": "grouped",  # grouped, combined, or both
        # Verified design for TREE_COUNT expressed as trees per acre.
        "basal_area_factor": -1,
        "inventory_plot_size": 1,
        "breakpoint_dbh": 0,
        "number_of_plots": 1,
        "nonstockable_plots": 0,
    },
}
