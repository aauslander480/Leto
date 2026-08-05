"""Verified Maine test-area configuration.

Edit paths here. The LETO source files should not need path edits.
"""

CONFIG = {
    "region_id": "ME",
    "run_name": "Maine_LETO_Baseline",
    "working_gdb": (
        r"C:\Users\aca28742\OneDrive - University of Georgia\Documents"
        r"\ArcGIS\Projects\Unified_leto_test\Unified_leto_test.gdb"
    ),
    "runs_root_folder": r"D:\Misc\Unified_Leto_test(ME)",

    "study_area": r"D:\Misc\Maine_5_co\Maine_5_co.gdb\ME_5_Co_Dissolve",
    "treemap_raster": (
        r"D:\Misc\Whole_treemap\RDS-2025-0032 (1)"
        r"\Data\TreeMap2022_CONUS.tif"
    ),

    "fia_tree_files": [
        r"D:\Misc\Maine_5_co\ME_FVS_files\FIA_Tree\MA_TREE.csv",
        r"D:\Misc\Maine_5_co\ME_FVS_files\FIA_Tree\ME_TREE.csv",
        r"D:\Misc\Maine_5_co\ME_FVS_files\FIA_Tree\NH_TREE.csv",
        r"D:\Misc\Maine_5_co\ME_FVS_files\FIA_Tree\NY_TREE.csv",
    ],
    "fia_cond_files": [
        r"D:\Misc\Maine_5_co\ME_FVS_files\FIA_Cond\MA_COND.csv",
        r"D:\Misc\Maine_5_co\ME_FVS_files\FIA_Cond\ME_COND.csv",
        r"D:\Misc\Maine_5_co\ME_FVS_files\FIA_Cond\NH_COND.csv",
        r"D:\Misc\Maine_5_co\ME_FVS_files\FIA_Cond\NY_COND.csv",
    ],

    "species_crosswalk": r"D:\forest_condition_2026\FVS\FVS_SpeciesCrosswalk.xls",
    "species_crosswalk_sheet": "EasternSpeciesTranslator",
    "blank_fvs_database": r"D:\forest_condition_2026\FVS\FVS_Data.db",
    "fvs_variant": "NE",
    "fvs_state": "ME",
    "fvs_state_code": 23,
    "inventory_year": 2022,
    "stand_group_name": "All_Stands",

    "use_ownership_attributes": True,
    "use_ownership_hard_boundaries": True,
    "ownership_source_raster": (
        r"D:\forest_condition_2026\Stand_delineation"
        r"\RDS-2025-0045 (1)\Data\US_forest_ownership.tif"
    ),
    "collapse_nonforest_and_water": True,

    "segmentation": {
        "initial_seed_acres": 100.0,
        "minimum_stand_acres": 5.0,
        "maximum_stand_acres": 300.0,
        "maximum_iterations": 40,
        "convergence_threshold": 0.01,
        "minimum_score_improvement": 0.01,
        "use_eight_neighbors": False,
        "shared_edge_bonus": 0.2,
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
        # Baseline behavior. Proximity merging will be added after this
        # project reproduces the verified Maine result.
        "remaining_small_stand_action": "force_merge",
    },

    "riparian": {
        "create_riparian_stands": True,
        "classify_wetlands": False,
        "nwi_wetland_inputs": [],
        "wetland_majority_threshold": 50.0,
        "nhd_flowline": (
            r"D:\Misc\Stream_buffer_data\NHDPlus_H_National_Release_2_GDB"
            r"\NHDPlus_H_National_Release_2_GDB.gdb\NHDPlus_H_National"
            r"\NetworkNHDFlowline"
        ),
        "nhd_area": (
            r"D:\Misc\Stream_buffer_data\NHDPlus_H_National_Release_2_GDB"
            r"\NHDPlus_H_National_Release_2_GDB.gdb\NHDPlus_H_National"
            r"\NHDArea"
        ),
        "nhd_waterbody": (
            r"D:\Misc\Stream_buffer_data\NHDPlus_H_National_Release_2_GDB"
            r"\NHDPlus_H_National_Release_2_GDB.gdb\NHDPlus_H_National"
            r"\NHDWaterbody"
        ),
    },

    "inventory": {
        "minimum_plot_weight": 0.05,
        "impute_missing_stands": True,
    },

    "database": {
        "mode": "grouped",  # grouped, combined, or both
    },
}
