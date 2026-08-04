"""
04_build_fvs_inputs_from_segments.py

Convert final TreeMap stand segments and their weighted FIA donor plots into
FVS_StandInit.csv and FVS_TreeInit.csv.

Inputs created by 03_build_segment_fvs_crosswalk_v2.py
------------------------------------------------------
- Segment_PLT_CN_Weights.csv
- Segment_Stand_List.csv

Other inputs
------------
- Final stand polygon feature class
- FIA TREE.csv files
- FVS species crosswalk workbook

Main method
-----------
1. Remove donor plots below MIN_PLT_WEIGHT.
2. Re-normalize retained donor weights within each stand.
3. Join retained PLT_CN donors to live FIA trees.
4. Multiply FIA TPA_UNADJ by donor weight.
5. Translate FIA species codes to Southern FVS species codes.
6. Use segment STDAGE_MEAN-derived AGE from Segment_Stand_List.csv.
7. Impute a nearest runnable stand only if a segment has no usable live trees.
8. Write FVS-ready stand and tree CSV files.

Run this script in the ArcGIS Pro Python environment.
"""

import os

import arcpy
import numpy as np
import pandas as pd


# =============================================================================
# USER INPUTS
# =============================================================================

segment_folder = (
    r"D:\forest_condition_2026\Cellular_automata"
    r"\TreeMap_Segmentation_Output"
)

stand_polygons = (
    r"D:\forest_condition_2026\Cellular_automata"
    r"\TreeMap_Segmentation_Output\TreeMap_Segmentation.gdb"
    r"\ManagementUnits_Final"
)

segment_gdb = (
    r"D:\forest_condition_2026\Cellular_automata"
    r"\TreeMap_Segmentation_Output\TreeMap_Segmentation.gdb"
)

weights_csv = os.path.join(
    segment_folder,
    "Segment_PLT_CN_Weights.csv",
)

stand_list_csv = os.path.join(
    segment_folder,
    "Segment_Stand_List.csv",
)

# Final FVS CSV output folder
out_folder = r"D:\forest_condition_2026\FVS\FVS_Inputs"
os.makedirs(out_folder, exist_ok=True)

fvs_stands_csv = os.path.join(
    out_folder,
    "FVS_StandInit.csv",
)

fvs_trees_csv = os.path.join(
    out_folder,
    "FVS_TreeInit.csv",
)

missing_csv = os.path.join(
    out_folder,
    "Segment_FVS_Stands_No_Live_Trees.csv",
)

filtered_weights_csv = os.path.join(
    out_folder,
    "Segment_PLT_CN_Weights_Filtered.csv",
)

qa_missing_species_csv = os.path.join(
    out_folder,
    "QA_Missing_FVS_Species.csv",
)

qa_missing_species_by_stand_csv = os.path.join(
    out_folder,
    "QA_Missing_FVS_Species_By_Stand.csv",
)

qa_run_summary_txt = os.path.join(
    out_folder,
    "QA_FVS_Run_Summary.txt",
)

species_crosswalk_xls = (
    r"D:\forest_condition_2026\FVS\FVS_SpeciesCrosswalk.xls"
)

species_crosswalk_sheet = "EasternSpeciesTranslator"

fia_tree_files = [
    r"D:\forest_condition_2026\FVS\FIA_tree_CSV\FL_TREE.csv",
    r"D:\forest_condition_2026\FVS\FIA_tree_CSV\GA_TREE.csv",
    r"D:\forest_condition_2026\FVS\FIA_tree_CSV\AL_TREE.csv",
    r"D:\forest_condition_2026\FVS\FIA_tree_CSV\SC_TREE.csv",
]

# Donors below this proportion are removed, then the retained weights are
# re-normalized to sum to 1 within each stand.
MIN_PLT_WEIGHT = 0.05

INVENTORY_YEAR = 2022
FVS_VARIANT = "SN"
STATE = "FL"

# If True, a segment with no usable live trees borrows the complete weighted
# tree list from its nearest segment that does have usable trees.
IMPUTE_MISSING_STANDS = True

near_table = os.path.join(
    segment_gdb,
    "missing_segment_to_nearest_runnable",
)


# =============================================================================
# HELPERS
# =============================================================================

def message(text=""):
    print(text)
    try:
        arcpy.AddMessage(str(text))
    except Exception:
        pass


def require_file(path, label):
    if not os.path.exists(path):
        raise FileNotFoundError(f"{label} was not found:\n{path}")


def clean_identifier(series):
    """
    Normalize identifiers read from CSV files without converting them to
    floating-point scientific notation.
    """
    cleaned = series.astype(str).str.strip()
    cleaned = cleaned.str.replace(r"\.0$", "", regex=True)
    cleaned = cleaned.replace(
        {
            "nan": pd.NA,
            "None": pd.NA,
            "": pd.NA,
        }
    )
    return cleaned


def sql_in(field, values, numeric=True):
    values = [
        value for value in values
        if pd.notna(value)
    ]

    if not values:
        return "1 = 0"

    clauses = []

    for start in range(0, len(values), 900):
        chunk = values[start:start + 900]

        if numeric:
            formatted = ",".join(
                str(int(float(value)))
                for value in chunk
            )
        else:
            formatted = ",".join(
                "'" + str(value).replace("'", "''") + "'"
                for value in chunk
            )

        clauses.append(f"{field} IN ({formatted})")

    return " OR ".join(clauses)


def first_existing_column(frame, candidates):
    lookup = {
        column.upper(): column
        for column in frame.columns
    }

    for candidate in candidates:
        actual = lookup.get(candidate.upper())
        if actual:
            return actual

    return None


# =============================================================================
# READ AND FILTER DONOR WEIGHTS
# =============================================================================

def read_weights():
    require_file(weights_csv, "Segment donor-weight CSV")

    weights = pd.read_csv(weights_csv, dtype=str)

    required = {
        "SEG_ID",
        "STAND_ID",
        "PLT_CN",
        "WEIGHT",
    }

    missing = required - set(weights.columns)

    if missing:
        raise ValueError(
            "Segment_PLT_CN_Weights.csv is missing fields: "
            f"{sorted(missing)}"
        )

    weights["SEG_ID"] = pd.to_numeric(
        weights["SEG_ID"],
        errors="raise",
    ).astype("int64")

    weights["STAND_ID"] = clean_identifier(
        weights["STAND_ID"]
    )

    weights["PLT_CN"] = clean_identifier(
        weights["PLT_CN"]
    )

    weights["WEIGHT_ORIGINAL"] = pd.to_numeric(
        weights["WEIGHT"],
        errors="coerce",
    )

    weights = weights.dropna(
        subset=[
            "SEG_ID",
            "STAND_ID",
            "PLT_CN",
            "WEIGHT_ORIGINAL",
        ]
    ).copy()

    weights = weights[
        weights["WEIGHT_ORIGINAL"] > 0
    ].copy()

    message(f"Original donor rows: {len(weights):,}")
    message(
        f"Original stands: {weights['SEG_ID'].nunique():,}"
    )
    message(
        f"Original unique PLT_CN donors: "
        f"{weights['PLT_CN'].nunique():,}"
    )

    retained = weights[
        weights["WEIGHT_ORIGINAL"] >= MIN_PLT_WEIGHT
    ].copy()

    # A stand can theoretically have every donor below the cutoff. Preserve
    # its largest donor so no stand disappears solely because of filtering.
    retained_stands = set(retained["SEG_ID"])
    all_stands = set(weights["SEG_ID"])
    lost_stands = all_stands - retained_stands

    if lost_stands:
        largest_for_lost = (
            weights[weights["SEG_ID"].isin(lost_stands)]
            .sort_values(
                ["SEG_ID", "WEIGHT_ORIGINAL"],
                ascending=[True, False],
            )
            .drop_duplicates("SEG_ID")
        )

        retained = pd.concat(
            [retained, largest_for_lost],
            ignore_index=True,
        )

        message(
            "Stands with no donor meeting the threshold; "
            f"largest donor retained instead: {len(lost_stands):,}"
        )

    retained["WEIGHT"] = (
        retained["WEIGHT_ORIGINAL"]
        / retained.groupby("SEG_ID")[
            "WEIGHT_ORIGINAL"
        ].transform("sum")
    )

    retained = retained.sort_values(
        ["SEG_ID", "WEIGHT"],
        ascending=[True, False],
    ).reset_index(drop=True)

    weight_check = retained.groupby("SEG_ID")["WEIGHT"].sum()

    message(f"Filtered donor rows: {len(retained):,}")
    message(
        f"Filtered unique PLT_CN donors: "
        f"{retained['PLT_CN'].nunique():,}"
    )
    message(
        f"Average retained donors/stand: "
        f"{retained.groupby('SEG_ID').size().mean():.2f}"
    )
    message(
        f"Maximum retained donors/stand: "
        f"{retained.groupby('SEG_ID').size().max():,}"
    )
    message(
        f"Weight sums — minimum: {weight_check.min():.12f}, "
        f"mean: {weight_check.mean():.12f}, "
        f"maximum: {weight_check.max():.12f}"
    )

    retained.to_csv(
        filtered_weights_csv,
        index=False,
    )

    return retained


# =============================================================================
# READ STAND LIST
# =============================================================================

def read_stands():
    require_file(stand_list_csv, "Segment stand-list CSV")

    stands = pd.read_csv(
        stand_list_csv,
        dtype=str,
    )

    required = {
        "SEG_ID",
        "STAND_ID",
        "Acres",
    }

    missing = required - set(stands.columns)

    if missing:
        raise ValueError(
            "Segment_Stand_List.csv is missing fields: "
            f"{sorted(missing)}"
        )

    stands["SEG_ID"] = pd.to_numeric(
        stands["SEG_ID"],
        errors="raise",
    ).astype("int64")

    stands["STAND_ID"] = clean_identifier(
        stands["STAND_ID"]
    )

    stands["Acres"] = pd.to_numeric(
        stands["Acres"],
        errors="coerce",
    )

    # AGE was created from rounded STDAGE_MEAN in the previous script.
    if "AGE" in stands.columns:
        stands["AGE"] = pd.to_numeric(
            stands["AGE"],
            errors="coerce",
        )
    elif "STDAGE_MEAN" in stands.columns:
        stands["AGE"] = np.floor(
            pd.to_numeric(
                stands["STDAGE_MEAN"],
                errors="coerce",
            ) + 0.5
        )
    else:
        raise ValueError(
            "The stand list has neither AGE nor STDAGE_MEAN."
        )

    stands["AGE"] = stands["AGE"].round().astype("Int64")

    if stands["SEG_ID"].duplicated().any():
        raise ValueError(
            "Segment_Stand_List.csv contains duplicate SEG_ID values."
        )

    stands["Stand_ID"] = stands["SEG_ID"].astype(str)

    # Preserve MU_ID for compatibility with the existing database-creation
    # script. In this workflow, MU_ID is simply the numeric SEG_ID.
    stands["MU_ID"] = stands["SEG_ID"].astype(str)

    stands["INV_YEAR"] = INVENTORY_YEAR
    stands["VARIANT"] = FVS_VARIANT
    stands["STATE"] = STATE

    message(f"Stand records read: {len(stands):,}")
    message(
        f"Stands with AGE: {stands['AGE'].notna().sum():,}"
    )

    return stands


# =============================================================================
# SPECIES TRANSLATION
# =============================================================================

def read_species_lookup():
    require_file(
        species_crosswalk_xls,
        "FVS species crosswalk workbook",
    )

    crosswalk = pd.read_excel(
        species_crosswalk_xls,
        sheet_name=species_crosswalk_sheet,
        dtype=str,
    )

    fia_field = first_existing_column(
        crosswalk,
        ["FIA CODE", "FIA_CODE", "FIACODE"],
    )

    fvs_field = first_existing_column(
        crosswalk,
        ["SN_Mapped_To", "SN MAPPED TO"],
    )

    if fia_field is None or fvs_field is None:
        raise ValueError(
            "Could not find the FIA-code or SN species field in the "
            "species crosswalk workbook."
        )

    crosswalk["FIA_CODE_CLEAN"] = (
        pd.to_numeric(
            crosswalk[fia_field],
            errors="coerce",
        )
        .astype("Int64")
        .astype(str)
        .str.zfill(3)
    )

    lookup = dict(
        zip(
            crosswalk["FIA_CODE_CLEAN"],
            crosswalk[fvs_field],
        )
    )

    message(f"Species crosswalk entries: {len(lookup):,}")

    return lookup


# =============================================================================
# READ ONLY NEEDED FIA TREE RECORDS
# =============================================================================

def read_fia_trees(required_plt_cns):
    use_columns = [
        "CN",
        "PLT_CN",
        "STATUSCD",
        "SPCD",
        "DIA",
        "HT",
        "ACTUALHT",
        "CR",
        "TPA_UNADJ",
        "VOLCFNET",
    ]

    frames = []

    for path in fia_tree_files:
        require_file(path, "FIA TREE.csv")

        header = pd.read_csv(
            path,
            nrows=0,
        ).columns.tolist()

        available = [
            column for column in use_columns
            if column in header
        ]

        required_columns = {
            "CN",
            "PLT_CN",
            "STATUSCD",
            "SPCD",
            "DIA",
            "CR",
            "TPA_UNADJ",
        }

        missing_required = (
            required_columns - set(available)
        )

        if missing_required:
            raise ValueError(
                f"{path} is missing required FIA TREE fields: "
                f"{sorted(missing_required)}"
            )

        state_frames = []

        for chunk in pd.read_csv(
            path,
            dtype=str,
            usecols=available,
            chunksize=250_000,
            low_memory=False,
        ):
            chunk["PLT_CN"] = clean_identifier(
                chunk["PLT_CN"]
            )

            chunk = chunk[
                chunk["PLT_CN"].isin(required_plt_cns)
            ].copy()

            if not chunk.empty:
                state_frames.append(chunk)

        if state_frames:
            state_trees = pd.concat(
                state_frames,
                ignore_index=True,
            )
            frames.append(state_trees)
            message(
                f"Matched FIA trees from {os.path.basename(path)}: "
                f"{len(state_trees):,}"
            )
        else:
            message(
                f"Matched FIA trees from {os.path.basename(path)}: 0"
            )

    if not frames:
        raise RuntimeError(
            "None of the retained PLT_CN donors were found in the "
            "provided FIA TREE.csv files."
        )

    trees = pd.concat(
        frames,
        ignore_index=True,
    )

    trees["PLT_CN"] = clean_identifier(
        trees["PLT_CN"]
    )

    message(f"Combined matched FIA tree records: {len(trees):,}")
    message(
        f"Matched unique PLT_CN: {trees['PLT_CN'].nunique():,}"
    )

    return trees


# =============================================================================
# BUILD WEIGHTED FVS TREE LIST
# =============================================================================

def build_tree_list(weights, species_lookup):
    required_plt_cns = set(
        weights["PLT_CN"].dropna().astype(str)
    )

    trees = read_fia_trees(
        required_plt_cns
    )

    joined = weights.merge(
        trees,
        on="PLT_CN",
        how="left",
        indicator=True,
        validate="many_to_many",
    )

    unmatched = joined[
        joined["_merge"] == "left_only"
    ].copy()

    message(
        f"Unmatched retained donor rows: {len(unmatched):,}"
    )
    message(
        f"Unique unmatched PLT_CN: "
        f"{unmatched['PLT_CN'].nunique():,}"
    )

    if not unmatched.empty:
        message("Most frequent unmatched PLT_CN values:")
        message(
            unmatched["PLT_CN"]
            .value_counts()
            .head(20)
            .to_string()
        )

    joined = joined[
        joined["_merge"] == "both"
    ].copy()

    joined = joined[
        joined["STATUSCD"].astype(str) == "1"
    ].copy()

    joined["Species_FIA"] = joined["SPCD"]

    joined["Species_FIA_CLEAN"] = (
        pd.to_numeric(
            joined["Species_FIA"],
            errors="coerce",
        )
        .astype("Int64")
        .astype(str)
        .str.zfill(3)
    )

    joined["SPECIES"] = (
        joined["Species_FIA_CLEAN"]
        .map(species_lookup)
    )

    for column in [
        "DIA",
        "HT",
        "ACTUALHT",
        "CR",
        "TPA_UNADJ",
        "WEIGHT",
        "VOLCFNET",
    ]:
        if column in joined.columns:
            joined[column] = pd.to_numeric(
                joined[column],
                errors="coerce",
            )

    joined["DIAMETER"] = joined["DIA"]
    joined["CRRATIO"] = joined["CR"]
    joined["TREE_COUNT"] = (
        joined["TPA_UNADJ"]
        * joined["WEIGHT"]
    )

    if "HT" not in joined.columns:
        joined["HT"] = np.nan

    if "ACTUALHT" in joined.columns:
        joined["HT"] = joined["HT"].fillna(
            joined["ACTUALHT"]
        )

    joined["TREE_ID_FIA"] = clean_identifier(
        joined["CN"]
    )

    unusable_species = joined["SPECIES"].isna().sum()
    message(
        f"Live trees missing an FVS species translation: "
        f"{unusable_species:,}"
    )

    # ---------------------------------------------------------------------
    # Missing FVS species diagnostics
    # ---------------------------------------------------------------------

    message("=" * 72)
    message("MISSING FVS SPECIES DIAGNOSTICS")
    message("=" * 72)

    missing_species = joined[
        joined["SPECIES"].isna()
    ].copy()

    total_weighted_tpa_before_species_filter = (
        joined["TREE_COUNT"].sum()
    )

    missing_weighted_tpa = (
        missing_species["TREE_COUNT"].sum()
        if not missing_species.empty
        else 0.0
    )

    if total_weighted_tpa_before_species_filter > 0:
        missing_tpa_percent = (
            100.0
            * missing_weighted_tpa
            / total_weighted_tpa_before_species_filter
        )
    else:
        missing_tpa_percent = 0.0

    if not missing_species.empty:
        species_summary = (
            missing_species.groupby(
                "Species_FIA",
                dropna=False,
            )
            .agg(
                Tree_Records=("Species_FIA", "size"),
                Weighted_TPA=("TREE_COUNT", "sum"),
                Mean_Diameter=("DIAMETER", "mean"),
                Mean_Height=("HT", "mean"),
                Donor_Plots=("PLT_CN", "nunique"),
                Stands=("SEG_ID", "nunique"),
            )
            .reset_index()
            .sort_values(
                "Weighted_TPA",
                ascending=False,
            )
        )

        species_summary.to_csv(
            qa_missing_species_csv,
            index=False,
        )

        stand_summary = (
            missing_species.groupby("STAND_ID")
            .agg(
                Missing_TPA=("TREE_COUNT", "sum"),
                Missing_Trees=("TREE_COUNT", "size"),
                Missing_Species=("Species_FIA", "nunique"),
            )
            .reset_index()
            .sort_values(
                "Missing_TPA",
                ascending=False,
            )
        )

        stand_summary["Percent_Stand_TPA_Missing"] = (
            stand_summary["Missing_TPA"]
            / joined.groupby("STAND_ID")[
                "TREE_COUNT"
            ].sum().reindex(
                stand_summary["STAND_ID"]
            ).to_numpy()
            * 100.0
        )

        stand_summary.to_csv(
            qa_missing_species_by_stand_csv,
            index=False,
        )

        message(
            f"Missing weighted TPA: "
            f"{missing_weighted_tpa:,.2f}"
        )
        message(
            f"Percent of total weighted TPA: "
            f"{missing_tpa_percent:.2f}%"
        )
        message(
            f"Unique FIA species missing: "
            f"{missing_species['Species_FIA'].nunique():,}"
        )
        message("")
        message("Top missing FIA species:")
        message(
            species_summary.head(15).to_string(
                index=False
            )
        )
        message(
            f"Wrote: {qa_missing_species_csv}"
        )
        message(
            f"Wrote: "
            f"{qa_missing_species_by_stand_csv}"
        )

    else:
        # Write empty files with headers so every run has the same QA outputs.
        pd.DataFrame(
            columns=[
                "Species_FIA",
                "Tree_Records",
                "Weighted_TPA",
                "Mean_Diameter",
                "Mean_Height",
                "Donor_Plots",
                "Stands",
            ]
        ).to_csv(
            qa_missing_species_csv,
            index=False,
        )

        pd.DataFrame(
            columns=[
                "STAND_ID",
                "Missing_TPA",
                "Missing_Trees",
                "Missing_Species",
                "Percent_Stand_TPA_Missing",
            ]
        ).to_csv(
            qa_missing_species_by_stand_csv,
            index=False,
        )

        message(
            "No live FIA trees were missing an FVS species translation."
        )

    joined = joined.dropna(
        subset=[
            "STAND_ID",
            "SEG_ID",
            "PLT_CN",
            "SPECIES",
            "DIAMETER",
            "TREE_COUNT",
        ]
    ).copy()

    joined = joined[
        joined["TREE_COUNT"] > 0
    ].copy()

    # Create unique, sequential TREE_ID values within each final stand.
    joined = joined.sort_values(
        [
            "SEG_ID",
            "PLT_CN",
            "TREE_ID_FIA",
        ]
    ).reset_index(drop=True)

    joined["TREE_ID"] = (
        joined.groupby("STAND_ID").cumcount() + 1
    )

    final_columns = [
        "STAND_ID",
        "TREE_ID",
        "SPECIES",
        "DIAMETER",
        "HT",
        "CRRATIO",
        "TREE_COUNT",
        "MU_ID",
        "SEG_ID",
        "PLT_CN",
        "WEIGHT",
        "WEIGHT_ORIGINAL",
        "Species_FIA",
        "TREE_ID_FIA",
    ]

    # MU_ID is included for compatibility with the established LETO/FVS
    # database workflow and is equivalent to SEG_ID.
    joined["MU_ID"] = joined["SEG_ID"].astype(str)

    if "VOLCFNET" in joined.columns:
        final_columns.append("VOLCFNET")

    final = joined[final_columns].copy()

    final["TREE_SOURCE"] = "FIA_WEIGHTED_DIRECT"
    final["DONOR_STAND_ID"] = ""
    final["NEAR_DIST"] = np.nan

    message(f"Usable weighted tree records: {len(final):,}")
    message(
        f"Stands with direct live trees: "
        f"{final['STAND_ID'].nunique():,}"
    )

    qa_metrics = {
        "fia_tree_records_matched": int(len(trees)),
        "fia_donor_plots_matched": int(
            trees["PLT_CN"].nunique()
        ),
        "untranslated_live_tree_records": int(
            unusable_species
        ),
        "untranslated_fia_species": int(
            missing_species["Species_FIA"].nunique()
            if not missing_species.empty
            else 0
        ),
        "weighted_tpa_before_species_filter": float(
            total_weighted_tpa_before_species_filter
        ),
        "missing_weighted_tpa": float(
            missing_weighted_tpa
        ),
        "missing_weighted_tpa_percent": float(
            missing_tpa_percent
        ),
        "usable_weighted_tree_records": int(
            len(final)
        ),
        "direct_stands": int(
            final["STAND_ID"].nunique()
        ),
    }

    return final, qa_metrics


# =============================================================================
# NEAREST-STAND IMPUTATION
# =============================================================================

def impute_missing_tree_lists(
    stands,
    tree_list,
):
    direct_stands = set(
        tree_list["STAND_ID"].astype(str)
    )

    missing = stands[
        ~stands["STAND_ID"].astype(str).isin(
            direct_stands
        )
    ].copy()

    message(
        f"Stands lacking a direct usable live-tree list: "
        f"{len(missing):,}"
    )

    if (
        missing.empty
        or not IMPUTE_MISSING_STANDS
    ):
        return tree_list, missing

    runnable = stands[
        stands["STAND_ID"].astype(str).isin(
            direct_stands
        )
    ].copy()

    if runnable.empty:
        raise RuntimeError(
            "No runnable stands are available for nearest-neighbor "
            "tree-list imputation."
        )

    arcpy.management.MakeFeatureLayer(
        stand_polygons,
        "missing_segment_lyr",
    )

    arcpy.management.MakeFeatureLayer(
        stand_polygons,
        "runnable_segment_lyr",
    )

    arcpy.management.SelectLayerByAttribute(
        "missing_segment_lyr",
        "NEW_SELECTION",
        sql_in(
            "MU_ID",
            missing["SEG_ID"].tolist(),
            numeric=True,
        ),
    )

    arcpy.management.SelectLayerByAttribute(
        "runnable_segment_lyr",
        "NEW_SELECTION",
        sql_in(
            "MU_ID",
            runnable["SEG_ID"].tolist(),
            numeric=True,
        ),
    )

    if arcpy.Exists(near_table):
        arcpy.management.Delete(near_table)

    missing_selected = int(
        arcpy.management.GetCount(
            "missing_segment_lyr"
        )[0]
    )

    runnable_selected = int(
        arcpy.management.GetCount(
            "runnable_segment_lyr"
        )[0]
    )

    message(
        f"Missing polygons selected: {missing_selected:,}"
    )
    message(
        f"Runnable polygons selected: {runnable_selected:,}"
    )

    arcpy.analysis.GenerateNearTable(
        in_features="missing_segment_lyr",
        near_features="runnable_segment_lyr",
        out_table=near_table,
        search_radius="",
        location="NO_LOCATION",
        angle="NO_ANGLE",
        closest="CLOSEST",
        closest_count=1,
        method="PLANAR",
    )

    oid_to_segment = {}

    with arcpy.da.SearchCursor(
        stand_polygons,
        ["OID@", "MU_ID"],
    ) as cursor:
        for oid, seg_id in cursor:
            oid_to_segment[int(oid)] = int(seg_id)

    near_rows = []

    with arcpy.da.SearchCursor(
        near_table,
        [
            "IN_FID",
            "NEAR_FID",
            "NEAR_DIST",
        ],
    ) as cursor:
        for in_fid, near_fid, distance in cursor:
            near_rows.append(
                {
                    "MISSING_SEG_ID": (
                        oid_to_segment[int(in_fid)]
                    ),
                    "DONOR_SEG_ID": (
                        oid_to_segment[int(near_fid)]
                    ),
                    "NEAR_DIST": distance,
                }
            )

    nearest = pd.DataFrame(near_rows)

    if nearest.empty:
        message(
            "No nearest-neighbor pairs were generated."
        )
        return tree_list, missing

    stand_id_lookup = dict(
        zip(
            stands["SEG_ID"],
            stands["STAND_ID"],
        )
    )

    nearest["MISSING_STAND_ID"] = (
        nearest["MISSING_SEG_ID"]
        .map(stand_id_lookup)
    )

    nearest["DONOR_STAND_ID"] = (
        nearest["DONOR_SEG_ID"]
        .map(stand_id_lookup)
    )

    imputed_frames = []

    for row in nearest.itertuples(index=False):
        donor_trees = tree_list[
            tree_list["STAND_ID"].astype(str)
            == str(row.DONOR_STAND_ID)
        ].copy()

        if donor_trees.empty:
            continue

        donor_trees["DONOR_STAND_ID"] = (
            row.DONOR_STAND_ID
        )
        donor_trees["TREE_SOURCE"] = (
            "IMPUTED_NEAREST"
        )
        donor_trees["NEAR_DIST"] = row.NEAR_DIST
        donor_trees["STAND_ID"] = (
            row.MISSING_STAND_ID
        )
        donor_trees["SEG_ID"] = row.MISSING_SEG_ID
        donor_trees["MU_ID"] = str(
            row.MISSING_SEG_ID
        )

        donor_trees["TREE_ID"] = range(
            1,
            len(donor_trees) + 1,
        )

        imputed_frames.append(donor_trees)

    if imputed_frames:
        imputed = pd.concat(
            imputed_frames,
            ignore_index=True,
        )

        tree_list = pd.concat(
            [tree_list, imputed],
            ignore_index=True,
        )

        message(
            f"Nearest-imputed stands: "
            f"{imputed['STAND_ID'].nunique():,}"
        )

    final_stand_ids = set(
        tree_list["STAND_ID"].astype(str)
    )

    still_missing = stands[
        ~stands["STAND_ID"].astype(str).isin(
            final_stand_ids
        )
    ].copy()

    return tree_list, still_missing


# =============================================================================
# BUILD FINAL STAND TABLE
# =============================================================================

def build_final_stands(
    stands,
    weights,
    tree_list,
):
    # The largest retained donor is included as the stand's representative
    # PLT_CN for QA and backwards compatibility. Tree expansion still uses
    # every retained donor and its normalized weight.
    representative = (
        weights.sort_values(
            ["SEG_ID", "WEIGHT"],
            ascending=[True, False],
        )
        .drop_duplicates("SEG_ID")
        [["SEG_ID", "PLT_CN"]]
        .rename(
            columns={
                "PLT_CN": "PLT_CN",
            }
        )
    )

    stands = stands.merge(
        representative,
        on="SEG_ID",
        how="left",
        validate="one_to_one",
    )

    runnable_stand_ids = set(
        tree_list["STAND_ID"].astype(str)
    )

    final = stands[
        stands["STAND_ID"].astype(str).isin(
            runnable_stand_ids
        )
    ].copy()

    preferred = [
        "STAND_ID",
        "VARIANT",
        "INV_YEAR",
        "STATE",
        "Stand_ID",
        "MU_ID",
        "SEG_ID",
        "ParentSEG_ID",
        "Acres",
        "OWN_CODE",
        "OWN_TYPE",
        "OWN_RAW",
        "OWN_PCT",
        "OWN_MIXED",
        "MGMT_CLASS",
        "MGMT_TYPE",
        "DB_GROUP",
        "PLT_CN",
        "AGE",
    ]

    optional = [
        "STDAGE_MEAN",
        "STDAGE_SD",
        "STDAGE_MIN",
        "STDAGE_MAX",
        "STDAGE_HOM",
        "FORTYPCD_DOM",
        "FORTYPE_PCT",
        "BALIVE_MEAN",
        "BALIVE_SD",
        "BALIVE_HOM",
        "QMD_MEAN",
        "QMD_SD",
        "QMD_HOM",
        "TPA_MEAN",
        "TPA_SD",
        "TPA_HOM",
        "HOMOGENEITY",
    ]

    columns = [
        column
        for column in preferred + optional
        if column in final.columns
    ]

    final = final[columns].copy()

    return final


# =============================================================================
# MAIN
# =============================================================================

def main():
    arcpy.env.overwriteOutput = True

    if not arcpy.Exists(stand_polygons):
        raise FileNotFoundError(
            "Final segment polygon feature class was not found:\n"
            f"{stand_polygons}"
        )

    message("=" * 72)
    message("READING FINAL SEGMENTS AND DONOR WEIGHTS")
    message("=" * 72)

    stands = read_stands()
    weights = read_weights()

    stand_ids = set(stands["SEG_ID"])
    weight_ids = set(weights["SEG_ID"])

    if stand_ids != weight_ids:
        raise RuntimeError(
            "SEG_ID mismatch after donor filtering.\n"
            f"Only in stand list: "
            f"{sorted(stand_ids - weight_ids)[:20]}\n"
            f"Only in weights: "
            f"{sorted(weight_ids - stand_ids)[:20]}"
        )

    message("=" * 72)
    message("READING SPECIES CROSSWALK")
    message("=" * 72)

    species_lookup = read_species_lookup()

    message("=" * 72)
    message("BUILDING WEIGHTED FIA TREE LISTS")
    message("=" * 72)

    tree_list, qa_metrics = build_tree_list(
        weights,
        species_lookup,
    )

    message("=" * 72)
    message("CHECKING FOR STANDS WITHOUT USABLE TREES")
    message("=" * 72)

    tree_list, still_missing = (
        impute_missing_tree_lists(
            stands,
            tree_list,
        )
    )

    final_stands = build_final_stands(
        stands,
        weights,
        tree_list,
    )

    message("=" * 72)
    message("WRITING FVS CSV FILES")
    message("=" * 72)

    final_stands.to_csv(
        fvs_stands_csv,
        index=False,
    )

    tree_list.to_csv(
        fvs_trees_csv,
        index=False,
    )

    still_missing.to_csv(
        missing_csv,
        index=False,
    )

    message(f"Wrote: {fvs_stands_csv}")
    message(f"Wrote: {fvs_trees_csv}")
    message(f"Wrote: {missing_csv}")
    message(f"Wrote: {filtered_weights_csv}")

    # ---------------------------------------------------------------------
    # Write permanent QA summary
    # ---------------------------------------------------------------------

    retained_donors_per_stand = (
        weights.groupby("SEG_ID").size()
    )

    direct_stands_count = int(
        tree_list.loc[
            tree_list["TREE_SOURCE"]
            == "FIA_WEIGHTED_DIRECT",
            "STAND_ID",
        ].nunique()
    )

    imputed_stands_count = int(
        tree_list.loc[
            tree_list["TREE_SOURCE"]
            == "IMPUTED_NEAREST",
            "STAND_ID",
        ].nunique()
    )

    stand_total_tpa = (
        tree_list.groupby("STAND_ID")[
            "TREE_COUNT"
        ].sum()
    )

    qa_lines = [
        "LETO / FVS INPUT QA SUMMARY",
        "=" * 72,
        "",
        "STANDS AND DONORS",
        f"Final stands: {len(final_stands):,}",
                f"Retained donor rows: {len(weights):,}",
        f"Retained unique PLT_CN donors: {weights['PLT_CN'].nunique():,}",
        f"Average retained donors per stand: {retained_donors_per_stand.mean():.2f}",
        f"Maximum retained donors per stand: {retained_donors_per_stand.max():,}",
        "",
        "FIA TREE MATCHING",
        f"Matched FIA tree records: {qa_metrics['fia_tree_records_matched']:,}",
        f"Matched FIA donor plots: {qa_metrics['fia_donor_plots_matched']:,}",
        f"Untranslated live tree records: {qa_metrics['untranslated_live_tree_records']:,}",
        f"Unique untranslated FIA species: {qa_metrics['untranslated_fia_species']:,}",
        "",
        "MISSING-SPECIES IMPACT",
        f"Weighted TPA before species filtering: {qa_metrics['weighted_tpa_before_species_filter']:,.2f}",
        f"Missing weighted TPA: {qa_metrics['missing_weighted_tpa']:,.2f}",
        f"Percent weighted TPA missing: {qa_metrics['missing_weighted_tpa_percent']:.4f}%",
        "",
        "FINAL FVS TREE LIST",
        f"Final tree records: {len(tree_list):,}",
        f"Unique stands in tree file: {tree_list['STAND_ID'].nunique():,}",
        f"Direct stands: {direct_stands_count:,}",
        f"Nearest-imputed stands: {imputed_stands_count:,}",
        f"Stands still lacking trees: {len(still_missing):,}",
        "",
        "STAND AGE",
        f"Mean age: {final_stands['AGE'].mean():.2f}",
        f"Median age: {final_stands['AGE'].median():.2f}",
        f"Minimum age: {final_stands['AGE'].min():.0f}",
        f"Maximum age: {final_stands['AGE'].max():.0f}",
        "",
        "STAND-LEVEL WEIGHTED TPA",
        f"Mean stand TPA: {stand_total_tpa.mean():.2f}",
        f"Median stand TPA: {stand_total_tpa.median():.2f}",
        f"Minimum stand TPA: {stand_total_tpa.min():.2f}",
        f"Maximum stand TPA: {stand_total_tpa.max():.2f}",
        "",
        "OUTPUT FILES",
        fvs_stands_csv,
        fvs_trees_csv,
        filtered_weights_csv,
        qa_missing_species_csv,
        qa_missing_species_by_stand_csv,
        missing_csv,
    ]

    with open(
        qa_run_summary_txt,
        "w",
        encoding="utf-8",
    ) as qa_file:
        qa_file.write("\n".join(qa_lines))

    message(f"Wrote: {qa_run_summary_txt}")

    message("=" * 72)
    message("COMPLETE")
    message("=" * 72)

    message(f"Final FVS stands: {len(final_stands):,}")
    message(f"Final FVS tree records: {len(tree_list):,}")
    message(
        f"Unique stands in tree file: "
        f"{tree_list['STAND_ID'].nunique():,}"
    )
    message(
        f"Direct stands: "
        f"{tree_list.loc[tree_list['TREE_SOURCE'] == 'FIA_WEIGHTED_DIRECT', 'STAND_ID'].nunique():,}"
    )
    message(
        f"Nearest-imputed stands: "
        f"{tree_list.loc[tree_list['TREE_SOURCE'] == 'IMPUTED_NEAREST', 'STAND_ID'].nunique():,}"
    )
    message(
        f"Stands still lacking trees: {len(still_missing):,}"
    )

    if not final_stands.empty:
        message(
            "Final AGE summary:\n"
            + final_stands["AGE"]
            .describe()
            .round(2)
            .to_string()
        )

    message(
        "Final TREE_COUNT (weighted TPA) summary:\n"
        + tree_list["TREE_COUNT"]
        .describe()
        .round(4)
        .to_string()
    )


if __name__ == "__main__":
    main()