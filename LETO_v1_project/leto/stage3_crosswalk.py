"""
Build stand-level and stand-to-TreeMap donor CSVs from the final segmentation.

Outputs
-------
1. Segment_PLT_CN_Weights.csv
   One row per SEG_ID × PLT_CN combination, including the number of cells,
   donor acres, and proportional stand weight.

2. Segment_Stand_List.csv
   One row per final stand polygon with stable FVS-ready identifiers and all
   available stand-summary attributes.

Requirements
------------
- ArcGIS Pro Python environment
- arcpy
- pandas
- numpy

Important
---------
The final segment raster and final stand polygon layer must describe the same
stand boundaries and use the same SEG_ID values.
"""

import os
import math
from pathlib import Path

import arcpy
import numpy as np
import pandas as pd

from .small_stands import aggregate_management_unit_weights


# ---------------------------------------------------------------------------
# USER INPUTS
# ---------------------------------------------------------------------------

# Final management units written by Script 2.
stand_polygons = (
    r"D:\Misc\Maine_5_co\ME_Output_stands"
    r"\TreeMap_Segmentation.gdb\ManagementUnits_Final"
)

# Biological parent-segment raster written before riparian splitting.
segment_raster = (
    r"D:\Misc\Maine_5_co\ME_Output_stands"
    r"\TreeMap_Segments_Initial.tif"
)

# Spatial TreeMap raster clipped in Script 1. This must align exactly with
# segment_raster and is used for cell-by-cell VALUE pairing.
treemap_spatial_raster = (
    r"D:\Misc\Maine_5_co\ME_Output_stands"
    r"\Segmentation_Working.gdb\TreeMap_Clipped"
)

# Full CONUS TreeMap raster containing VALUE and PLT_CN in its raster
# attribute table. This raster is used only for VALUE-to-PLT_CN lookup.
treemap_lookup_raster = (
    r"D:\Misc\Whole_treemap\RDS-2025-0032 (1)"
    r"\Data\TreeMap2022_CONUS.tif"
)

output_folder = (
    r"D:\Misc\Maine_5_co\ME_Output_stands"
)

weights_csv = os.path.join(
    output_folder,
    "Segment_PLT_CN_Weights.csv",
)

stand_list_csv = os.path.join(
    output_folder,
    "Segment_Stand_List.csv",
)

# FVS metadata carried forward.
inventory_year = 2022
variant = "NE"
state = "ME"


# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

def message(text):
    print(text)
    try:
        arcpy.AddMessage(text)
    except Exception:
        pass


def require_exists(path, label):
    if not arcpy.Exists(path):
        raise FileNotFoundError(
            f"{label} does not exist or is not supported:\n{path}"
        )


def raster_properties(path):
    raster = arcpy.Raster(path)

    return {
        "width": int(raster.width),
        "height": int(raster.height),
        "cell_width": float(raster.meanCellWidth),
        "cell_height": float(raster.meanCellHeight),
        "xmin": float(raster.extent.XMin),
        "ymin": float(raster.extent.YMin),
        "xmax": float(raster.extent.XMax),
        "ymax": float(raster.extent.YMax),
        "spatial_reference": raster.spatialReference,
    }


def nearly_equal(a, b, tolerance=1e-6):
    return abs(a - b) <= tolerance


def validate_raster_alignment(segment_path, treemap_path):
    seg = raster_properties(segment_path)
    tree = raster_properties(treemap_path)

    checks = {
        "width": seg["width"] == tree["width"],
        "height": seg["height"] == tree["height"],
        "cell_width": nearly_equal(
            seg["cell_width"],
            tree["cell_width"],
        ),
        "cell_height": nearly_equal(
            seg["cell_height"],
            tree["cell_height"],
        ),
        "xmin": nearly_equal(seg["xmin"], tree["xmin"]),
        "ymin": nearly_equal(seg["ymin"], tree["ymin"]),
        "xmax": nearly_equal(seg["xmax"], tree["xmax"]),
        "ymax": nearly_equal(seg["ymax"], tree["ymax"]),
    }

    failed = [name for name, passed in checks.items() if not passed]

    if failed:
        details = "\n".join(
            [
                f"Segment raster: {seg}",
                f"Clipped spatial TreeMap raster: {tree}",
            ]
        )
        raise RuntimeError(
            "The segment raster and clipped spatial TreeMap raster are not perfectly "
            f"aligned. Failed checks: {failed}\n{details}"
        )

    message("Raster alignment check passed.")

    return seg


def get_field_lookup(dataset):
    """
    Return a case-insensitive dictionary mapping uppercase field names to
    the actual stored field names.
    """
    return {
        field.name.upper(): field.name
        for field in arcpy.ListFields(dataset)
    }


def find_field(dataset, candidates, required=True):
    lookup = get_field_lookup(dataset)

    for candidate in candidates:
        actual = lookup.get(candidate.upper())
        if actual:
            return actual

    if required:
        raise ValueError(
            f"None of these fields were found in {dataset}: {candidates}"
        )

    return None


def read_treemap_value_lookup(
    lookup_raster,
    required_values=None,
):
    """
    Read VALUE-to-PLT_CN records from the full TreeMap raster.

    When required_values is supplied, only VALUE codes occurring in the
    clipped study-area TreeMap raster are retained.
    """
    value_field = find_field(
        lookup_raster,
        ["VALUE"],
        required=True,
    )

    plt_cn_field = find_field(
        lookup_raster,
        ["PLT_CN", "PLTCN", "PLOT_CN"],
        required=True,
    )

    if required_values is not None:
        required_values = {
            int(value)
            for value in required_values
        }

    lookup = {}

    message(
        "Reading VALUE-to-PLT_CN records from the full TreeMap "
        "raster attribute table..."
    )

    with arcpy.da.SearchCursor(
        lookup_raster,
        [value_field, plt_cn_field],
    ) as cursor:
        for raster_value, plt_cn in cursor:
            if raster_value is None:
                continue

            raster_value = int(raster_value)

            if (
                required_values is not None
                and raster_value not in required_values
            ):
                continue

            if plt_cn is None:
                continue

            plt_cn_text = str(plt_cn).strip()

            if plt_cn_text.endswith(".0"):
                plt_cn_text = plt_cn_text[:-2]

            lookup[raster_value] = plt_cn_text

    if not lookup:
        raise RuntimeError(
            "No VALUE-to-PLT_CN records were read from the full "
            "TreeMap raster attribute table."
        )

    if required_values is not None:
        missing_values = sorted(
            required_values - set(lookup)
        )

        if missing_values:
            raise RuntimeError(
                f"{len(missing_values):,} TreeMap VALUE codes in the "
                "clipped study raster do not have PLT_CN records in "
                "the full TreeMap raster attribute table. "
                f"Examples: {missing_values[:20]}"
            )

    message(
        f"TreeMap VALUE-to-PLT_CN lookup records retained: "
        f"{len(lookup):,}"
    )

    return lookup


def read_required_treemap_values(treemap_path):
    """Return VALUE codes present in the clipped spatial TreeMap raster."""
    nodata = -999999999

    message(
        "Reading clipped TreeMap VALUE codes for lookup filtering..."
    )

    array = arcpy.RasterToNumPyArray(
        treemap_path,
        nodata_to_value=nodata,
    )

    values = {
        int(value)
        for value in np.unique(array[array != nodata])
    }

    if not values:
        raise RuntimeError(
            "No valid TreeMap VALUE codes were found in the clipped "
            "spatial raster."
        )

    message(
        f"Unique TreeMap VALUE codes in study area: {len(values):,}"
    )

    return values

def cell_area_acres(raster_info):
    sr = raster_info["spatial_reference"]

    if sr is None or sr.type != "Projected":
        raise RuntimeError(
            "The segment raster must use a projected coordinate system so "
            "cell acreage can be calculated reliably."
        )

    meters_per_unit = float(sr.metersPerUnit)

    cell_area_square_meters = (
        abs(raster_info["cell_width"])
        * abs(raster_info["cell_height"])
        * meters_per_unit
        * meters_per_unit
    )

    return cell_area_square_meters / 4046.8564224


def build_weights_table(
    segment_path,
    treemap_spatial_path,
    value_to_plt_cn,
    cell_acres,
):
    """
    Pair final SEG_ID cells with TreeMap raster VALUE cells and aggregate
    each SEG_ID × PLT_CN combination.
    """
    segment_nodata = -999999999
    treemap_nodata = -999999999

    message("Reading final segment raster...")
    segment_array = arcpy.RasterToNumPyArray(
        segment_path,
        nodata_to_value=segment_nodata,
    )

    message("Reading TreeMap raster...")
    treemap_array = arcpy.RasterToNumPyArray(
        treemap_spatial_path,
        nodata_to_value=treemap_nodata,
    )

    if segment_array.shape != treemap_array.shape:
        raise RuntimeError(
            "Raster arrays have different shapes after reading: "
            f"{segment_array.shape} vs {treemap_array.shape}"
        )

    valid = (
        (segment_array != segment_nodata)
        & (treemap_array != treemap_nodata)
        & (segment_array > 0)
    )

    segment_values = segment_array[valid].astype(np.int64)
    treemap_values = treemap_array[valid].astype(np.int64)

    message(
        f"Valid paired raster cells: {len(segment_values):,}"
    )

    pairs = pd.DataFrame(
        {
            "SEG_ID": segment_values,
            "TREEMAP_VALUE": treemap_values,
        }
    )

    counts = (
        pairs.groupby(
            ["SEG_ID", "TREEMAP_VALUE"],
            as_index=False,
            sort=True,
        )
        .size()
        .rename(columns={"size": "PIXEL_COUNT"})
    )

    counts["PLT_CN"] = counts["TREEMAP_VALUE"].map(
        value_to_plt_cn
    )

    missing_lookup = counts["PLT_CN"].isna()

    if missing_lookup.any():
        missing_values = sorted(
            counts.loc[
                missing_lookup,
                "TREEMAP_VALUE",
            ].unique().tolist()
        )

        preview = missing_values[:20]

        raise RuntimeError(
            "Some TreeMap raster values did not have a PLT_CN in the "
            f"raster attribute table. Example values: {preview}"
        )

    totals = (
        counts.groupby("SEG_ID")["PIXEL_COUNT"]
        .sum()
        .rename("STAND_PIXEL_COUNT")
        .reset_index()
    )

    counts = counts.merge(
        totals,
        on="SEG_ID",
        how="left",
        validate="many_to_one",
    )

    counts["DONOR_ACRES"] = (
        counts["PIXEL_COUNT"] * cell_acres
    )

    counts["WEIGHT"] = (
        counts["PIXEL_COUNT"]
        / counts["STAND_PIXEL_COUNT"]
    )

    counts["STAND_ID"] = (
        "MU_" + counts["SEG_ID"].astype(str)
    )

    counts["FVS_Stand_ID"] = counts["STAND_ID"]

    ordered = [
        "FVS_Stand_ID",
        "STAND_ID",
        "SEG_ID",
        "PLT_CN",
        "TREEMAP_VALUE",
        "PIXEL_COUNT",
        "STAND_PIXEL_COUNT",
        "DONOR_ACRES",
        "WEIGHT",
    ]

    counts = counts[ordered].sort_values(
        ["SEG_ID", "WEIGHT"],
        ascending=[True, False],
    )

    # Diagnostic: weights should sum to 1 for every stand.
    weight_sums = counts.groupby("SEG_ID")["WEIGHT"].sum()
    maximum_error = float((weight_sums - 1.0).abs().max())

    message(
        f"Maximum stand-weight sum error: {maximum_error:.12f}"
    )

    if maximum_error > 1e-9:
        raise RuntimeError(
            "Stand donor weights do not sum to 1 within tolerance."
        )

    return counts


def read_stand_polygons(polygon_path, return_components=False):
    """Read polygon components and aggregate one FVS stand per final MU_ID."""
    parent_field = find_field(polygon_path, ["SEG_ID", "GRIDCODE", "VALUE"], True)
    mu_field = find_field(polygon_path, ["MU_ID", "FINAL_ID"], True)
    excluded_types = {"OID", "Geometry", "Blob", "Raster", "GUID", "GlobalID"}
    attribute_fields = [
        field.name
        for field in arcpy.ListFields(polygon_path)
        if field.type not in excluded_types
        and field.name.upper() not in {"SHAPE_LENGTH", "SHAPE_AREA"}
    ]
    rows = []
    with arcpy.da.SearchCursor(polygon_path, attribute_fields) as cursor:
        for row in cursor:
            rows.append(dict(zip(attribute_fields, row)))
    components = pd.DataFrame(rows)
    if components.empty:
        raise RuntimeError("No final stand polygon records were read.")

    components["SEG_ID"] = pd.to_numeric(
        components[parent_field], errors="raise"
    ).astype("int64")
    components["MU_ID"] = pd.to_numeric(
        components[mu_field], errors="raise"
    ).astype("int64")
    components["Acres"] = pd.to_numeric(components["Acres"], errors="raise")
    components["ParentSEG_ID"] = components["SEG_ID"]
    if (components["Acres"] <= 0).any():
        raise RuntimeError("Final management-unit components must have positive acreage.")

    stand_rows = []
    for mu_id, group in components.groupby("MU_ID", sort=True):
        row = group.iloc[0].to_dict()
        acres = group["Acres"].sum()
        row["MU_ID"] = int(mu_id)
        row["SEG_ID"] = int(mu_id)
        row["ParentSEG_ID"] = int(group["SEG_ID"].min())
        row["SEG_IDS"] = ",".join(
            str(value) for value in sorted(group["SEG_ID"].astype(int).unique())
        )
        row["SEGMENT_COUNT"] = int(group["SEG_ID"].nunique())
        row["Acres"] = float(acres)
        if "STDAGE_MEAN" in group.columns:
            ages = pd.to_numeric(group["STDAGE_MEAN"], errors="coerce")
            valid_age = ages.notna()
            if valid_age.any():
                row["STDAGE_MEAN"] = float(
                    (ages[valid_age] * group.loc[valid_age, "Acres"]).sum()
                    / group.loc[valid_age, "Acres"].sum()
                )
        if "MGMT_TYPE" in group.columns:
            management_types = set(group["MGMT_TYPE"].dropna().astype(str))
            if "Riparian" in management_types:
                row["MGMT_TYPE"] = "Riparian"
                row["DB_GROUP"] = "Riparian"
            elif "Wetland" in management_types:
                row["MGMT_TYPE"] = "Wetland"
        stand_rows.append(row)

    stands = pd.DataFrame(stand_rows)
    stands["MU_ID"] = stands["MU_ID"].astype("int64")
    stands["SEG_ID"] = stands["MU_ID"]
    stands["STAND_ID"] = "MU_" + stands["MU_ID"].astype(str)
    stands["FVS_Stand_ID"] = stands["STAND_ID"]
    stands["Inventory_Year"] = inventory_year
    stands["Variant"] = variant
    stands["State"] = state
    if "STDAGE_MEAN" in stands.columns:
        stands["AGE"] = pd.to_numeric(
            stands["STDAGE_MEAN"], errors="coerce"
        ).round().astype("Int64")
    front = [
        "FVS_Stand_ID", "STAND_ID", "SEG_ID", "MU_ID", "ParentSEG_ID",
        "SEG_IDS", "SEGMENT_COUNT", "Inventory_Year", "Variant", "State",
    ]
    if "AGE" in stands.columns:
        front.append("AGE")
    remaining = [
        column for column in stands.columns
        if column not in front and column not in {parent_field, mu_field}
    ]
    stands = stands[front + remaining].sort_values("MU_ID").reset_index(drop=True)
    message(
        f"Final management-unit polygon records: {len(components):,}; "
        f"unique MU_ID values: {len(stands):,}"
    )
    if return_components:
        return stands, components
    return stands


def expand_parent_weights_to_management_units(
    parent_weights,
    stands,
    components=None,
):
    """Acreage-weight parent donors across every SEG_ID contributing to MU_ID."""
    if components is None:
        components = stands.rename(columns={"ParentSEG_ID": "SEG_ID"}).copy()
        if "MU_ID" not in components.columns:
            components["MU_ID"] = components["SEG_ID"]
    expanded = aggregate_management_unit_weights(parent_weights, components)
    stand_mapping = stands[["MU_ID", "STAND_ID", "FVS_Stand_ID"]].copy()
    stand_mapping["MU_ID"] = pd.to_numeric(
        stand_mapping["MU_ID"], errors="raise"
    ).astype("int64")
    expanded = expanded.drop(
        columns=["STAND_ID", "FVS_Stand_ID"], errors="ignore"
    ).merge(
        stand_mapping, on="MU_ID", how="left", validate="many_to_one"
    )
    if expanded["STAND_ID"].isna().any():
        raise RuntimeError("Donor weights reference an unknown MU_ID.")
    ordered = [
        "FVS_Stand_ID", "STAND_ID", "SEG_ID", "MU_ID", "PLT_CN",
        "SOURCE_SEG_IDS", "DONOR_ACRES", "MU_ACRES", "WEIGHT",
        "WEIGHT_SUM_ERROR",
    ]
    return expanded[
        [column for column in ordered if column in expanded.columns]
    ].sort_values(
        ["MU_ID", "WEIGHT"], ascending=[True, False]
    ).reset_index(drop=True)


def validate_segment_ids(parent_weights, stands, components=None):
    raster_ids = set(parent_weights["SEG_ID"].unique().tolist())
    if components is None:
        parent_ids = set(stands["ParentSEG_ID"].unique().tolist())
    else:
        parent_ids = set(components["SEG_ID"].unique().tolist())
    missing = sorted(parent_ids - raster_ids)
    if missing:
        raise RuntimeError(
            "Some final management units reference parent segments absent "
            f"from the raster weights. Examples: {missing[:20]}"
        )
    message(
        f"Parent-segment donor coverage passed for {len(parent_ids):,} "
        "parent stands."
    )


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    arcpy.env.overwriteOutput = True

    os.makedirs(output_folder, exist_ok=True)

    require_exists(
        stand_polygons,
        "Final stand polygon feature class",
    )
    require_exists(
        segment_raster,
        "Final segment raster",
    )
    require_exists(
        treemap_spatial_raster,
        "Clipped spatial TreeMap raster",
    )
    require_exists(
        treemap_lookup_raster,
        "Full TreeMap lookup raster",
    )

    message("=" * 70)
    message("VALIDATING INPUTS")
    message("=" * 70)

    raster_info = validate_raster_alignment(
        segment_raster,
        treemap_spatial_raster,
    )

    acres_per_cell = cell_area_acres(raster_info)

    message(
        f"Raster cell area: {acres_per_cell:.6f} acres"
    )

    required_values = read_required_treemap_values(
        treemap_spatial_raster
    )

    value_to_plt_cn = read_treemap_value_lookup(
        treemap_lookup_raster,
        required_values=required_values,
    )

    message("=" * 70)
    message("BUILDING SEGMENT-TO-PLT_CN WEIGHTS")
    message("=" * 70)

    parent_weights = build_weights_table(
        segment_raster,
        treemap_spatial_raster,
        value_to_plt_cn,
        acres_per_cell,
    )

    message("=" * 70)
    message("READING FINAL STAND POLYGONS")
    message("=" * 70)

    stands = read_stand_polygons(
        stand_polygons
    )

    validate_segment_ids(
        parent_weights,
        stands,
    )

    weights = expand_parent_weights_to_management_units(
        parent_weights,
        stands,
    )

    message("=" * 70)
    message("WRITING OUTPUTS")
    message("=" * 70)

    weights.to_csv(
        weights_csv,
        index=False,
    )

    stands.to_csv(
        stand_list_csv,
        index=False,
    )

    message(f"Wrote: {weights_csv}")
    message(f"Wrote: {stand_list_csv}")

    message("=" * 70)
    message("COMPLETE")
    message("=" * 70)
    message(
        f"Stands exported: {len(stands):,}"
    )
    message(
        f"Stand-to-donor combinations: {len(weights):,}"
    )
    message(
        f"Unique PLT_CN donors used: "
        f"{weights['PLT_CN'].nunique():,}"
    )


if __name__ == "__main__":
    main()