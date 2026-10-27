import os
import arcpy
import numpy as np
import pandas as pd
from scipy import ndimage


# ============================================================
# USER INPUTS
# ============================================================

# TreeMap raster containing VALUE and PLT_CN in its raster
# attribute table.
treemap_raster = r"C:\path\to\data.gdb\TreeMap_2022_Clip"

fia_tree_files = [
    r"C:\path\to\FL_TREE.csv",
    r"C:\path\to\GA_TREE.csv",
    r"C:\path\to\AL_TREE.csv",
    r"C:\path\to\SC_TREE.csv",
]

fia_cond_files = [
    r"C:\path\to\FL_COND.csv",
    r"C:\path\to\GA_COND.csv",
    r"C:\path\to\AL_COND.csv",
    r"C:\path\to\SC_COND.csv",
]

output_folder = (
    r"C:\path\to\data"
)

# Field in the TreeMap raster attribute table.
treemap_plt_field = "PLT_CN"

# Optional clip boundary.
#
# Leave as None to process the full TreeMap raster.
# For initial testing, using a smaller study boundary is strongly
# recommended.
clip_boundary = None

# Example:
# clip_boundary = (
#     r"C:\path\to\data.gdb"
#     r"\Study_Area"
# )


# ============================================================
# OUTPUTS
# ============================================================

working_gdb = os.path.join(
    output_folder,
    "Segmentation_Working.gdb",
)

lookup_csv = os.path.join(
    output_folder,
    "TreeMap_PLT_CN_Segmentation_Lookup.csv",
)

age_raster = os.path.join(
    output_folder,
    "TreeMap_STDAGE.tif",
)

forest_type_raster = os.path.join(
    output_folder,
    "TreeMap_FORTYPCD.tif",
)

basal_area_raster = os.path.join(
    output_folder,
    "TreeMap_BALIVE.tif",
)

qmd_raster = os.path.join(
    output_folder,
    "TreeMap_QMD.tif",
)

tpa_raster = os.path.join(
    output_folder,
    "TreeMap_TPA.tif",
)


# ============================================================
# SETTINGS
# ============================================================

arcpy.env.overwriteOutput = True

os.makedirs(
    output_folder,
    exist_ok=True,
)

if not arcpy.Exists(working_gdb):
    arcpy.management.CreateFileGDB(
        output_folder,
        os.path.basename(working_gdb),
    )


# ============================================================
# HELPERS
# ============================================================

def message(text):
    print(text)

    try:
        arcpy.AddMessage(text)
    except Exception:
        pass


def normalize_cn(series):
    """
    Normalize FIA control numbers without converting them to
    floating-point values.
    """

    return (
        series
        .astype("string")
        .str.strip()
        .str.replace(r"\.0$", "", regex=True)
    )


def read_condition_files(file_paths):
    """
    Read only the FIA condition fields needed for raster creation.
    """

    fields = [
        "PLT_CN",
        "CONDID",
        "COND_STATUS_CD",
        "CONDPROP_UNADJ",
        "STDAGE",
        "FORTYPCD",
        "BALIVE",
    ]

    frames = []

    for file_path in file_paths:
        message(
            f"Reading condition file: "
            f"{os.path.basename(file_path)}"
        )

        frame = pd.read_csv(
            file_path,
            usecols=fields,
            dtype={
                "PLT_CN": "string",
            },
            low_memory=False,
        )

        frames.append(frame)

    conditions = pd.concat(
        frames,
        ignore_index=True,
    )

    conditions["PLT_CN"] = normalize_cn(
        conditions["PLT_CN"]
    )

    numeric_fields = [
        "CONDID",
        "COND_STATUS_CD",
        "CONDPROP_UNADJ",
        "STDAGE",
        "FORTYPCD",
        "BALIVE",
    ]

    for field in numeric_fields:
        conditions[field] = pd.to_numeric(
            conditions[field],
            errors="coerce",
        )

    return conditions


def select_dominant_conditions(conditions):
    """
    Select the sampled forest condition having the largest
    CONDPROP_UNADJ for each PLT_CN.
    """

    forest = conditions[
        conditions["COND_STATUS_CD"] == 1
    ].copy()

    forest = forest[
        forest["PLT_CN"].notna()
        & forest["CONDID"].notna()
    ].copy()

    # Sort so the greatest condition proportion is first.
    #
    # CONDID is included as a deterministic tie-breaker.
    forest = forest.sort_values(
        by=[
            "PLT_CN",
            "CONDPROP_UNADJ",
            "CONDID",
        ],
        ascending=[
            True,
            False,
            True,
        ],
        na_position="last",
    )

    dominant = (
        forest
        .drop_duplicates(
            subset="PLT_CN",
            keep="first",
        )
        .copy()
    )

    dominant["CONDID"] = (
        dominant["CONDID"]
        .round()
        .astype("Int64")
    )

    return dominant


def read_and_summarize_trees(
    file_paths,
    dominant_conditions,
):
    """
    Calculate live-tree TPA and QMD for the dominant condition
    of each PLT_CN.

    QMD =
        sqrt(
            sum(TPA_UNADJ * DIA^2)
            /
            sum(TPA_UNADJ)
        )
    """

    dominant_keys = (
        dominant_conditions[
            ["PLT_CN", "CONDID"]
        ]
        .dropna()
        .copy()
    )

    dominant_keys["CONDID"] = (
        dominant_keys["CONDID"]
        .astype("int64")
    )

    frames = []

    tree_fields = [
        "PLT_CN",
        "CONDID",
        "STATUSCD",
        "DIA",
        "TPA_UNADJ",
    ]

    for file_path in file_paths:
        message(
            f"Reading tree file: "
            f"{os.path.basename(file_path)}"
        )

        trees = pd.read_csv(
            file_path,
            usecols=tree_fields,
            dtype={
                "PLT_CN": "string",
            },
            low_memory=False,
        )

        trees["PLT_CN"] = normalize_cn(
            trees["PLT_CN"]
        )

        for field in [
            "CONDID",
            "STATUSCD",
            "DIA",
            "TPA_UNADJ",
        ]:
            trees[field] = pd.to_numeric(
                trees[field],
                errors="coerce",
            )

        # Keep live trees with usable diameter and expansion
        # factors.
        trees = trees[
            (trees["STATUSCD"] == 1)
            & trees["PLT_CN"].notna()
            & trees["CONDID"].notna()
            & trees["DIA"].notna()
            & trees["TPA_UNADJ"].notna()
            & (trees["DIA"] > 0)
            & (trees["TPA_UNADJ"] > 0)
        ].copy()

        trees["CONDID"] = (
            trees["CONDID"]
            .round()
            .astype("int64")
        )

        # Retain only trees belonging to each plot's selected
        # dominant condition.
        trees = trees.merge(
            dominant_keys,
            on=[
                "PLT_CN",
                "CONDID",
            ],
            how="inner",
            validate="many_to_one",
        )

        if trees.empty:
            continue

        trees["TPA_DIA2"] = (
            trees["TPA_UNADJ"]
            * np.square(trees["DIA"])
        )

        summary = (
            trees
            .groupby(
                ["PLT_CN", "CONDID"],
                as_index=False,
            )
            .agg(
                TPA=(
                    "TPA_UNADJ",
                    "sum",
                ),
                TPA_DIA2=(
                    "TPA_DIA2",
                    "sum",
                ),
                LIVE_TREE_RECORDS=(
                    "DIA",
                    "size",
                ),
            )
        )

        frames.append(summary)

    if not frames:
        return pd.DataFrame(
            columns=[
                "PLT_CN",
                "CONDID",
                "TPA",
                "QMD",
                "LIVE_TREE_RECORDS",
            ]
        )

    # A PLT_CN should normally occur in only one supplied state
    # file. This extra aggregation safely handles duplicates.
    summaries = pd.concat(
        frames,
        ignore_index=True,
    )

    summaries = (
        summaries
        .groupby(
            ["PLT_CN", "CONDID"],
            as_index=False,
        )
        .agg(
            TPA=("TPA", "sum"),
            TPA_DIA2=("TPA_DIA2", "sum"),
            LIVE_TREE_RECORDS=(
                "LIVE_TREE_RECORDS",
                "sum",
            ),
        )
    )

    summaries["QMD"] = np.sqrt(
        summaries["TPA_DIA2"]
        / summaries["TPA"]
    )

    return summaries[
        [
            "PLT_CN",
            "CONDID",
            "TPA",
            "QMD",
            "LIVE_TREE_RECORDS",
        ]
    ]


def prepare_source_raster():
    """
    Optionally clip TreeMap to a test boundary while retaining
    the raster's VALUE field and attribute table.
    """

    if clip_boundary is None:
        return treemap_raster

    clipped = os.path.join(
        working_gdb,
        "TreeMap_Clipped",
    )

    message("Clipping TreeMap to the test boundary...")

    arcpy.management.Clip(
        in_raster=treemap_raster,
        rectangle="#",
        out_raster=clipped,
        in_template_dataset=clip_boundary,
        nodata_value="#",
        clipping_geometry="ClippingGeometry",
        maintain_clipping_extent="NO_MAINTAIN_EXTENT",
    )

    try:
        arcpy.management.BuildRasterAttributeTable(
            clipped,
            "Overwrite",
        )
    except arcpy.ExecuteError:
        message(
            "The clipped raster attribute table could not be "
            "rebuilt. Existing attributes will be used."
        )

    return clipped


def read_raster_value_lookup(
    source_raster,
    lookup_table,
):
    """
    Create a record for every VALUE in the TreeMap raster attribute table.

    FIA-derived values are preferred when available. Existing TreeMap raster
    attributes are used as fallbacks so valid TreeMap VALUE codes are not
    dropped merely because a PLT_CN is absent from the supplied FIA CSV files
    or one FIA-derived field is null.
    """

    fia_records = (
        lookup_table
        .dropna(subset=["PLT_CN"])
        .copy()
    )

    fia_records["PLT_CN"] = normalize_cn(
        fia_records["PLT_CN"]
    )

    fia_records = (
        fia_records
        .drop_duplicates(
            subset="PLT_CN",
            keep="first",
        )
        .set_index("PLT_CN")
        .to_dict(orient="index")
    )

    available_fields = {
        field.name.upper(): field.name
        for field in arcpy.ListFields(source_raster)
    }

    # TreeMap field names used as fallbacks.
    fallback_candidates = {
        "STDAGE": ["STDAGE"],
        "FORTYPCD": ["FORTYPCD"],
        "BALIVE": ["BALIVE"],
        "QMD": ["QMD"],
        "TPA": ["TPA_LIVE", "TPA"],
    }

    fallback_fields = {}

    for output_field, candidates in fallback_candidates.items():
        for candidate in candidates:
            actual = available_fields.get(candidate.upper())

            if actual is not None:
                fallback_fields[output_field] = actual
                break

    cursor_fields = [
        available_fields.get("VALUE", "VALUE"),
        treemap_plt_field,
    ]

    for actual_field in fallback_fields.values():
        if actual_field not in cursor_fields:
            cursor_fields.append(actual_field)

    field_positions = {
        field_name: index
        for index, field_name in enumerate(cursor_fields)
    }

    value_records = {}
    fia_matches = 0
    rat_fallback_values = 0
    unmatched_plots = set()

    message(
        "Reading every VALUE and PLT_CN from the TreeMap raster "
        "attribute table..."
    )

    with arcpy.da.SearchCursor(
        source_raster,
        cursor_fields,
    ) as cursor:

        for row in cursor:
            raster_value = row[0]
            plot_cn = row[1]

            if raster_value is None:
                continue

            normalized_plot = None

            if plot_cn is not None:
                normalized_plot = str(plot_cn).strip()

                if normalized_plot.endswith(".0"):
                    normalized_plot = normalized_plot[:-2]

            record = {}

            if (
                normalized_plot is not None
                and normalized_plot in fia_records
            ):
                record.update(
                    fia_records[normalized_plot]
                )
                fia_matches += 1
            elif normalized_plot is not None:
                unmatched_plots.add(normalized_plot)

            # Fill absent FIA fields from the existing TreeMap RAT.
            for output_field, actual_field in fallback_fields.items():
                current_value = record.get(output_field)

                if pd.isna(current_value):
                    fallback_value = row[
                        field_positions[actual_field]
                    ]

                    if fallback_value is not None:
                        record[output_field] = fallback_value
                        rat_fallback_values += 1

            record["PLT_CN"] = normalized_plot
            value_records[int(raster_value)] = record

    message(
        f"TreeMap VALUE records retained: "
        f"{len(value_records):,}"
    )
    message(
        f"VALUE records matched to FIA lookup: "
        f"{fia_matches:,}"
    )
    message(
        f"TreeMap RAT fallback values used: "
        f"{rat_fallback_values:,}"
    )
    message(
        f"Unique PLT_CN absent from FIA lookup: "
        f"{len(unmatched_plots):,}"
    )

    return value_records

def write_lookup_raster(
    source_array,
    value_records,
    field,
    output_path,
    lower_left,
    cell_width,
    cell_height,
    spatial_reference,
    output_dtype,
    output_nodata,
    decimal_multiplier=1.0,
):
    """
    Write one segmentation raster while preserving the full valid TreeMap
    footprint.

    Values are assigned by TreeMap VALUE. Any remaining missing cells inside
    the valid TreeMap footprint are filled from the nearest assigned cell so
    all segmentation rasters have identical coverage.
    """

    message(
        f"Building {field} raster..."
    )

    source_valid = source_array != source_nodata

    working_array = np.full(
        source_array.shape,
        np.nan,
        dtype=np.float64,
    )

    matched_values = 0

    for raster_value, record in value_records.items():
        field_value = record.get(field)

        if pd.isna(field_value):
            continue

        stored_value = (
            float(field_value)
            * decimal_multiplier
        )

        cell_mask = source_array == raster_value

        if np.any(cell_mask):
            working_array[cell_mask] = stored_value
            matched_values += 1

    assigned_before_fill = (
        source_valid
        & np.isfinite(working_array)
    )

    missing_mask = (
        source_valid
        & ~np.isfinite(working_array)
    )

    missing_before_fill = int(
        missing_mask.sum()
    )

    if missing_before_fill > 0:
        donor_mask = assigned_before_fill

        if not np.any(donor_mask):
            raise RuntimeError(
                f"{field} has no assigned donor cells."
            )

        nearest_indices = ndimage.distance_transform_edt(
            ~donor_mask,
            return_distances=False,
            return_indices=True,
        )

        nearest_values = working_array[
            tuple(nearest_indices)
        ]

        working_array[missing_mask] = (
            nearest_values[missing_mask]
        )

    missing_after_fill = int(
        np.count_nonzero(
            source_valid
            & ~np.isfinite(working_array)
        )
    )

    if missing_after_fill > 0:
        raise RuntimeError(
            f"{field} still has {missing_after_fill:,} missing "
            "valid TreeMap cells after nearest-neighbor filling."
        )

    rounded = np.rint(
        working_array
    )

    output_array = np.full(
        source_array.shape,
        output_nodata,
        dtype=output_dtype,
    )

    output_array[source_valid] = rounded[
        source_valid
    ].astype(output_dtype)

    matched_cells = int(
        np.count_nonzero(
            source_valid
        )
    )

    message(
        f"  TreeMap VALUE codes with direct values: "
        f"{matched_values:,}"
    )
    message(
        f"  Cells filled from nearest valid value: "
        f"{missing_before_fill:,}"
    )
    message(
        f"  Valid TreeMap cells assigned: "
        f"{matched_cells:,}"
    )

    output_raster = arcpy.NumPyArrayToRaster(
        output_array,
        lower_left,
        cell_width,
        cell_height,
        value_to_nodata=output_nodata,
    )

    output_raster.save(output_path)

    arcpy.management.DefineProjection(
        output_path,
        spatial_reference,
    )

    return matched_cells


# ============================================================
# READ FIA CONDITIONS
# ============================================================

message("=" * 60)
message("READING FIA CONDITION DATA")
message("=" * 60)

conditions = read_condition_files(
    fia_cond_files
)

message(
    f"Total condition records read: "
    f"{len(conditions):,}"
)

dominant_conditions = select_dominant_conditions(
    conditions
)

message(
    f"Unique plots with a dominant forest condition: "
    f"{dominant_conditions['PLT_CN'].nunique():,}"
)


# ============================================================
# CALCULATE DOMINANT-CONDITION TREE STATISTICS
# ============================================================

message("=" * 60)
message("CALCULATING TPA AND QMD")
message("=" * 60)

tree_summary = read_and_summarize_trees(
    fia_tree_files,
    dominant_conditions,
)

message(
    f"Dominant conditions with live-tree summaries: "
    f"{len(tree_summary):,}"
)


# ============================================================
# BUILD FINAL PLT_CN LOOKUP
# ============================================================

lookup = dominant_conditions.merge(
    tree_summary,
    on=[
        "PLT_CN",
        "CONDID",
    ],
    how="left",
    validate="one_to_one",
)

lookup = lookup.rename(
    columns={
        "CONDPROP_UNADJ": "DOM_COND_PROP",
    }
)

lookup_fields = [
    "PLT_CN",
    "CONDID",
    "DOM_COND_PROP",
    "STDAGE",
    "FORTYPCD",
    "BALIVE",
    "TPA",
    "QMD",
    "LIVE_TREE_RECORDS",
]

lookup = lookup[lookup_fields].copy()

lookup.to_csv(
    lookup_csv,
    index=False,
)

message(
    f"Lookup table written:\n{lookup_csv}"
)

message(
    f"Plots with STDAGE: "
    f"{lookup['STDAGE'].notna().sum():,}"
)

message(
    f"Plots with FORTYPCD: "
    f"{lookup['FORTYPCD'].notna().sum():,}"
)

message(
    f"Plots with BALIVE: "
    f"{lookup['BALIVE'].notna().sum():,}"
)

message(
    f"Plots with QMD: "
    f"{lookup['QMD'].notna().sum():,}"
)

message(
    f"Plots with TPA: "
    f"{lookup['TPA'].notna().sum():,}"
)


# ============================================================
# PREPARE TREEMAP
# ============================================================

message("=" * 60)
message("PREPARING TREEMAP RASTER")
message("=" * 60)

source_raster = prepare_source_raster()

arcpy.env.snapRaster = source_raster
arcpy.env.cellSize = source_raster
arcpy.env.extent = source_raster
arcpy.env.outputCoordinateSystem = source_raster

description = arcpy.Describe(
    source_raster
)

lower_left = arcpy.Point(
    description.extent.XMin,
    description.extent.YMin,
)

cell_width = description.meanCellWidth
cell_height = description.meanCellHeight
spatial_reference = description.spatialReference

source_nodata = -999999

message("Reading TreeMap raster into a NumPy array...")

source_array = arcpy.RasterToNumPyArray(
    source_raster,
    nodata_to_value=source_nodata,
)

message(
    f"Raster dimensions: "
    f"{source_array.shape[0]:,} rows by "
    f"{source_array.shape[1]:,} columns"
)

value_records = read_raster_value_lookup(
    source_raster,
    lookup,
)

message(
    f"TreeMap VALUE codes linked to FIA attributes: "
    f"{len(value_records):,}"
)


# ============================================================
# CREATE SEGMENTATION RASTERS
# ============================================================

message("=" * 60)
message("CREATING SEGMENTATION RASTERS")
message("=" * 60)

# Age in whole years.
write_lookup_raster(
    source_array=source_array,
    value_records=value_records,
    field="STDAGE",
    output_path=age_raster,
    lower_left=lower_left,
    cell_width=cell_width,
    cell_height=cell_height,
    spatial_reference=spatial_reference,
    output_dtype=np.int16,
    output_nodata=-9999,
    decimal_multiplier=1.0,
)

# FIA forest type code.
write_lookup_raster(
    source_array=source_array,
    value_records=value_records,
    field="FORTYPCD",
    output_path=forest_type_raster,
    lower_left=lower_left,
    cell_width=cell_width,
    cell_height=cell_height,
    spatial_reference=spatial_reference,
    output_dtype=np.int16,
    output_nodata=-9999,
    decimal_multiplier=1.0,
)

# BALIVE is stored in tenths of square feet per acre.
#
# Raster value 287 means 28.7 ft2/ac.
write_lookup_raster(
    source_array=source_array,
    value_records=value_records,
    field="BALIVE",
    output_path=basal_area_raster,
    lower_left=lower_left,
    cell_width=cell_width,
    cell_height=cell_height,
    spatial_reference=spatial_reference,
    output_dtype=np.int32,
    output_nodata=-9999,
    decimal_multiplier=10.0,
)

# QMD is stored in tenths of an inch.
#
# Raster value 112 means 11.2 inches.
write_lookup_raster(
    source_array=source_array,
    value_records=value_records,
    field="QMD",
    output_path=qmd_raster,
    lower_left=lower_left,
    cell_width=cell_width,
    cell_height=cell_height,
    spatial_reference=spatial_reference,
    output_dtype=np.int16,
    output_nodata=-9999,
    decimal_multiplier=10.0,
)

# TPA is stored in tenths of a tree per acre.
#
# Raster value 1247 means 124.7 trees/ac.
write_lookup_raster(
    source_array=source_array,
    value_records=value_records,
    field="TPA",
    output_path=tpa_raster,
    lower_left=lower_left,
    cell_width=cell_width,
    cell_height=cell_height,
    spatial_reference=spatial_reference,
    output_dtype=np.int32,
    output_nodata=-9999,
    decimal_multiplier=10.0,
)


# ============================================================
# VERIFY IDENTICAL RASTER COVERAGE
# ============================================================

message("=" * 60)
message("VERIFYING SEGMENTATION RASTER COVERAGE")
message("=" * 60)

source_valid_count = int(
    np.count_nonzero(
        source_array != source_nodata
    )
)

for raster_path in [
    age_raster,
    forest_type_raster,
    basal_area_raster,
    qmd_raster,
    tpa_raster,
]:
    check_array = arcpy.RasterToNumPyArray(
        raster_path,
        nodata_to_value=-999999999,
    )

    assigned_count = int(
        np.count_nonzero(
            check_array != -999999999
        )
    )

    missing_count = (
        source_valid_count - assigned_count
    )

    message(
        f"{os.path.basename(raster_path)}: "
        f"assigned={assigned_count:,}; "
        f"missing={missing_count:,}"
    )

    if missing_count != 0:
        raise RuntimeError(
            f"Coverage QA failed for {raster_path}"
        )

# ============================================================
# BUILD PYRAMIDS AND STATISTICS
# ============================================================

message("=" * 60)
message("BUILDING RASTER STATISTICS")
message("=" * 60)

created_rasters = [
    age_raster,
    forest_type_raster,
    basal_area_raster,
    qmd_raster,
    tpa_raster,
]

for raster_path in created_rasters:
    message(
        f"Calculating statistics: "
        f"{os.path.basename(raster_path)}"
    )

    arcpy.management.CalculateStatistics(
        raster_path,
    )

    try:
        arcpy.management.BuildPyramids(
            raster_path,
        )
    except arcpy.ExecuteError:
        message(
            "  Pyramids were not built, but the raster "
            "was created successfully."
        )


# ============================================================
# FINISH
# ============================================================

message("=" * 60)
message("SEGMENTATION RASTERS COMPLETE")
message("=" * 60)

for raster_path in created_rasters:
    message(raster_path)

message(
    "\nStored-value scales:"
)

message(
    "  STDAGE: value = years"
)

message(
    "  FORTYPCD: value = FIA forest type code"
)

message(
    "  BALIVE: divide raster value by 10 for ft2/ac"
)

message(
    "  QMD: divide raster value by 10 for inches"
)

message(
    "  TPA: divide raster value by 10 for trees/ac"
)