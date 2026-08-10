import os
import arcpy
import numpy as np
import pandas as pd
from scipy import ndimage




# Configured by run_stage().
treemap_raster = None
treemap_lookup_raster = None
fia_tree_files = []
fia_cond_files = []
output_folder = None
treemap_plt_field = "PLT_CN"
clip_boundary = None
working_gdb = None
lookup_csv = None
age_raster = None
forest_type_raster = None
basal_area_raster = None
qmd_raster = None
tpa_raster = None
source_nodata = None
forest_condition_status_code = None
live_tree_status_code = None
output_nodata = None
coverage_check_nodata = None
stored_value_multipliers = {}
stage1_qa_summary_csv = None
stage1_qa_landscape_csv = None
stage1_qa_provenance_rasters = True

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
        conditions["COND_STATUS_CD"] == forest_condition_status_code
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
            (trees["STATUSCD"] == live_tree_status_code)
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
    Optionally clip TreeMap to the study boundary.

    The clipped raster controls the spatial footprint, cell alignment,
    and output geometry. Its extended raster attribute fields are not
    required because those are read separately from treemap_lookup_raster.
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
    lookup_raster,
    lookup_table,
    required_values=None,
):
    """
    Build records for TreeMap VALUE codes from the full TreeMap raster
    attribute table.

    Parameters
    ----------
    lookup_raster : str
        Full TreeMap raster containing VALUE, PLT_CN, and available
        fallback attributes.
    lookup_table : pandas.DataFrame
        FIA-derived PLT_CN lookup.
    required_values : set[int] or None
        Optional set of raster VALUE codes present in the clipped study
        raster. When supplied, national TreeMap records outside the study
        area are skipped.

    FIA-derived values are preferred. Existing TreeMap raster attributes
    are used only as fallbacks when a supplied FIA record or field is absent.
    """

    if not arcpy.Exists(lookup_raster):
        raise FileNotFoundError(
            "TreeMap lookup raster was not found:\n"
            f"{lookup_raster}"
        )

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
        for field in arcpy.ListFields(lookup_raster)
    }

    value_field = available_fields.get("VALUE")

    if value_field is None:
        raise RuntimeError(
            "The TreeMap lookup raster does not contain a VALUE field:\n"
            f"{lookup_raster}"
        )

    plt_field = available_fields.get(
        treemap_plt_field.upper()
    )

    if plt_field is None:
        alternatives = [
            available_fields.get("PLTCN"),
            available_fields.get("PLOT_CN"),
        ]
        plt_field = next(
            (field for field in alternatives if field is not None),
            None,
        )

    if plt_field is None:
        field_names = ", ".join(
            sorted(available_fields.values())
        )
        raise RuntimeError(
            "The TreeMap lookup raster does not contain PLT_CN "
            "or a recognized alternative.\n"
            f"Raster: {lookup_raster}\n"
            f"Available fields: {field_names}"
        )

    # TreeMap fields used as fallbacks when FIA-derived values are absent.
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

    missing_fallbacks = sorted(
        set(fallback_candidates) - set(fallback_fields)
    )

    if missing_fallbacks:
        message(
            "TreeMap fallback fields not found for: "
            + ", ".join(missing_fallbacks)
        )

    cursor_fields = [
        value_field,
        plt_field,
    ]

    for actual_field in fallback_fields.values():
        if actual_field not in cursor_fields:
            cursor_fields.append(actual_field)

    field_positions = {
        field_name: index
        for index, field_name in enumerate(cursor_fields)
    }

    if required_values is not None:
        required_values = {
            int(value)
            for value in required_values
        }

    value_records = {}
    fia_matches = 0
    rat_fallback_values = 0
    unmatched_plots = set()

    message(
        "Reading VALUE, PLT_CN, and fallback attributes from "
        "the full TreeMap raster attribute table..."
    )

    with arcpy.da.SearchCursor(
        lookup_raster,
        cursor_fields,
    ) as cursor:

        for row in cursor:
            raster_value = row[0]
            plot_cn = row[1]

            if raster_value is None:
                continue

            raster_value = int(raster_value)

            if (
                required_values is not None
                and raster_value not in required_values
            ):
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
                for source_field in fallback_candidates:
                    if not pd.isna(record.get(source_field)):
                        record[f"_SOURCE_{source_field}"] = 1
            elif normalized_plot is not None:
                unmatched_plots.add(normalized_plot)

            # Fill missing FIA-derived fields from the original TreeMap RAT.
            for output_field, actual_field in fallback_fields.items():
                current_value = record.get(output_field)

                if pd.isna(current_value):
                    fallback_value = row[
                        field_positions[actual_field]
                    ]

                    if fallback_value is not None:
                        record[output_field] = fallback_value
                        record[f"_SOURCE_{output_field}"] = 2
                        rat_fallback_values += 1

            record["_FIA_PLOT_FOUND"] = bool(
                normalized_plot is not None and normalized_plot in fia_records
            )
            for provenance_field in fallback_candidates:
                source_code = record.get(f"_SOURCE_{provenance_field}")
                record[f"_FIA_ATTRIBUTE_AVAILABLE_{provenance_field}"] = (
                    source_code == 1
                )
                record[f"_RAT_ATTRIBUTE_AVAILABLE_{provenance_field}"] = (
                    source_code == 2
                )
            record["PLT_CN"] = normalized_plot
            value_records[raster_value] = record

    if required_values is not None:
        missing_values = sorted(
            required_values - set(value_records)
        )
    else:
        missing_values = []

    message(
        f"TreeMap VALUE records retained: "
        f"{len(value_records):,}"
    )
    message(
        f"VALUE records matched to supplied FIA lookup: "
        f"{fia_matches:,}"
    )
    message(
        f"TreeMap RAT fallback values used: "
        f"{rat_fallback_values:,}"
    )
    message(
        f"Unique PLT_CN absent from supplied FIA lookup: "
        f"{len(unmatched_plots):,}"
    )

    if missing_values:
        preview = missing_values[:20]
        raise RuntimeError(
            f"{len(missing_values):,} VALUE codes present in the "
            "clipped TreeMap raster were not found in the full "
            "TreeMap attribute table.\n"
            f"Example missing VALUE codes: {preview}"
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
    provenance = np.zeros(source_array.shape, dtype=np.uint8)
    reason_masks = {
        "MISSING_FIA_PLT_CN": np.zeros(source_array.shape, dtype=bool),
        "FIA_ATTRIBUTE_MISSING": np.zeros(source_array.shape, dtype=bool),
        "NO_LIVE_TREES": np.zeros(source_array.shape, dtype=bool),
        "TREEMAP_RAT_ATTRIBUTE_MISSING": np.zeros(source_array.shape, dtype=bool),
    }

    for raster_value, record in value_records.items():
        field_value = record.get(field)

        if pd.isna(field_value):
            cell_mask = source_array == raster_value
            fia_found = bool(record.get("_FIA_PLOT_FOUND", False))
            reason_masks["MISSING_FIA_PLT_CN"][cell_mask] = not fia_found
            reason_masks["FIA_ATTRIBUTE_MISSING"][cell_mask] = (
                fia_found
                and not bool(record.get(f"_FIA_ATTRIBUTE_AVAILABLE_{field}", False))
            )
            reason_masks["TREEMAP_RAT_ATTRIBUTE_MISSING"][cell_mask] = (
                not bool(record.get(f"_RAT_ATTRIBUTE_AVAILABLE_{field}", False))
            )
            if field in {"QMD", "TPA"} and fia_found:
                live_records = record.get("LIVE_TREE_RECORDS")
                reason_masks["NO_LIVE_TREES"][cell_mask] = (
                    pd.isna(live_records) or float(live_records or 0) <= 0
                )
            continue

        stored_value = (
            float(field_value)
            * decimal_multiplier
        )

        cell_mask = source_array == raster_value

        if np.any(cell_mask):
            working_array[cell_mask] = stored_value
            provenance[cell_mask] = int(record.get(f"_SOURCE_{field}", 1))
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
        provenance[missing_mask] = 3

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

    distance_map_units = np.full(source_array.shape, np.nan, dtype=np.float64)
    donor_cell_id = np.full(source_array.shape, -1, dtype=np.int64)
    own_ids = np.arange(source_array.size, dtype=np.int64).reshape(source_array.shape)
    donor_cell_id[assigned_before_fill] = own_ids[assigned_before_fill]
    if missing_before_fill > 0:
        rows, cols = np.indices(source_array.shape)
        donor_rows = nearest_indices[0]
        donor_cols = nearest_indices[1]
        distances = np.sqrt(
            np.square((rows - donor_rows) * abs(cell_height))
            + np.square((cols - donor_cols) * abs(cell_width))
        )
        distance_map_units[missing_mask] = distances[missing_mask]
        nearest_flat = np.ravel_multi_index(
            (donor_rows, donor_cols), source_array.shape
        )
        donor_cell_id[missing_mask] = nearest_flat[missing_mask]

    return {
        "matched_cells": matched_cells,
        "provenance": provenance,
        "nearest_distance": distance_map_units,
        "donor_cell_id": donor_cell_id,
        "source_valid": source_valid,
        "reason_counts": {
            reason: int((mask & missing_mask).sum())
            for reason, mask in reason_masks.items()
        },
    }




def write_stage1_imputation_qa(
    qa_results, lower_left, cell_width, cell_height, spatial_reference
):
    """Write cell-source, distance, landscape, and donor-consistency QA."""
    if spatial_reference is None or spatial_reference.type != "Projected":
        raise RuntimeError("Stage 1 acreage and distance QA requires a projected raster.")
    meters_per_unit = float(spatial_reference.metersPerUnit)
    cell_acres = (
        abs(cell_width * cell_height) * meters_per_unit * meters_per_unit
        / 4046.8564224
    )
    units = str(spatial_reference.linearUnitName or "Map units")
    variable_rows = []
    provenance_stack, donor_stack = [], []
    for variable, result in qa_results.items():
        provenance = result["provenance"]
        valid = result["source_valid"]
        distances = result["nearest_distance"][provenance == 3]
        row = {
            "VARIABLE": variable,
            "TOTAL_VALID_TREEMAP_CELLS": int(valid.sum()),
            "DISTANCE_UNITS": units,
        }
        for code, label in (
            (1, "FIA_DIRECT"), (2, "TREEMAP_RAT_FALLBACK"),
            (3, "NEAREST_IMPUTED"),
        ):
            count = int((provenance == code).sum())
            row[f"{label}_CELLS"] = count
            row[f"{label}_ACRES"] = count * cell_acres
            row[f"PCT_{label}"] = 100.0 * count / int(valid.sum())
        if distances.size:
            row.update({
                "MIN_DISTANCE": float(np.min(distances)),
                "MEDIAN_DISTANCE": float(np.median(distances)),
                "MEAN_DISTANCE": float(np.mean(distances)),
                "P90_DISTANCE": float(np.percentile(distances, 90)),
                "P95_DISTANCE": float(np.percentile(distances, 95)),
                "MAX_DISTANCE": float(np.max(distances)),
            })
        else:
            for name in ("MIN", "MEDIAN", "MEAN", "P90", "P95", "MAX"):
                row[f"{name}_DISTANCE"] = 0.0
        for low, high, label in (
            (0, 30, "0_30"), (30, 60, "30_60"), (60, 120, "60_120"),
            (120, 300, "120_300"), (300, 600, "300_600"),
            (600, 1000, "600_1000"), (1000, np.inf, "GT_1000"),
        ):
            count = int(((distances > low) & (distances <= high)).sum())
            row[f"DIST_{label}_CELLS"] = count
            row[f"DIST_{label}_ACRES"] = count * cell_acres
        row.update({
            f"{reason}_CELLS": count
            for reason, count in result["reason_counts"].items()
        })
        variable_rows.append(row)
        provenance_stack.append(provenance)
        donor_stack.append(result["donor_cell_id"])
        if stage1_qa_provenance_rasters:
            output = os.path.join(output_folder, f"QA_{variable}_Source.tif")
            raster = arcpy.NumPyArrayToRaster(
                provenance, lower_left, cell_width, cell_height, value_to_nodata=0
            )
            raster.save(output)
            arcpy.management.DefineProjection(output, spatial_reference)

    pd.DataFrame(variable_rows).to_csv(stage1_qa_summary_csv, index=False)
    sources = np.stack(provenance_stack)
    donors = np.stack(donor_stack)
    valid = qa_results[next(iter(qa_results))]["source_valid"]
    nearest_count = (sources == 3).sum(axis=0)
    all_fia = valid & np.all(sources == 1, axis=0)
    any_rat = valid & np.any(sources == 2, axis=0)
    any_nearest = valid & (nearest_count > 0)
    multi_nearest = nearest_count >= 2
    same_donor = np.zeros(valid.shape, dtype=bool)
    different_donor = np.zeros(valid.shape, dtype=bool)
    for row, col in np.argwhere(multi_nearest):
        ids = donors[:, row, col][sources[:, row, col] == 3]
        if len(set(ids.tolist())) == 1:
            same_donor[row, col] = True
        else:
            different_donor[row, col] = True
    landscape_rows = []
    total = int(valid.sum())
    for metric, mask in (
        ("TOTAL_VALID_TREEMAP_CELLS", valid),
        ("CELLS_WITH_ALL_FIA_DIRECT_ATTRIBUTES", all_fia),
        ("CELLS_WITH_ANY_TREEMAP_RAT_FALLBACK", any_rat),
        ("CELLS_WITH_ANY_NEAREST_CELL_IMPUTATION", any_nearest),
        ("MULTI_IMPUTED_SAME_DONOR_CELL", same_donor),
        ("MULTI_IMPUTED_DIFFERENT_DONOR_CELLS", different_donor),
    ):
        count = int(mask.sum())
        landscape_rows.append({
            "METRIC": metric, "CELL_COUNT": count, "ACRES": count * cell_acres,
            "PERCENT_OF_VALID_TREEMAP_ACRES": 100.0 * count / total,
        })
    for count_value in range(1, 6):
        count = int((nearest_count == count_value).sum())
        landscape_rows.append({
            "METRIC": f"CELLS_WITH_{count_value}_NEAREST_IMPUTED_VARIABLES",
            "CELL_COUNT": count, "ACRES": count * cell_acres,
            "PERCENT_OF_VALID_TREEMAP_ACRES": 100.0 * count / total,
        })
    pd.DataFrame(landscape_rows).to_csv(stage1_qa_landscape_csv, index=False)
    if stage1_qa_provenance_rasters:
        output = os.path.join(output_folder, "QA_Any_Nearest_Imputation.tif")
        raster = arcpy.NumPyArrayToRaster(
            nearest_count.astype(np.uint8), lower_left, cell_width, cell_height,
            value_to_nodata=0,
        )
        raster.save(output)
        arcpy.management.DefineProjection(output, spatial_reference)
    message("PRE-SEGMENTATION IMPUTATION QA")
    message("-" * 50)
    message(f"Valid TreeMap cells: {total:,}")
    message(f"Valid TreeMap acres: {total * cell_acres:,.2f}")
    message(
        f"Cells with any nearest fill: {int(any_nearest.sum()):,} "
        f"({100.0 * any_nearest.sum() / total:.2f}%)"
    )
    message(
        "Cells using multiple donor cells for imputed attributes: "
        f"{int(different_donor.sum()):,}"
    )




def configure(config):
    global treemap_raster, treemap_lookup_raster, fia_tree_files, fia_cond_files
    global output_folder, treemap_plt_field, clip_boundary, working_gdb
    global lookup_csv, age_raster, forest_type_raster, basal_area_raster
    global qmd_raster, tpa_raster, source_nodata
    global forest_condition_status_code, live_tree_status_code
    global output_nodata, coverage_check_nodata, stored_value_multipliers
    global stage1_qa_summary_csv, stage1_qa_landscape_csv
    global stage1_qa_provenance_rasters

    treemap_raster = config["treemap_raster"]
    treemap_lookup_raster = config.get("treemap_lookup_raster") or treemap_raster
    fia_tree_files = list(config["fia_tree_files"])
    fia_cond_files = list(config["fia_cond_files"])
    output_folder = config["output_folder"]
    treemap_plt_field = config.get("treemap_plt_field", "PLT_CN")
    clip_boundary = config.get("clip_boundary")
    working_gdb = config["working_gdb"]
    forest_condition_status_code = config["forest_condition_status_code"]
    live_tree_status_code = config["live_tree_status_code"]
    source_nodata = config["source_nodata"]
    output_nodata = config["output_nodata"]
    coverage_check_nodata = config["coverage_check_nodata"]
    stored_value_multipliers = dict(config["stored_value_multipliers"])
    stage1_qa_provenance_rasters = bool(
        config.get("write_qa_provenance_rasters", True)
    )

    lookup_csv = os.path.join(output_folder, "TreeMap_PLT_CN_Segmentation_Lookup.csv")
    age_raster = os.path.join(output_folder, "TreeMap_STDAGE.tif")
    forest_type_raster = os.path.join(output_folder, "TreeMap_FORTYPCD.tif")
    basal_area_raster = os.path.join(output_folder, "TreeMap_BALIVE.tif")
    qmd_raster = os.path.join(output_folder, "TreeMap_QMD.tif")
    tpa_raster = os.path.join(output_folder, "TreeMap_TPA.tif")
    stage1_qa_summary_csv = os.path.join(
        output_folder, "QA_STAGE1_Raster_Imputation_Summary.csv"
    )
    stage1_qa_landscape_csv = os.path.join(
        output_folder, "QA_STAGE1_Landscape_Imputation_Summary.csv"
    )

    arcpy.env.overwriteOutput = True
    os.makedirs(output_folder, exist_ok=True)
    if not arcpy.Exists(working_gdb):
        parent = os.path.dirname(working_gdb)
        os.makedirs(parent, exist_ok=True)
        arcpy.management.CreateFileGDB(parent, os.path.basename(working_gdb))


def run_stage(config):
    configure(config)
    # ============================================================
    # VALIDATE CORE INPUTS
    # ============================================================

    message("=" * 60)
    message("VALIDATING CORE INPUTS")
    message("=" * 60)

    for path, label in [
        (treemap_raster, "TreeMap spatial raster"),
        (treemap_lookup_raster, "TreeMap attribute lookup raster"),
        (clip_boundary, "Maine clip boundary"),
    ]:
        if path is not None and not arcpy.Exists(path):
            raise FileNotFoundError(
                f"{label} was not found:\n{path}"
            )

    for path in fia_tree_files + fia_cond_files:
        if not os.path.exists(path):
            raise FileNotFoundError(
                f"FIA CSV file was not found:\n{path}"
            )

    message("Core input validation passed.")


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

    # Identify only the TreeMap VALUE codes that occur inside the Maine
    # study-area raster. The full CONUS attribute table is then filtered
    # to this set.
    source_values = set(
        np.unique(
            source_array[source_array != source_nodata]
        ).astype(np.int64)
    )

    message(
        f"Unique TreeMap VALUE codes in clipped study area: "
        f"{len(source_values):,}"
    )

    value_records = read_raster_value_lookup(
        treemap_lookup_raster,
        lookup,
        required_values=source_values,
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

    qa_results = {}

    # Age in whole years.
    qa_results["STDAGE"] = write_lookup_raster(
        source_array=source_array,
        value_records=value_records,
        field="STDAGE",
        output_path=age_raster,
        lower_left=lower_left,
        cell_width=cell_width,
        cell_height=cell_height,
        spatial_reference=spatial_reference,
        output_dtype=np.int16,
        output_nodata=output_nodata,
        decimal_multiplier=stored_value_multipliers["STDAGE"],
    )

    # FIA forest type code.
    qa_results["FORTYPCD"] = write_lookup_raster(
        source_array=source_array,
        value_records=value_records,
        field="FORTYPCD",
        output_path=forest_type_raster,
        lower_left=lower_left,
        cell_width=cell_width,
        cell_height=cell_height,
        spatial_reference=spatial_reference,
        output_dtype=np.int16,
        output_nodata=output_nodata,
        decimal_multiplier=stored_value_multipliers["FORTYPCD"],
    )

    # BALIVE is stored in tenths of square feet per acre.
    #
    # Raster value 287 means 28.7 ft2/ac.
    qa_results["BALIVE"] = write_lookup_raster(
        source_array=source_array,
        value_records=value_records,
        field="BALIVE",
        output_path=basal_area_raster,
        lower_left=lower_left,
        cell_width=cell_width,
        cell_height=cell_height,
        spatial_reference=spatial_reference,
        output_dtype=np.int32,
        output_nodata=output_nodata,
        decimal_multiplier=stored_value_multipliers["BALIVE"],
    )

    # QMD is stored in tenths of an inch.
    #
    # Raster value 112 means 11.2 inches.
    qa_results["QMD"] = write_lookup_raster(
        source_array=source_array,
        value_records=value_records,
        field="QMD",
        output_path=qmd_raster,
        lower_left=lower_left,
        cell_width=cell_width,
        cell_height=cell_height,
        spatial_reference=spatial_reference,
        output_dtype=np.int16,
        output_nodata=output_nodata,
        decimal_multiplier=stored_value_multipliers["QMD"],
    )

    # TPA is stored in tenths of a tree per acre.
    #
    # Raster value 1247 means 124.7 trees/ac.
    qa_results["TPA"] = write_lookup_raster(
        source_array=source_array,
        value_records=value_records,
        field="TPA",
        output_path=tpa_raster,
        lower_left=lower_left,
        cell_width=cell_width,
        cell_height=cell_height,
        spatial_reference=spatial_reference,
        output_dtype=np.int32,
        output_nodata=output_nodata,
        decimal_multiplier=stored_value_multipliers["TPA"],
    )

    write_stage1_imputation_qa(
        qa_results, lower_left, cell_width, cell_height, spatial_reference
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
            nodata_to_value=coverage_check_nodata,
        )

        assigned_count = int(
            np.count_nonzero(
                check_array != coverage_check_nodata
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
    return {
        "source_raster": source_raster,
        "working_gdb": working_gdb,
        "lookup_csv": lookup_csv,
        "age_raster": age_raster,
        "forest_type_raster": forest_type_raster,
        "basal_area_raster": basal_area_raster,
        "qmd_raster": qmd_raster,
        "tpa_raster": tpa_raster,
    }
