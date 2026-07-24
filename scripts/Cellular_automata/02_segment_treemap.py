"""
02_segment_treemap.py

Prototype TreeMap stand segmentation for LETO.

Inputs
------
Five aligned rasters created by 01_build_segmentation_rasters.py:
    TreeMap_STDAGE.tif
    TreeMap_FORTYPCD.tif
    TreeMap_BALIVE.tif   (stored value = actual ft2/ac * 10)
    TreeMap_QMD.tif      (stored value = actual inches * 10)
    TreeMap_TPA.tif      (stored value = actual trees/ac * 10)

Outputs
-------
    TreeMap_Segments_Initial.tif
    TreeMap_Segments_Initial  (polygon feature class)
    TreeMap_Segment_Summary.csv

Method
------
1. Read aligned raster arrays.
2. Standardize continuous variables.
3. Create regular block seeds of approximately INITIAL_SEED_ACRES.
4. Iteratively reassign boundary cells to the most similar adjacent segment.
5. Separate disconnected segment pieces.
6. Merge segments smaller than MIN_STAND_ACRES.
7. Export raster, polygons, and segment summaries.

This is an intentionally simple prototype. It is meant to test whether the
TreeMap/FIA raster stack produces believable stand-like units before adding
more complex shape, ownership, and maximum-area controls.
"""

import math
import os
from collections import Counter, defaultdict

import arcpy
import numpy as np
import pandas as pd
from scipy import ndimage


# ============================================================
# USER INPUTS
# ============================================================

input_folder = (
    r"D:\forest_condition_2026\Cellular_automata"
    r"\TreeMap_Segmentation_Rasters"
)

age_raster = os.path.join(input_folder, "TreeMap_STDAGE.tif")
forest_type_raster = os.path.join(input_folder, "TreeMap_FORTYPCD.tif")
basal_area_raster = os.path.join(input_folder, "TreeMap_BALIVE.tif")
qmd_raster = os.path.join(input_folder, "TreeMap_QMD.tif")
tpa_raster = os.path.join(input_folder, "TreeMap_TPA.tif")

output_folder = os.path.join(
    r"D:\forest_condition_2026\Cellular_automata",
    "TreeMap_Segmentation_Output",
)

output_gdb = os.path.join(
    output_folder,
    "TreeMap_Segmentation.gdb",
)

segment_raster = os.path.join(
    output_folder,
    "TreeMap_Segments_Initial.tif",
)

segment_polygons = os.path.join(
    output_gdb,
    "TreeMap_Segments_Initial",
)

segment_summary_csv = os.path.join(
    output_folder,
    "TreeMap_Segment_Summary.csv",
)


# ============================================================
# SEGMENTATION SETTINGS
# ============================================================

INITIAL_SEED_ACRES = 20.0
MIN_STAND_ACRES = 10.0
MAX_STAND_ACRES = 120.0

MAX_ITERATIONS = 20
CONVERGENCE_THRESHOLD = 0.001

# A cell can move only if the candidate segment improves the score by at
# least this amount. A small positive threshold reduces oscillation.
MINIMUM_SCORE_IMPROVEMENT = 0.01

# Four-neighbor adjacency is preferable for stand delineation because
# diagonal contact alone should not join two stands.
USE_EIGHT_NEIGHBORS = False

# Relative importance of each segmentation variable.
WEIGHTS = {
    "FORTYPCD": 0.30,
    "STDAGE": 0.25,
    "BALIVE": 0.20,
    "QMD": 0.15,
    "TPA": 0.10,
}

# Small preference for a segment sharing more cell edges with the candidate
# cell. This encourages smoother boundaries.
SHARED_EDGE_BONUS = 0.05

# Continuous-variable outlier clipping after standardization.
Z_CLIP = 4.0

# Raster NoData sentinel used when reading arrays.
READ_NODATA = -999999


# ============================================================
# GENERAL SETTINGS
# ============================================================

arcpy.env.overwriteOutput = True

os.makedirs(output_folder, exist_ok=True)

if not arcpy.Exists(output_gdb):
    arcpy.management.CreateFileGDB(
        output_folder,
        os.path.basename(output_gdb),
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


def validate_inputs():
    missing = [
        path
        for path in [
            age_raster,
            forest_type_raster,
            basal_area_raster,
            qmd_raster,
            tpa_raster,
        ]
        if not arcpy.Exists(path)
    ]

    if missing:
        raise FileNotFoundError(
            "The following input rasters were not found:\n"
            + "\n".join(missing)
        )


def raster_geometry(raster_path):
    desc = arcpy.Describe(raster_path)

    return {
        "rows": int(desc.height),
        "cols": int(desc.width),
        "cell_width": float(desc.meanCellWidth),
        "cell_height": float(desc.meanCellHeight),
        "xmin": float(desc.extent.XMin),
        "ymin": float(desc.extent.YMin),
        "xmax": float(desc.extent.XMax),
        "ymax": float(desc.extent.YMax),
        "spatial_reference": desc.spatialReference,
    }


def validate_alignment(reference_path, other_paths):
    reference = raster_geometry(reference_path)

    for raster_path in other_paths:
        current = raster_geometry(raster_path)

        checks = [
            reference["rows"] == current["rows"],
            reference["cols"] == current["cols"],
            math.isclose(
                reference["cell_width"],
                current["cell_width"],
                rel_tol=0,
                abs_tol=1e-8,
            ),
            math.isclose(
                reference["cell_height"],
                current["cell_height"],
                rel_tol=0,
                abs_tol=1e-8,
            ),
            math.isclose(
                reference["xmin"],
                current["xmin"],
                rel_tol=0,
                abs_tol=1e-6,
            ),
            math.isclose(
                reference["ymin"],
                current["ymin"],
                rel_tol=0,
                abs_tol=1e-6,
            ),
        ]

        if not all(checks):
            raise ValueError(
                f"Raster is not aligned with the reference raster:\n"
                f"{raster_path}"
            )

    return reference


def read_raster(path, scale=1.0):
    array = arcpy.RasterToNumPyArray(
        path,
        nodata_to_value=READ_NODATA,
    ).astype(np.float64)

    array[array == READ_NODATA] = np.nan

    if scale != 1.0:
        array /= scale

    return array


def robust_standardize(array, valid_mask):
    """
    Standardize with median and interquartile range.

    This is less sensitive than mean/standard deviation to unusually dense
    FIA plots with very high TPA or basal area.
    """
    values = array[valid_mask]

    median = float(np.nanmedian(values))
    q1 = float(np.nanpercentile(values, 25))
    q3 = float(np.nanpercentile(values, 75))
    iqr = q3 - q1

    # Fall back to standard deviation if the IQR is zero.
    if not np.isfinite(iqr) or iqr == 0:
        scale = float(np.nanstd(values))
    else:
        scale = iqr / 1.349

    if not np.isfinite(scale) or scale == 0:
        scale = 1.0

    standardized = (array - median) / scale
    standardized = np.clip(standardized, -Z_CLIP, Z_CLIP)

    return standardized, median, scale


def neighbor_offsets():
    if USE_EIGHT_NEIGHBORS:
        return [
            (-1, -1), (-1, 0), (-1, 1),
            (0, -1),           (0, 1),
            (1, -1),  (1, 0),  (1, 1),
        ]

    return [
        (-1, 0),
        (0, -1),
        (0, 1),
        (1, 0),
    ]


def create_initial_segments(valid_mask, cells_per_seed):
    """
    Create regular square block seeds, while preserving raster gaps.
    """
    rows, cols = valid_mask.shape
    block_side = max(1, int(round(math.sqrt(cells_per_seed))))

    labels = np.zeros((rows, cols), dtype=np.int32)
    next_id = 1

    for row_start in range(0, rows, block_side):
        row_end = min(rows, row_start + block_side)

        for col_start in range(0, cols, block_side):
            col_end = min(cols, col_start + block_side)

            block_mask = valid_mask[
                row_start:row_end,
                col_start:col_end,
            ]

            if not np.any(block_mask):
                continue

            # A block can contain several disconnected forest pieces because
            # roads and water appear as NoData. Give each piece its own ID.
            local_labels, component_count = ndimage.label(
                block_mask,
                structure=np.array(
                    [
                        [0, 1, 0],
                        [1, 1, 1],
                        [0, 1, 0],
                    ],
                    dtype=np.int8,
                ),
            )

            target = labels[
                row_start:row_end,
                col_start:col_end,
            ]

            for component_id in range(1, component_count + 1):
                target[local_labels == component_id] = next_id
                next_id += 1

    return labels


def segment_statistics(labels, feature_arrays, forest_type):
    """
    Calculate continuous means and dominant forest type for every segment.
    """
    segment_ids = np.unique(labels)
    segment_ids = segment_ids[segment_ids > 0]

    stats = {}

    for segment_id in segment_ids:
        mask = labels == segment_id
        pixel_count = int(mask.sum())

        record = {
            "PIXEL_COUNT": pixel_count,
        }

        for name, array in feature_arrays.items():
            record[name] = float(np.nanmean(array[mask]))

        types = forest_type[mask]
        types = types[np.isfinite(types)]

        if len(types):
            type_counts = Counter(types.astype(np.int32))
            dominant_type, dominant_count = type_counts.most_common(1)[0]
            record["FORTYPCD"] = int(dominant_type)
            record["FORTYPE_COUNT"] = int(dominant_count)
        else:
            record["FORTYPCD"] = None
            record["FORTYPE_COUNT"] = 0

        stats[int(segment_id)] = record

    return stats


def candidate_segments(labels, row, col, offsets):
    rows, cols = labels.shape
    counts = Counter()

    for dr, dc in offsets:
        rr = row + dr
        cc = col + dc

        if rr < 0 or rr >= rows or cc < 0 or cc >= cols:
            continue

        candidate = int(labels[rr, cc])

        if candidate > 0:
            counts[candidate] += 1

    return counts


def assignment_cost(
    row,
    col,
    candidate_id,
    candidate_edge_count,
    stats,
    feature_arrays,
    forest_type,
):
    candidate = stats[candidate_id]

    cost = 0.0

    for name in ["STDAGE", "BALIVE", "QMD", "TPA"]:
        difference = (
            feature_arrays[name][row, col]
            - candidate[name]
        )

        cost += WEIGHTS[name] * difference * difference

    cell_type = int(forest_type[row, col])
    type_penalty = (
        0.0
        if candidate["FORTYPCD"] == cell_type
        else WEIGHTS["FORTYPCD"]
    )

    cost += type_penalty

    # Lower cost for candidates sharing more edges with the cell.
    cost -= SHARED_EDGE_BONUS * candidate_edge_count

    return float(cost)


def identify_boundary_cells(labels, valid_mask):
    """
    Return a Boolean mask identifying cells touching a different segment.
    """
    boundary = np.zeros(labels.shape, dtype=bool)

    boundary[:-1, :] |= (
        valid_mask[:-1, :]
        & valid_mask[1:, :]
        & (labels[:-1, :] != labels[1:, :])
    )
    boundary[1:, :] |= (
        valid_mask[1:, :]
        & valid_mask[:-1, :]
        & (labels[1:, :] != labels[:-1, :])
    )
    boundary[:, :-1] |= (
        valid_mask[:, :-1]
        & valid_mask[:, 1:]
        & (labels[:, :-1] != labels[:, 1:])
    )
    boundary[:, 1:] |= (
        valid_mask[:, 1:]
        & valid_mask[:, :-1]
        & (labels[:, 1:] != labels[:, :-1])
    )

    return boundary


def iterative_reassignment(
    labels,
    valid_mask,
    feature_arrays,
    forest_type,
):
    """
    Reassign boundary cells using synchronous updates.

    Synchronous updating makes the result less dependent on traversal order.
    """
    offsets = neighbor_offsets()
    total_valid_cells = int(valid_mask.sum())

    for iteration in range(1, MAX_ITERATIONS + 1):
        stats = segment_statistics(
            labels,
            feature_arrays,
            forest_type,
        )

        boundary_mask = identify_boundary_cells(
            labels,
            valid_mask,
        )

        boundary_rows, boundary_cols = np.where(boundary_mask)

        new_labels = labels.copy()
        changed = 0

        for row, col in zip(boundary_rows, boundary_cols):
            current_id = int(labels[row, col])

            candidates = candidate_segments(
                labels,
                row,
                col,
                offsets,
            )

            if not candidates:
                continue

            # Include the current segment even if none of the sampled
            # neighbors happens to have the current ID.
            if current_id not in candidates:
                candidates[current_id] = 0

            current_cost = assignment_cost(
                row,
                col,
                current_id,
                candidates[current_id],
                stats,
                feature_arrays,
                forest_type,
            )

            best_id = current_id
            best_cost = current_cost

            for candidate_id, edge_count in candidates.items():
                if candidate_id == current_id:
                    continue

                candidate_cost = assignment_cost(
                    row,
                    col,
                    candidate_id,
                    edge_count,
                    stats,
                    feature_arrays,
                    forest_type,
                )

                if (
                    candidate_cost
                    < best_cost - MINIMUM_SCORE_IMPROVEMENT
                ):
                    best_cost = candidate_cost
                    best_id = candidate_id

            if best_id != current_id:
                new_labels[row, col] = best_id
                changed += 1

        labels = new_labels

        changed_fraction = (
            changed / total_valid_cells
            if total_valid_cells
            else 0.0
        )

        message(
            f"Iteration {iteration:02d}: "
            f"{changed:,} cells changed "
            f"({changed_fraction:.4%})"
        )

        if changed_fraction <= CONVERGENCE_THRESHOLD:
            message("Convergence threshold reached.")
            break

    return labels


def split_disconnected_segments(labels):
    """
    Ensure every segment ID represents one four-neighbor-connected region.
    """
    output = np.zeros(labels.shape, dtype=np.int32)
    next_id = 1

    structure = np.array(
        [
            [0, 1, 0],
            [1, 1, 1],
            [0, 1, 0],
        ],
        dtype=np.int8,
    )

    for old_id in np.unique(labels):
        if old_id <= 0:
            continue

        pieces, piece_count = ndimage.label(
            labels == old_id,
            structure=structure,
        )

        for piece_id in range(1, piece_count + 1):
            output[pieces == piece_id] = next_id
            next_id += 1

    return output


def segment_adjacency(labels):
    """
    Return shared-edge counts for each adjacent segment pair.
    """
    adjacency = defaultdict(Counter)

    left = labels[:, :-1]
    right = labels[:, 1:]

    mask = (
        (left > 0)
        & (right > 0)
        & (left != right)
    )

    for first, second in zip(left[mask], right[mask]):
        first = int(first)
        second = int(second)
        adjacency[first][second] += 1
        adjacency[second][first] += 1

    upper = labels[:-1, :]
    lower = labels[1:, :]

    mask = (
        (upper > 0)
        & (lower > 0)
        & (upper != lower)
    )

    for first, second in zip(upper[mask], lower[mask]):
        first = int(first)
        second = int(second)
        adjacency[first][second] += 1
        adjacency[second][first] += 1

    return adjacency


def merge_small_segments(
    labels,
    feature_arrays,
    forest_type,
    minimum_cells,
    maximum_cells=None,
):
    """
    Iteratively merge undersized segments into the most similar adjacent
    segment. Shared-border length is used as a small tie-breaker.

    When maximum_cells is supplied, a merge is skipped if the combined
    segment would exceed the maximum stand size. Any undersized segments
    that cannot be merged safely are left for the final deletion step.
    """
    pass_number = 0

    while True:
        pass_number += 1

        stats = segment_statistics(
            labels,
            feature_arrays,
            forest_type,
        )

        small_ids = [
            segment_id
            for segment_id, record in stats.items()
            if record["PIXEL_COUNT"] < minimum_cells
        ]

        if not small_ids:
            break

        adjacency = segment_adjacency(labels)
        merged_this_pass = 0

        # Process smallest first.
        small_ids.sort(
            key=lambda segment_id: stats[segment_id]["PIXEL_COUNT"]
        )

        for small_id in small_ids:
            if small_id not in stats:
                continue

            neighbors = adjacency.get(small_id)

            if not neighbors:
                continue

            source = stats[small_id]

            best_neighbor = None
            best_cost = float("inf")

            for neighbor_id, shared_edges in neighbors.items():
                if neighbor_id not in stats:
                    continue

                target = stats[neighbor_id]

                if (
                    maximum_cells is not None
                    and (
                        source["PIXEL_COUNT"]
                        + target["PIXEL_COUNT"]
                    ) > maximum_cells
                ):
                    continue

                cost = 0.0

                for name in ["STDAGE", "BALIVE", "QMD", "TPA"]:
                    difference = source[name] - target[name]
                    cost += WEIGHTS[name] * difference * difference

                if source["FORTYPCD"] != target["FORTYPCD"]:
                    cost += WEIGHTS["FORTYPCD"]

                cost -= SHARED_EDGE_BONUS * shared_edges

                if cost < best_cost:
                    best_cost = cost
                    best_neighbor = neighbor_id

            if best_neighbor is None:
                continue

            labels[labels == small_id] = best_neighbor
            merged_this_pass += 1

        message(
            f"Small-segment merge pass {pass_number}: "
            f"{merged_this_pass:,} segments merged"
        )

        if merged_this_pass == 0:
            break

        labels = split_disconnected_segments(labels)

    return labels


def renumber_segments(labels):
    output = np.zeros(labels.shape, dtype=np.int32)

    valid_ids = np.unique(labels)
    valid_ids = valid_ids[valid_ids > 0]

    for new_id, old_id in enumerate(valid_ids, start=1):
        output[labels == old_id] = new_id

    return output


def _connected_components(mask):
    structure = np.array(
        [
            [0, 1, 0],
            [1, 1, 1],
            [0, 1, 0],
        ],
        dtype=np.int8,
    )
    return ndimage.label(mask, structure=structure)


def _bisect_connected_mask(mask):
    """
    Split one connected region near its spatial median.

    The longer row/column extent is used first. Any disconnected pieces
    created by the cut are returned separately and are handled recursively.
    """
    rows, cols = np.where(mask)

    if len(rows) < 2:
        return [mask]

    row_span = int(rows.max() - rows.min())
    col_span = int(cols.max() - cols.min())

    if col_span >= row_span:
        coordinate = cols
        cut = int(np.median(coordinate))
        first = mask.copy()
        second = mask.copy()

        grid_cols = np.indices(mask.shape)[1]
        first[grid_cols > cut] = False
        second[grid_cols <= cut] = False
    else:
        coordinate = rows
        cut = int(np.median(coordinate))
        first = mask.copy()
        second = mask.copy()

        grid_rows = np.indices(mask.shape)[0]
        first[grid_rows > cut] = False
        second[grid_rows <= cut] = False

    pieces = []

    for part in [first, second]:
        if not np.any(part):
            continue

        component_labels, component_count = _connected_components(part)

        for component_id in range(1, component_count + 1):
            component = component_labels == component_id
            if np.any(component):
                pieces.append(component)

    # Rare fallback when an unusual geometry cannot be split by the median.
    if len(pieces) < 2:
        order = np.lexsort((cols, rows))
        half = len(order) // 2

        first = np.zeros(mask.shape, dtype=bool)
        second = np.zeros(mask.shape, dtype=bool)

        first[rows[order[:half]], cols[order[:half]]] = True
        second[rows[order[half:]], cols[order[half:]]] = True

        pieces = [first, second]

    return pieces


def split_oversized_segments(labels, maximum_cells):
    """
    Recursively spatially bisect segments larger than maximum_cells.

    This is intentionally simple and deterministic. It enforces the acreage
    ceiling without adding more biological variables to the segmentation
    objective.
    """
    output = np.zeros(labels.shape, dtype=np.int32)
    next_id = 1
    split_count = 0

    for segment_id in np.unique(labels):
        if segment_id <= 0:
            continue

        queue = [labels == segment_id]

        while queue:
            mask = queue.pop()
            cell_count = int(mask.sum())

            if cell_count <= maximum_cells or cell_count <= 1:
                output[mask] = next_id
                next_id += 1
                continue

            pieces = _bisect_connected_mask(mask)

            if len(pieces) < 2:
                output[mask] = next_id
                next_id += 1
                continue

            split_count += 1
            queue.extend(pieces)

    message(
        f"Oversized-segment split operations: "
        f"{split_count:,}"
    )

    return output


def delete_small_segments(labels, minimum_cells):
    """
    Delete isolated or undersized stands instead of forcing them into a
    dissimilar neighbor. Deleted cells become NoData in the final raster.
    """
    output = labels.copy()
    deleted_segments = 0
    deleted_cells = 0

    segment_ids, counts = np.unique(
        labels[labels > 0],
        return_counts=True,
    )

    for segment_id, cell_count in zip(segment_ids, counts):
        if int(cell_count) < minimum_cells:
            output[labels == segment_id] = 0
            deleted_segments += 1
            deleted_cells += int(cell_count)

    return output, deleted_segments, deleted_cells


def create_summary(
    labels,
    raw_arrays,
    forest_type,
    cell_acres,
    homogeneity_scales,
):
    """
    Create stand summaries and a 0-100 landscape-relative homogeneity score.

    For each continuous variable:

        homogeneity = 1 - (stand SD / landscape SD)

    Values are clipped to the range 0-1. Therefore:

        100 = no within-stand variation
          0 = stand variation equals or exceeds landscape variation

    Forest-type homogeneity is the proportion of stand cells belonging to
    the dominant FORTYPCD. The final score uses the same weights as the
    segmentation model.

    The raw within-stand means, standard deviations, minima, and maxima are
    retained so the custom score does not replace directly interpretable
    statistics.
    """
    records = []

    for segment_id in np.unique(labels):
        if segment_id <= 0:
            continue

        mask = labels == segment_id
        pixel_count = int(mask.sum())

        record = {
            "SEG_ID": int(segment_id),
            "PIXEL_COUNT": pixel_count,
            "Acres": pixel_count * cell_acres,
        }

        component_scores = {}

        for name, array in raw_arrays.items():
            values = array[mask]

            mean_value = float(np.nanmean(values))
            sd_value = float(np.nanstd(values))

            record[f"{name}_MEAN"] = mean_value
            record[f"{name}_SD"] = sd_value
            record[f"{name}_MIN"] = float(np.nanmin(values))
            record[f"{name}_MAX"] = float(np.nanmax(values))

            landscape_sd = float(homogeneity_scales[name])

            if not np.isfinite(landscape_sd) or landscape_sd <= 0:
                component = 1.0 if sd_value == 0 else 0.0
            else:
                component = 1.0 - (sd_value / landscape_sd)
                component = float(np.clip(component, 0.0, 1.0))

            component_scores[name] = component
            record[f"{name}_HOM"] = 100.0 * component

        types = forest_type[mask]
        types = types[np.isfinite(types)].astype(np.int32)

        type_counts = Counter(types)
        dominant_type, dominant_count = type_counts.most_common(1)[0]

        forest_type_fraction = dominant_count / pixel_count

        record["FORTYPCD_DOM"] = int(dominant_type)
        record["FORTYPE_PCT"] = 100.0 * forest_type_fraction
        record["FORTYPCD_HOM"] = 100.0 * forest_type_fraction

        weighted_score = (
            WEIGHTS["STDAGE"] * component_scores["STDAGE"]
            + WEIGHTS["BALIVE"] * component_scores["BALIVE"]
            + WEIGHTS["QMD"] * component_scores["QMD"]
            + WEIGHTS["TPA"] * component_scores["TPA"]
            + WEIGHTS["FORTYPCD"] * forest_type_fraction
        )

        record["HOMOGENEITY"] = 100.0 * weighted_score

        records.append(record)

    return pd.DataFrame(records)


def write_segment_raster(labels, geometry):
    output_nodata = -9999

    export_array = labels.astype(np.int32)
    export_array[export_array == 0] = output_nodata

    lower_left = arcpy.Point(
        geometry["xmin"],
        geometry["ymin"],
    )

    raster = arcpy.NumPyArrayToRaster(
        export_array,
        lower_left,
        geometry["cell_width"],
        geometry["cell_height"],
        value_to_nodata=output_nodata,
    )

    raster.save(segment_raster)

    arcpy.management.DefineProjection(
        segment_raster,
        geometry["spatial_reference"],
    )

    arcpy.management.CalculateStatistics(
        segment_raster,
    )

    try:
        arcpy.management.BuildRasterAttributeTable(
            segment_raster,
            "Overwrite",
        )
    except arcpy.ExecuteError:
        message(
            "Raster attribute table was not built, but the "
            "segment raster was created."
        )


def write_segment_polygons(summary):
    if arcpy.Exists(segment_polygons):
        arcpy.management.Delete(segment_polygons)

    arcpy.conversion.RasterToPolygon(
        in_raster=segment_raster,
        out_polygon_features=segment_polygons,
        simplify="NO_SIMPLIFY",
        raster_field="Value",
        create_multipart_features="SINGLE_OUTER_PART",
    )

    existing_fields = {
        field.name.upper()
        for field in arcpy.ListFields(segment_polygons)
    }

    if "SEG_ID" not in existing_fields:
        arcpy.management.AddField(
            segment_polygons,
            "SEG_ID",
            "LONG",
        )

    value_field = None

    for field in arcpy.ListFields(segment_polygons):
        if field.name.upper() in {"GRIDCODE", "GRID_CODE"}:
            value_field = field.name
            break

    if value_field is None:
        raise RuntimeError(
            "Could not identify the raster-value field in the "
            "polygon output."
        )

    arcpy.management.CalculateField(
        segment_polygons,
        "SEG_ID",
        f"!{value_field}!",
        "PYTHON3",
    )

    # Join summary attributes through an intermediate geodatabase table.
    summary_table = os.path.join(
        output_gdb,
        "TreeMap_Segment_Summary",
    )

    if arcpy.Exists(summary_table):
        arcpy.management.Delete(summary_table)

    arcpy.conversion.TableToTable(
        segment_summary_csv,
        output_gdb,
        os.path.basename(summary_table),
    )

    join_fields = [
        field
        for field in summary.columns
        if field != "SEG_ID"
    ]

    arcpy.management.JoinField(
        segment_polygons,
        "SEG_ID",
        summary_table,
        "SEG_ID",
        join_fields,
    )


# ============================================================
# MAIN
# ============================================================

def main():
    message("=" * 60)
    message("VALIDATING INPUT RASTERS")
    message("=" * 60)

    validate_inputs()

    geometry = validate_alignment(
        age_raster,
        [
            forest_type_raster,
            basal_area_raster,
            qmd_raster,
            tpa_raster,
        ],
    )

    arcpy.env.snapRaster = age_raster
    arcpy.env.cellSize = age_raster
    arcpy.env.extent = age_raster
    arcpy.env.outputCoordinateSystem = age_raster

    cell_area_map_units = (
        geometry["cell_width"]
        * geometry["cell_height"]
    )

    # The current TreeMap data are in a meter-based projected coordinate
    # system. Convert square meters to acres.
    cell_acres = cell_area_map_units / 4046.8564224

    initial_seed_cells = max(
        1,
        int(round(INITIAL_SEED_ACRES / cell_acres)),
    )

    minimum_stand_cells = max(
        1,
        int(math.ceil(MIN_STAND_ACRES / cell_acres)),
    )

    message(
        f"Cell size: "
        f"{geometry['cell_width']:.3f} × "
        f"{geometry['cell_height']:.3f} map units"
    )
    message(f"Cell area: {cell_acres:.6f} acres")
    message(
        f"Initial seed size: approximately "
        f"{initial_seed_cells:,} cells"
    )
    message(
        f"Minimum stand size: "
        f"{minimum_stand_cells:,} cells"
    )
    message(
        f"Maximum stand size: approximately "
        f"{int(math.floor(MAX_STAND_ACRES / cell_acres)):,} cells"
    )

    message("=" * 60)
    message("READING SEGMENTATION RASTERS")
    message("=" * 60)

    age = read_raster(age_raster)
    forest_type = read_raster(forest_type_raster)
    basal_area = read_raster(
        basal_area_raster,
        scale=10.0,
    )
    qmd = read_raster(
        qmd_raster,
        scale=10.0,
    )
    tpa = read_raster(
        tpa_raster,
        scale=10.0,
    )

    valid_mask = (
        np.isfinite(age)
        & np.isfinite(forest_type)
        & np.isfinite(basal_area)
        & np.isfinite(qmd)
        & np.isfinite(tpa)
    )

    valid_cells = int(valid_mask.sum())

    if valid_cells == 0:
        raise RuntimeError(
            "No cells contained complete data across all five rasters."
        )

    message(
        f"Cells with complete segmentation data: "
        f"{valid_cells:,}"
    )

    message("=" * 60)
    message("STANDARDIZING CONTINUOUS VARIABLES")
    message("=" * 60)

    age_z, age_center, age_scale = robust_standardize(
        age,
        valid_mask,
    )
    ba_z, ba_center, ba_scale = robust_standardize(
        basal_area,
        valid_mask,
    )
    qmd_z, qmd_center, qmd_scale = robust_standardize(
        qmd,
        valid_mask,
    )
    tpa_z, tpa_center, tpa_scale = robust_standardize(
        tpa,
        valid_mask,
    )

    message(
        f"STDAGE center={age_center:.3f}, scale={age_scale:.3f}"
    )
    message(
        f"BALIVE center={ba_center:.3f}, scale={ba_scale:.3f}"
    )
    message(
        f"QMD center={qmd_center:.3f}, scale={qmd_scale:.3f}"
    )
    message(
        f"TPA center={tpa_center:.3f}, scale={tpa_scale:.3f}"
    )

    feature_arrays = {
        "STDAGE": age_z,
        "BALIVE": ba_z,
        "QMD": qmd_z,
        "TPA": tpa_z,
    }

    raw_arrays = {
        "STDAGE": age,
        "BALIVE": basal_area,
        "QMD": qmd,
        "TPA": tpa,
    }

    message("=" * 60)
    message("CREATING INITIAL SEEDS")
    message("=" * 60)

    labels = create_initial_segments(
        valid_mask,
        initial_seed_cells,
    )

    initial_count = len(
        np.unique(labels[labels > 0])
    )

    message(
        f"Initial connected seed regions: "
        f"{initial_count:,}"
    )

    message("=" * 60)
    message("ITERATIVE CELL REASSIGNMENT")
    message("=" * 60)

    labels = iterative_reassignment(
        labels,
        valid_mask,
        feature_arrays,
        forest_type,
    )

    message("=" * 60)
    message("SEPARATING DISCONNECTED PIECES")
    message("=" * 60)

    labels = split_disconnected_segments(labels)

    pre_merge_count = len(
        np.unique(labels[labels > 0])
    )

    message(
        f"Connected segments before minimum-size merging: "
        f"{pre_merge_count:,}"
    )

    maximum_stand_cells = max(
        minimum_stand_cells,
        int(math.floor(MAX_STAND_ACRES / cell_acres)),
    )

    message("=" * 60)
    message("SPLITTING OVERSIZED SEGMENTS")
    message("=" * 60)

    labels = split_oversized_segments(
        labels,
        maximum_stand_cells,
    )

    labels = split_disconnected_segments(labels)

    count_after_split = len(
        np.unique(labels[labels > 0])
    )

    message(
        f"Connected segments after maximum-size splitting: "
        f"{count_after_split:,}"
    )

    message("=" * 60)
    message("MERGING UNDERSIZED SEGMENTS")
    message("=" * 60)

    labels = merge_small_segments(
        labels,
        feature_arrays,
        forest_type,
        minimum_stand_cells,
        maximum_cells=maximum_stand_cells,
    )

    labels = split_disconnected_segments(labels)

    remaining_small_before_delete = sum(
        1
        for count in np.unique(
            labels[labels > 0],
            return_counts=True,
        )[1]
        if int(count) < minimum_stand_cells
    )

    message(
        f"Undersized segments remaining after merge: "
        f"{remaining_small_before_delete:,}"
    )

    message("=" * 60)
    message("DELETING REMAINING UNDERSIZED SEGMENTS")
    message("=" * 60)

    labels, deleted_segments, deleted_cells = delete_small_segments(
        labels,
        minimum_stand_cells,
    )

    labels = split_disconnected_segments(labels)
    labels = renumber_segments(labels)

    message(
        f"Deleted undersized segments: "
        f"{deleted_segments:,}"
    )
    message(
        f"Deleted forest cells: "
        f"{deleted_cells:,} "
        f"({deleted_cells * cell_acres:,.2f} acres)"
    )

    final_count = len(
        np.unique(labels[labels > 0])
    )

    message(
        f"Final segment count: "
        f"{final_count:,}"
    )

    message("=" * 60)
    message("CREATING SEGMENT SUMMARY")
    message("=" * 60)

    # Use the actual standard deviation across all valid landscape cells
    # as the reference for the landscape-relative homogeneity calculation.
    homogeneity_scales = {
        "STDAGE": float(np.nanstd(age[valid_mask])),
        "BALIVE": float(np.nanstd(basal_area[valid_mask])),
        "QMD": float(np.nanstd(qmd[valid_mask])),
        "TPA": float(np.nanstd(tpa[valid_mask])),
    }

    message("Landscape standard deviations used for homogeneity:")
    for variable, landscape_sd in homogeneity_scales.items():
        message(
            f"  {variable}: {landscape_sd:.3f}"
        )

    summary = create_summary(
        labels,
        raw_arrays,
        forest_type,
        cell_acres,
        homogeneity_scales,
    )

    summary.to_csv(
        segment_summary_csv,
        index=False,
    )

    message(
        f"Segment summary written:\n"
        f"{segment_summary_csv}"
    )

    message(
        f"Mean stand acreage: "
        f"{summary['Acres'].mean():.2f}"
    )
    message(
        f"Median stand acreage: "
        f"{summary['Acres'].median():.2f}"
    )
    message(
        f"Minimum stand acreage: "
        f"{summary['Acres'].min():.2f}"
    )
    message(
        f"Maximum stand acreage: "
        f"{summary['Acres'].max():.2f}"
    )

    mean_homogeneity = float(
        summary["HOMOGENEITY"].mean()
    )

    area_weighted_homogeneity = float(
        np.average(
            summary["HOMOGENEITY"],
            weights=summary["Acres"],
        )
    )

    message(
        f"Mean stand homogeneity score: "
        f"{mean_homogeneity:.2f} / 100"
    )
    message(
        f"Area-weighted landscape homogeneity score: "
        f"{area_weighted_homogeneity:.2f} / 100"
    )

    message("Mean within-stand standard deviations:")
    for variable in ["STDAGE", "BALIVE", "QMD", "TPA"]:
        sd_field = f"{variable}_SD"

        unweighted_mean_sd = float(
            summary[sd_field].mean()
        )

        area_weighted_mean_sd = float(
            np.average(
                summary[sd_field],
                weights=summary["Acres"],
            )
        )

        message(
            f"  {variable}: "
            f"mean={unweighted_mean_sd:.3f}, "
            f"area-weighted={area_weighted_mean_sd:.3f}"
        )

    message("=" * 60)
    message("WRITING SEGMENT RASTER")
    message("=" * 60)

    write_segment_raster(
        labels,
        geometry,
    )

    message(
        f"Segment raster written:\n"
        f"{segment_raster}"
    )

    message("=" * 60)
    message("CONVERTING SEGMENTS TO POLYGONS")
    message("=" * 60)

    write_segment_polygons(summary)

    message(
        f"Segment polygons written:\n"
        f"{segment_polygons}"
    )

    message("=" * 60)
    message("SEGMENTATION COMPLETE")
    message("=" * 60)


if __name__ == "__main__":
    main()