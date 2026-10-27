"""
02_segment_treemap_v3.py

Boundary-vectorized TreeMap stand segmentation for LETO.

Version 3 preserves synchronous cell reassignment but evaluates only boundary
cells. Segment statistics are computed in linear time with NumPy bincount, and
neighbor candidates are stored as a compact boundary-cell matrix rather than
as several full-raster temporary arrays.
"""

import math
import os
import time
from collections import Counter, defaultdict

import arcpy
import numpy as np
import pandas as pd
from scipy import ndimage


# ============================================================
# USER INPUTS
# ============================================================

input_folder = (
    r"C:\path\to\data"
)

age_raster = os.path.join(input_folder, "TreeMap_STDAGE.tif")
forest_type_raster = os.path.join(input_folder, "TreeMap_FORTYPCD.tif")
basal_area_raster = os.path.join(input_folder, "TreeMap_BALIVE.tif")
qmd_raster = os.path.join(input_folder, "TreeMap_QMD.tif")
tpa_raster = os.path.join(input_folder, "TreeMap_TPA.tif")

output_folder = os.path.join(
    r"C:\path\to\data",
    "TreeMap_Segmentation_Output",
)

output_gdb = os.path.join(output_folder, "TreeMap_Segmentation.gdb")
segment_raster = os.path.join(output_folder, "TreeMap_Segments_Initial.tif")
segment_polygons = os.path.join(output_gdb, "TreeMap_Segments_Initial")
segment_summary_csv = os.path.join(output_folder, "TreeMap_Segment_Summary.csv")

# NHDPlus High Resolution national hydrography inputs.
nhd_flowline = (
    r"C:\path\to\data.gdb\NHDPlus_H_National\NetworkNHDFlowline"
)
nhd_area = (
    r"C:\path\to\data.gdb\NHDPlus_H_National\NHDArea"
)
nhd_waterbody = (
    r"C:\path\to\data.gdb\NHDPlus_H_National\NHDWaterbody"
)

# USFWS National Wetlands Inventory polygon inputs.
# Add one or more state feature classes or shapefiles to this list. The script
# clips, merges, and dissolves them automatically, so state boundaries do not
# affect wetland-percentage calculations.
NWI_WETLAND_INPUTS = [
   r"C:\path\to\data.gdb\FL_Wetlands"
]

# Riparian-processing outputs.
riparian_buffer = os.path.join(output_gdb, "Riparian_Buffer")
management_units_final = os.path.join(
    output_gdb,
    "ManagementUnits_Final",
)
# Backward-compatible internal alias used by the riparian functions.
riparian_split_stands = management_units_final
ownership_dissolved = os.path.join(output_gdb, "Ownership_Classes_Dissolved")
ownership_raster = os.path.join(output_folder, "Ownership_Classes.tif")
wetland_footprint = os.path.join(output_gdb, "NWI_Wetlands_Dissolved")


# ============================================================
# SEGMENTATION SETTINGS
# ============================================================

INITIAL_SEED_ACRES = 100.0
MIN_STAND_ACRES = 5.0
MAX_STAND_ACRES = 300.0

MAX_ITERATIONS = 40
CONVERGENCE_THRESHOLD = 0.01
MINIMUM_SCORE_IMPROVEMENT = 0.01
USE_EIGHT_NEIGHBORS = False

WEIGHTS = {
    "FORTYPCD": 0.30,
    "STDAGE": 0.25,
    "BALIVE": 0.20,
    "QMD": 0.15,
    "TPA": 0.10,
}

SHARED_EDGE_BONUS = 0.1
Z_CLIP = 4.0
READ_NODATA = -999999

# Print a progress message during long vectorized stages.
REPORT_TIMING = True

# Optional post-processing merge for biologically similar adjacent stands.
MERGE_SIMILAR_STANDS = True
SIMILAR_MERGE_MAX_AGE_DIFFERENCE = 10.0  # years
SIMILAR_MERGE_REQUIRE_SAME_FORTYPCD = True
SIMILAR_MERGE_MIN_SHARED_EDGES = 1
SIMILAR_MERGE_MAX_PASSES = 50


# ============================================================
# RIPARIAN BUFFER SETTINGS
# ============================================================

# Set to False to run only the original stand segmentation.
CREATE_RIPARIAN_STANDS = True

# These are editable first-pass assumptions. They are intentionally stored in
# one place so other states or later Florida BMP refinements can use different
# distances without changing the geometry-processing functions.
FLOWLINE_BUFFER_RULES = {
    46000: ("Unclassified Stream", 35.0),
    46003: ("Intermittent Stream", 35.0),
    46006: ("Perennial Stream", 75.0),
    46007: ("Ephemeral Stream", 0.0),
}

WATERBODY_BUFFER_RULES = {
    39000: ("Unclassified Lake/Pond", 35.0),
    39001: ("Intermittent Lake/Pond", 35.0),
    39005: ("Intermittent Lake/Pond", 35.0),
    39006: ("Intermittent Lake/Pond", 35.0),
    39004: ("Perennial Lake/Pond", 75.0),
    39009: ("Perennial Lake/Pond", 75.0),
    39010: ("Perennial Lake/Pond", 75.0),
    39011: ("Perennial Lake/Pond", 75.0),
    39012: ("Perennial Lake/Pond", 75.0),
}

AREA_BUFFER_RULES = {
    46000: ("Unclassified Stream/River Area", 35.0),
    46003: ("Intermittent Stream/River Area", 35.0),
    46006: ("Perennial Stream/River Area", 75.0),
    46007: ("Ephemeral Stream/River Area", 0.0),
    53700: ("Area of Complex Channels", 75.0),
}

# No minimum-area deletion is performed after the riparian split. Narrow
# riparian pieces are intentional management units rather than segmentation
# errors.


# ============================================================
# WETLAND AND MANAGEMENT-CLASS SETTINGS
# ============================================================

# Classify each original TreeMap stand from all polygons in the listed USFWS
# National Wetlands Inventory inputs. No attribute filtering is performed.
# The final riparian-split pieces inherit the original-stand classification.
CLASSIFY_WETLANDS = False
WETLAND_MAJORITY_THRESHOLD = 50.0

# Final management codes. Riparian always has first priority.
UPLAND_MANAGEMENT_CODE = 0
RIPARIAN_MANAGEMENT_CODE = 1
WETLAND_MANAGEMENT_CODE = 2




# ============================================================
# PARCEL OWNERSHIP SETTINGS
# ============================================================

# Update this path to the extracted parcel feature class on your computer.
PARCEL_OWNERSHIP_LAYER = (
    r"C:\path\to\data.gdb\FL_5_Co_Parcels"
)
PARCEL_OWNER_FIELD = "owner"  # Located case-insensitively.
USE_OWNERSHIP_HARD_BOUNDARIES = True

# Ownership is classified by the suffix after the first hyphen. This means
# both Forest-* and Other-* labels enter the same ownership class.
OWNERSHIP_SUFFIX_GROUPS = {
    "private": (1, "Private"),
    "company": (2, "Corporate"),
    "reit": (2, "Corporate"),
    "timo": (2, "Corporate"),
    "federal": (3, "Federal"),
    "state": (4, "State"),
    "county": (5, "County"),
    "ngo": (6, "NGO"),
    "misc": (7, "Other"),
}

UNKNOWN_OWNERSHIP_CODE = 0
UNKNOWN_OWNERSHIP_TYPE = "Unknown"


# ============================================================
# GENERAL SETTINGS
# ============================================================

arcpy.env.overwriteOutput = True
os.makedirs(output_folder, exist_ok=True)

if not arcpy.Exists(output_gdb):
    arcpy.management.CreateFileGDB(output_folder, os.path.basename(output_gdb))


# ============================================================
# HELPERS
# ============================================================

def message(text):
    print(text, flush=True)
    try:
        arcpy.AddMessage(text)
    except Exception:
        pass


def timed_message(label, start_time):
    if REPORT_TIMING:
        elapsed = time.perf_counter() - start_time
        message(f"{label} completed in {elapsed:,.1f} seconds")


def validate_inputs():
    paths = [
        age_raster,
        forest_type_raster,
        basal_area_raster,
        qmd_raster,
        tpa_raster,
    ]
    missing = [path for path in paths if not arcpy.Exists(path)]
    if missing:
        raise FileNotFoundError(
            "The following input rasters were not found:\n" + "\n".join(missing)
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
            math.isclose(reference["cell_width"], current["cell_width"], rel_tol=0, abs_tol=1e-8),
            math.isclose(reference["cell_height"], current["cell_height"], rel_tol=0, abs_tol=1e-8),
            math.isclose(reference["xmin"], current["xmin"], rel_tol=0, abs_tol=1e-6),
            math.isclose(reference["ymin"], current["ymin"], rel_tol=0, abs_tol=1e-6),
        ]
        if not all(checks):
            raise ValueError(f"Raster is not aligned with the reference raster:\n{raster_path}")
    return reference


def read_raster(path, scale=1.0):
    array = arcpy.RasterToNumPyArray(path, nodata_to_value=READ_NODATA).astype(np.float32)
    array[array == READ_NODATA] = np.nan
    if scale != 1.0:
        array /= np.float32(scale)
    return array


def robust_standardize(array, valid_mask):
    values = array[valid_mask]
    median = float(np.nanmedian(values))
    q1 = float(np.nanpercentile(values, 25))
    q3 = float(np.nanpercentile(values, 75))
    iqr = q3 - q1
    scale = float(np.nanstd(values)) if not np.isfinite(iqr) or iqr == 0 else iqr / 1.349
    if not np.isfinite(scale) or scale == 0:
        scale = 1.0
    standardized = ((array - median) / scale).astype(np.float32)
    np.clip(standardized, -Z_CLIP, Z_CLIP, out=standardized)
    return standardized, median, scale


def connectivity_structure():
    return np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=np.int8)


def create_initial_segments(valid_mask, cells_per_seed, ownership):
    rows, cols = valid_mask.shape
    block_side = max(1, int(round(math.sqrt(cells_per_seed))))
    labels = np.zeros((rows, cols), dtype=np.int32)
    next_id = 1
    structure = connectivity_structure()

    for row_start in range(0, rows, block_side):
        row_end = min(rows, row_start + block_side)
        for col_start in range(0, cols, block_side):
            col_end = min(cols, col_start + block_side)
            block_mask = valid_mask[row_start:row_end, col_start:col_end]
            if not np.any(block_mask):
                continue
            block_owner = ownership[row_start:row_end, col_start:col_end]
            target = labels[row_start:row_end, col_start:col_end]
            for owner_code in np.unique(block_owner[block_mask]):
                owner_mask = block_mask & (block_owner == owner_code)
                local_labels, component_count = ndimage.label(
                    owner_mask, structure=structure
                )
                positive = local_labels > 0
                target[positive] = local_labels[positive] + next_id - 1
                next_id += component_count

    return labels


# ============================================================
# FAST SEGMENT STATISTICS
# ============================================================

def segment_statistics_arrays(labels, feature_arrays, forest_type):
    """Return segment statistics as arrays indexed directly by segment ID."""
    valid = labels > 0
    ids = labels[valid].astype(np.int64, copy=False)
    max_id = int(ids.max()) if ids.size else 0

    counts = np.bincount(ids, minlength=max_id + 1).astype(np.int64)
    means = {}

    for name, array in feature_arrays.items():
        values = array[valid].astype(np.float64, copy=False)
        sums = np.bincount(ids, weights=values, minlength=max_id + 1)
        means[name] = np.divide(
            sums,
            counts,
            out=np.zeros(max_id + 1, dtype=np.float64),
            where=counts > 0,
        )

    types = forest_type[valid].astype(np.int64, copy=False)
    unique_types, compact_types = np.unique(types, return_inverse=True)
    n_types = len(unique_types)

    if n_types:
        combined = ids * n_types + compact_types
        pair_counts = np.bincount(
            combined,
            minlength=(max_id + 1) * n_types,
        ).reshape(max_id + 1, n_types)
        dominant_index = np.argmax(pair_counts, axis=1)
        dominant_type = unique_types[dominant_index].astype(np.int32)
        dominant_count = pair_counts[np.arange(max_id + 1), dominant_index].astype(np.int64)
    else:
        dominant_type = np.zeros(max_id + 1, dtype=np.int32)
        dominant_count = np.zeros(max_id + 1, dtype=np.int64)

    return {
        "counts": counts,
        "means": means,
        "dominant_type": dominant_type,
        "dominant_count": dominant_count,
    }


def identify_boundary_cells(labels, valid_mask):
    boundary = np.zeros(labels.shape, dtype=bool)
    boundary[:-1, :] |= valid_mask[:-1, :] & valid_mask[1:, :] & (labels[:-1, :] != labels[1:, :])
    boundary[1:, :] |= valid_mask[1:, :] & valid_mask[:-1, :] & (labels[1:, :] != labels[:-1, :])
    boundary[:, :-1] |= valid_mask[:, :-1] & valid_mask[:, 1:] & (labels[:, :-1] != labels[:, 1:])
    boundary[:, 1:] |= valid_mask[:, 1:] & valid_mask[:, :-1] & (labels[:, 1:] != labels[:, :-1])
    return boundary


def boundary_neighbor_matrix(labels, boundary_rows, boundary_cols):
    """
    Return adjacent segment IDs for boundary cells only.

    Shape is (number of boundary cells, number of neighbor directions). This
    avoids constructing four or eight complete shifted rasters each iteration.
    """
    rows, cols = labels.shape

    offsets = [
        (-1, 0),
        (0, -1),
        (0, 1),
        (1, 0),
    ]

    if USE_EIGHT_NEIGHBORS:
        offsets = [
            (-1, -1), (-1, 0), (-1, 1),
            (0, -1),            (0, 1),
            (1, -1),  (1, 0),   (1, 1),
        ]

    neighbors = np.zeros(
        (boundary_rows.size, len(offsets)),
        dtype=np.int32,
    )

    for index, (dr, dc) in enumerate(offsets):
        rr = boundary_rows + dr
        cc = boundary_cols + dc
        inside = (
            (rr >= 0)
            & (rr < rows)
            & (cc >= 0)
            & (cc < cols)
        )
        neighbors[inside, index] = labels[rr[inside], cc[inside]]

    return neighbors


def assignment_cost_vector(
    candidate_ids,
    edge_counts,
    stats,
    cell_features,
    cell_forest_type,
):
    """Calculate assignment cost for a vector of boundary-cell candidates."""
    cost = np.zeros(candidate_ids.size, dtype=np.float64)
    means = stats["means"]

    for name in ["STDAGE", "BALIVE", "QMD", "TPA"]:
        difference = cell_features[name] - means[name][candidate_ids]
        cost += WEIGHTS[name] * difference * difference

    cost += np.where(
        stats["dominant_type"][candidate_ids] == cell_forest_type,
        0.0,
        WEIGHTS["FORTYPCD"],
    )
    cost -= SHARED_EDGE_BONUS * edge_counts

    return cost


def segment_categorical_mode(labels, values):
    """Return the dominant categorical value for each segment ID."""
    valid = labels > 0
    ids = labels[valid].astype(np.int64, copy=False)
    max_id = int(ids.max()) if ids.size else 0
    categories = values[valid].astype(np.int64, copy=False)
    unique_values, compact = np.unique(categories, return_inverse=True)
    n_values = len(unique_values)
    if n_values == 0:
        return np.zeros(max_id + 1, dtype=np.int32)
    combined = ids * n_values + compact
    counts = np.bincount(
        combined, minlength=(max_id + 1) * n_values
    ).reshape(max_id + 1, n_values)
    return unique_values[np.argmax(counts, axis=1)].astype(np.int32)


def iterative_reassignment(labels, valid_mask, feature_arrays, forest_type, ownership):
    """
    Synchronously reassign all boundary cells using compact vector arrays.

    The Python loop is over neighbor directions (four or eight), not over
    individual cells. Every candidate direction is evaluated for all boundary
    cells simultaneously using NumPy array operations.
    """
    total_valid_cells = int(valid_mask.sum())
    forest_type_int = np.where(valid_mask, forest_type, 0).astype(np.int32)

    for iteration in range(1, MAX_ITERATIONS + 1):
        iteration_start = time.perf_counter()

        # O(number of valid cells), rather than one raster scan per segment.
        stats = segment_statistics_arrays(
            labels,
            feature_arrays,
            forest_type_int,
        )
        segment_owner = segment_categorical_mode(labels, ownership)

        boundary = identify_boundary_cells(labels, valid_mask)
        boundary_rows, boundary_cols = np.where(boundary)
        boundary_count = boundary_rows.size

        if boundary_count == 0:
            message("No segment boundaries remain.")
            break

        current_ids = labels[boundary_rows, boundary_cols]
        neighbor_ids = boundary_neighbor_matrix(
            labels,
            boundary_rows,
            boundary_cols,
        )

        # Values are extracted only for boundary cells.
        cell_features = {
            name: array[boundary_rows, boundary_cols].astype(
                np.float64,
                copy=False,
            )
            for name, array in feature_arrays.items()
        }
        cell_forest_type = forest_type_int[
            boundary_rows,
            boundary_cols,
        ]
        cell_owner = ownership[boundary_rows, boundary_cols]

        current_edge_counts = np.sum(
            neighbor_ids == current_ids[:, None],
            axis=1,
        ).astype(np.float64)

        best_ids = current_ids.copy()
        best_cost = assignment_cost_vector(
            current_ids,
            current_edge_counts,
            stats,
            cell_features,
            cell_forest_type,
        )

        # Only four iterations for rook adjacency, regardless of the number
        # of boundary cells.
        for direction in range(neighbor_ids.shape[1]):
            candidate_ids = neighbor_ids[:, direction]
            candidate_valid = (
                (candidate_ids > 0)
                & (candidate_ids != current_ids)
                & (segment_owner[candidate_ids] == cell_owner)
            )

            if not np.any(candidate_valid):
                continue

            candidate_edge_counts = np.sum(
                neighbor_ids == candidate_ids[:, None],
                axis=1,
            ).astype(np.float64)

            candidate_cost = assignment_cost_vector(
                candidate_ids,
                candidate_edge_counts,
                stats,
                cell_features,
                cell_forest_type,
            )

            improve = candidate_valid & (
                candidate_cost
                < best_cost - MINIMUM_SCORE_IMPROVEMENT
            )

            best_cost[improve] = candidate_cost[improve]
            best_ids[improve] = candidate_ids[improve]

        changed_mask = best_ids != current_ids
        changed = int(changed_mask.sum())

        # Synchronous update: all decisions used the same original labels.
        if changed:
            new_labels = labels.copy()
            new_labels[
                boundary_rows[changed_mask],
                boundary_cols[changed_mask],
            ] = best_ids[changed_mask]
            labels = new_labels

        changed_fraction = (
            changed / total_valid_cells
            if total_valid_cells
            else 0.0
        )
        elapsed = time.perf_counter() - iteration_start

        message(
            f"Iteration {iteration:02d}: "
            f"{boundary_count:,} boundary cells, "
            f"{changed:,} changed "
            f"({changed_fraction:.4%}), "
            f"{elapsed:,.1f} seconds"
        )

        if changed_fraction <= CONVERGENCE_THRESHOLD:
            message("Convergence threshold reached.")
            break

    return labels


# ============================================================
# CONNECTIVITY AND ID MANAGEMENT
# ============================================================

def split_disconnected_segments(labels):
    """Split disconnected pieces using per-segment bounding boxes, not full rasters."""
    output = np.zeros_like(labels, dtype=np.int32)
    max_id = int(labels.max())
    if max_id == 0:
        return output

    slices = ndimage.find_objects(labels, max_label=max_id)
    structure = connectivity_structure()
    next_id = 1

    for old_id, bbox in enumerate(slices, start=1):
        if bbox is None:
            continue
        local = labels[bbox]
        local_mask = local == old_id
        pieces, piece_count = ndimage.label(local_mask, structure=structure)
        target = output[bbox]
        positive = pieces > 0
        target[positive] = pieces[positive] + next_id - 1
        next_id += piece_count

    return output


def renumber_segments(labels):
    max_id = int(labels.max())
    if max_id == 0:
        return labels.astype(np.int32, copy=True)
    present = np.bincount(labels.ravel(), minlength=max_id + 1) > 0
    present[0] = False
    lookup = np.zeros(max_id + 1, dtype=np.int32)
    lookup[present] = np.arange(1, int(present.sum()) + 1, dtype=np.int32)
    return lookup[labels]


def segment_adjacency_arrays(labels):
    """Return unique adjacent segment pairs and their shared-edge counts."""
    pair_a = []
    pair_b = []

    left = labels[:, :-1]
    right = labels[:, 1:]
    mask = (left > 0) & (right > 0) & (left != right)
    pair_a.append(left[mask].astype(np.int64, copy=False))
    pair_b.append(right[mask].astype(np.int64, copy=False))

    upper = labels[:-1, :]
    lower = labels[1:, :]
    mask = (upper > 0) & (lower > 0) & (upper != lower)
    pair_a.append(upper[mask].astype(np.int64, copy=False))
    pair_b.append(lower[mask].astype(np.int64, copy=False))

    if not pair_a or sum(arr.size for arr in pair_a) == 0:
        return np.empty((0, 2), dtype=np.int64), np.empty(0, dtype=np.int64)

    a = np.concatenate(pair_a)
    b = np.concatenate(pair_b)
    lo = np.minimum(a, b)
    hi = np.maximum(a, b)
    max_id = int(labels.max()) + 1
    encoded = lo * max_id + hi
    unique, counts = np.unique(encoded, return_counts=True)
    pairs = np.column_stack((unique // max_id, unique % max_id)).astype(np.int64)
    return pairs, counts.astype(np.int64)


# ============================================================
# OVERSIZED AND UNDERSIZED SEGMENTS
# ============================================================

def _connected_components(mask):
    return ndimage.label(mask, structure=connectivity_structure())


def _bisect_connected_mask(mask):
    rows, cols = np.where(mask)
    if len(rows) < 2:
        return [mask]

    row_span = int(rows.max() - rows.min())
    col_span = int(cols.max() - cols.min())

    if col_span >= row_span:
        cut = int(np.median(cols))
        grid = np.indices(mask.shape)[1]
        first = mask & (grid <= cut)
        second = mask & (grid > cut)
    else:
        cut = int(np.median(rows))
        grid = np.indices(mask.shape)[0]
        first = mask & (grid <= cut)
        second = mask & (grid > cut)

    pieces = []
    for part in (first, second):
        if not np.any(part):
            continue
        component_labels, component_count = _connected_components(part)
        for component_id in range(1, component_count + 1):
            component = component_labels == component_id
            if np.any(component):
                pieces.append(component)

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
    output = np.zeros_like(labels, dtype=np.int32)
    max_id = int(labels.max())
    slices = ndimage.find_objects(labels, max_label=max_id)
    next_id = 1
    split_count = 0

    for segment_id, bbox in enumerate(slices, start=1):
        if bbox is None:
            continue
        local_labels = labels[bbox]
        queue = [local_labels == segment_id]
        local_output = output[bbox]

        while queue:
            mask = queue.pop()
            cell_count = int(mask.sum())
            if cell_count <= maximum_cells or cell_count <= 1:
                local_output[mask] = next_id
                next_id += 1
                continue
            pieces = _bisect_connected_mask(mask)
            if len(pieces) < 2:
                local_output[mask] = next_id
                next_id += 1
                continue
            split_count += 1
            queue.extend(pieces)

    message(f"Oversized-segment split operations: {split_count:,}")
    return output


def merge_small_segments(labels, feature_arrays, forest_type, ownership, minimum_cells, maximum_cells=None):
    """Merge many small segments per pass and remap the raster only once per pass."""
    pass_number = 0

    while True:
        pass_number += 1
        start = time.perf_counter()
        stats = segment_statistics_arrays(labels, feature_arrays, forest_type)
        segment_owner = segment_categorical_mode(labels, ownership)
        counts = stats["counts"].copy()
        active_ids = np.flatnonzero(counts > 0)
        active_ids = active_ids[active_ids > 0]
        small_ids = active_ids[counts[active_ids] < minimum_cells]

        if small_ids.size == 0:
            break

        pairs, shared_counts = segment_adjacency_arrays(labels)
        adjacency = defaultdict(list)
        for (a, b), shared in zip(pairs, shared_counts):
            adjacency[int(a)].append((int(b), int(shared)))
            adjacency[int(b)].append((int(a), int(shared)))

        parent = np.arange(len(counts), dtype=np.int64)
        current_size = counts.astype(np.int64, copy=True)

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return int(x)

        merged = 0
        order = small_ids[np.argsort(counts[small_ids])]

        for small_id in order:
            source_root = find(int(small_id))
            if source_root != int(small_id) or current_size[source_root] >= minimum_cells:
                continue

            best_root = None
            best_cost = float("inf")

            for neighbor_id, shared_edges in adjacency.get(int(small_id), []):
                target_root = find(neighbor_id)
                if target_root == source_root or current_size[target_root] == 0:
                    continue
                if segment_owner[source_root] != segment_owner[target_root]:
                    continue
                combined_size = current_size[source_root] + current_size[target_root]
                if maximum_cells is not None and combined_size > maximum_cells:
                    continue

                cost = 0.0
                for name in ["STDAGE", "BALIVE", "QMD", "TPA"]:
                    difference = (
                        stats["means"][name][source_root]
                        - stats["means"][name][target_root]
                    )
                    cost += WEIGHTS[name] * difference * difference

                if stats["dominant_type"][source_root] != stats["dominant_type"][target_root]:
                    cost += WEIGHTS["FORTYPCD"]
                cost -= SHARED_EDGE_BONUS * shared_edges

                if cost < best_cost:
                    best_cost = cost
                    best_root = target_root

            if best_root is None:
                continue

            # Update the continuous sufficient statistics immediately so
            # later decisions in this pass use the combined group means.
            source_size = current_size[source_root]
            target_size = current_size[best_root]
            combined_size = source_size + target_size

            for name in ["STDAGE", "BALIVE", "QMD", "TPA"]:
                combined_sum = (
                    stats["means"][name][source_root] * source_size
                    + stats["means"][name][best_root] * target_size
                )
                stats["means"][name][best_root] = (
                    combined_sum / combined_size
                    if combined_size > 0
                    else 0.0
                )

            # Dominant forest type is updated conservatively from the
            # dominant-type counts available for each group.
            if (
                stats["dominant_type"][source_root]
                == stats["dominant_type"][best_root]
            ):
                stats["dominant_count"][best_root] += (
                    stats["dominant_count"][source_root]
                )
            elif (
                stats["dominant_count"][source_root]
                > stats["dominant_count"][best_root]
            ):
                stats["dominant_type"][best_root] = (
                    stats["dominant_type"][source_root]
                )
                stats["dominant_count"][best_root] = (
                    stats["dominant_count"][source_root]
                )

            parent[source_root] = best_root
            current_size[best_root] = combined_size
            current_size[source_root] = 0
            merged += 1

        message(
            f"Small-segment merge pass {pass_number}: {merged:,} segments merged "
            f"in {time.perf_counter() - start:,.1f} seconds"
        )

        if merged == 0:
            break

        lookup = np.arange(len(parent), dtype=np.int64)
        for segment_id in active_ids:
            lookup[segment_id] = find(int(segment_id))
        labels = lookup[labels].astype(np.int32)
        labels = split_disconnected_segments(labels)

    return labels



def merge_similar_adjacent_stands(
    labels,
    raw_age,
    forest_type,
    ownership,
    maximum_cells,
    maximum_age_difference,
    require_same_forest_type=True,
    minimum_shared_edges=1,
    maximum_passes=50,
):
    """
    Greedily merge adjacent stands that are biologically similar.

    A candidate pair must remain under the maximum stand size, meet the age
    threshold, and (optionally) share the same dominant FORTYPCD. Each stand
    can participate in at most one merge per pass so all decisions in a pass
    use current, internally consistent stand statistics.
    """
    total_merged = 0

    for pass_number in range(1, maximum_passes + 1):
        start = time.perf_counter()

        stats = segment_statistics_arrays(
            labels,
            {"STDAGE": raw_age},
            forest_type,
        )
        counts = stats["counts"]
        age_means = stats["means"]["STDAGE"]
        dominant_type = stats["dominant_type"]
        segment_owner = segment_categorical_mode(labels, ownership)

        pairs, shared_edges = segment_adjacency_arrays(labels)
        if pairs.size == 0:
            break

        first = pairs[:, 0]
        second = pairs[:, 1]
        combined = counts[first] + counts[second]
        age_difference = np.abs(age_means[first] - age_means[second])

        eligible = (
            (combined <= maximum_cells)
            & (age_difference <= maximum_age_difference)
            & (shared_edges >= minimum_shared_edges)
            & (segment_owner[first] == segment_owner[second])
        )

        if require_same_forest_type:
            eligible &= dominant_type[first] == dominant_type[second]

        candidate_indices = np.flatnonzero(eligible)
        if candidate_indices.size == 0:
            message(
                f"Similar-stand merge pass {pass_number}: 0 segments merged "
                f"in {time.perf_counter() - start:,.1f} seconds"
            )
            break

        # Prefer the closest ages, then the longest shared boundary.
        order = candidate_indices[
            np.lexsort((
                -shared_edges[candidate_indices],
                age_difference[candidate_indices],
            ))
        ]

        used = np.zeros(len(counts), dtype=bool)
        parent = np.arange(len(counts), dtype=np.int64)
        merged_this_pass = 0

        for index in order:
            a = int(first[index])
            b = int(second[index])
            if used[a] or used[b]:
                continue

            # Keep the larger segment ID as the target only for deterministic
            # behavior; IDs are renumbered later.
            if counts[a] > counts[b] or (counts[a] == counts[b] and a < b):
                target, source = a, b
            else:
                target, source = b, a

            parent[source] = target
            used[source] = True
            used[target] = True
            merged_this_pass += 1

        message(
            f"Similar-stand merge pass {pass_number}: "
            f"{merged_this_pass:,} segments merged "
            f"in {time.perf_counter() - start:,.1f} seconds"
        )

        if merged_this_pass == 0:
            break

        labels = parent[labels].astype(np.int32)
        labels = split_disconnected_segments(labels)
        total_merged += merged_this_pass

    message(f"Total similar adjacent stands merged: {total_merged:,}")
    return labels


def diagnose_remaining_small_segments(labels, minimum_cells, maximum_cells):
    """Explain why undersized segments could not be merged."""
    counts = np.bincount(labels.ravel())
    small_ids = np.flatnonzero(
        (counts > 0) & (counts < minimum_cells)
    )
    small_ids = small_ids[small_ids > 0]

    pairs, _ = segment_adjacency_arrays(labels)
    adjacency = defaultdict(list)
    for a, b in pairs:
        adjacency[int(a)].append(int(b))
        adjacency[int(b)].append(int(a))

    no_neighbor = 0
    maximum_blocked = 0
    otherwise_unmerged = 0

    for segment_id in small_ids:
        neighbors = adjacency.get(int(segment_id), [])
        if not neighbors:
            no_neighbor += 1
            continue

        has_size_eligible_neighbor = any(
            counts[int(segment_id)] + counts[neighbor] <= maximum_cells
            for neighbor in neighbors
        )

        if not has_size_eligible_neighbor:
            maximum_blocked += 1
        else:
            otherwise_unmerged += 1

    message("Remaining undersized-segment diagnostic:")
    message(f"  No adjacent segment: {no_neighbor:,}")
    message(
        f"  All adjacent merges exceed maximum size: "
        f"{maximum_blocked:,}"
    )
    message(
        f"  At least one size-eligible neighbor remained: "
        f"{otherwise_unmerged:,}"
    )

    return {
        "no_neighbor": no_neighbor,
        "maximum_blocked": maximum_blocked,
        "otherwise_unmerged": otherwise_unmerged,
    }

def delete_small_segments(labels, minimum_cells):
    counts = np.bincount(labels.ravel())
    small = np.flatnonzero((counts > 0) & (counts < minimum_cells))
    small = small[small > 0]
    delete_lookup = np.zeros(len(counts), dtype=bool)
    delete_lookup[small] = True
    delete_mask = delete_lookup[labels]
    output = labels.copy()
    output[delete_mask] = 0
    return output, int(len(small)), int(delete_mask.sum())


# ============================================================
# FAST SUMMARY
# ============================================================

def create_summary(labels, raw_arrays, forest_type, ownership, cell_acres, homogeneity_scales):
    valid = labels > 0
    ids = labels[valid].astype(np.int64, copy=False)
    max_id = int(ids.max()) if ids.size else 0
    counts = np.bincount(ids, minlength=max_id + 1).astype(np.int64)
    segment_ids = np.flatnonzero(counts > 0)
    segment_ids = segment_ids[segment_ids > 0]

    data = {
        "SEG_ID": segment_ids.astype(np.int32),
        "PIXEL_COUNT": counts[segment_ids],
        "Acres": counts[segment_ids] * cell_acres,
    }
    component_scores = {}

    for name, array in raw_arrays.items():
        values = array[valid].astype(np.float64, copy=False)
        sums = np.bincount(ids, weights=values, minlength=max_id + 1)
        sums_sq = np.bincount(ids, weights=values * values, minlength=max_id + 1)
        means = np.divide(sums, counts, out=np.zeros(max_id + 1), where=counts > 0)
        variances = np.divide(sums_sq, counts, out=np.zeros(max_id + 1), where=counts > 0) - means * means
        variances = np.maximum(variances, 0.0)
        sds = np.sqrt(variances)

        mins = ndimage.minimum(array, labels=labels, index=segment_ids)
        maxs = ndimage.maximum(array, labels=labels, index=segment_ids)

        data[f"{name}_MEAN"] = means[segment_ids]
        data[f"{name}_SD"] = sds[segment_ids]
        data[f"{name}_MIN"] = mins
        data[f"{name}_MAX"] = maxs

        landscape_sd = float(homogeneity_scales[name])
        if not np.isfinite(landscape_sd) or landscape_sd <= 0:
            component = np.where(sds[segment_ids] == 0, 1.0, 0.0)
        else:
            component = np.clip(1.0 - sds[segment_ids] / landscape_sd, 0.0, 1.0)
        component_scores[name] = component
        data[f"{name}_HOM"] = 100.0 * component

    types = forest_type[valid].astype(np.int64, copy=False)
    unique_types, compact_types = np.unique(types, return_inverse=True)
    n_types = len(unique_types)
    combined = ids * n_types + compact_types
    pair_counts = np.bincount(combined, minlength=(max_id + 1) * n_types).reshape(max_id + 1, n_types)
    dominant_index = np.argmax(pair_counts, axis=1)
    dominant_type = unique_types[dominant_index]
    dominant_count = pair_counts[np.arange(max_id + 1), dominant_index]
    forest_fraction = dominant_count[segment_ids] / counts[segment_ids]

    data["FORTYPCD_DOM"] = dominant_type[segment_ids].astype(np.int32)
    owner_mode = segment_categorical_mode(labels, ownership)
    data["OWN_CODE"] = owner_mode[segment_ids].astype(np.int16)
    code_to_type = {code: name for code, name in OWNERSHIP_SUFFIX_GROUPS.values()}
    code_to_type[UNKNOWN_OWNERSHIP_CODE] = UNKNOWN_OWNERSHIP_TYPE
    data["OWN_TYPE"] = [code_to_type.get(int(code), UNKNOWN_OWNERSHIP_TYPE) for code in data["OWN_CODE"]]
    data["FORTYPE_PCT"] = 100.0 * forest_fraction
    data["FORTYPCD_HOM"] = 100.0 * forest_fraction

    weighted_score = (
        WEIGHTS["STDAGE"] * component_scores["STDAGE"]
        + WEIGHTS["BALIVE"] * component_scores["BALIVE"]
        + WEIGHTS["QMD"] * component_scores["QMD"]
        + WEIGHTS["TPA"] * component_scores["TPA"]
        + WEIGHTS["FORTYPCD"] * forest_fraction
    )
    data["HOMOGENEITY"] = 100.0 * weighted_score
    return pd.DataFrame(data)


# ============================================================
# OUTPUTS
# ============================================================

def write_segment_raster(labels, geometry):
    output_nodata = -9999
    export_array = labels.astype(np.int32)
    export_array[export_array == 0] = output_nodata
    lower_left = arcpy.Point(geometry["xmin"], geometry["ymin"])
    raster = arcpy.NumPyArrayToRaster(
        export_array,
        lower_left,
        geometry["cell_width"],
        geometry["cell_height"],
        value_to_nodata=output_nodata,
    )
    raster.save(segment_raster)
    arcpy.management.DefineProjection(segment_raster, geometry["spatial_reference"])
    arcpy.management.CalculateStatistics(segment_raster)
    try:
        arcpy.management.BuildRasterAttributeTable(segment_raster, "Overwrite")
    except arcpy.ExecuteError:
        message("Raster attribute table was not built, but the segment raster was created.")


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

    existing_fields = {field.name.upper() for field in arcpy.ListFields(segment_polygons)}
    if "SEG_ID" not in existing_fields:
        arcpy.management.AddField(segment_polygons, "SEG_ID", "LONG")

    value_field = None
    for field in arcpy.ListFields(segment_polygons):
        if field.name.upper() in {"GRIDCODE", "GRID_CODE"}:
            value_field = field.name
            break
    if value_field is None:
        raise RuntimeError("Could not identify the raster-value field in the polygon output.")

    arcpy.management.CalculateField(segment_polygons, "SEG_ID", f"!{value_field}!", "PYTHON3")

    summary_table = os.path.join(output_gdb, "TreeMap_Segment_Summary")
    if arcpy.Exists(summary_table):
        arcpy.management.Delete(summary_table)

    arcpy.conversion.TableToTable(segment_summary_csv, output_gdb, os.path.basename(summary_table))
    join_fields = [field for field in summary.columns if field != "SEG_ID"]
    arcpy.management.JoinField(segment_polygons, "SEG_ID", summary_table, "SEG_ID", join_fields)


# ============================================================
# RIPARIAN BUFFERING AND STAND SPLITTING
# ============================================================

def _field_name(feature_class, requested_name):
    """Return the actual field name using a case-insensitive match."""
    requested_upper = requested_name.upper()
    for field in arcpy.ListFields(feature_class):
        if field.name.upper() == requested_upper:
            return field.name
    raise RuntimeError(
        f"Field '{requested_name}' was not found in:\n{feature_class}"
    )


def _delete_if_exists(*paths):
    for path in paths:
        if arcpy.Exists(path):
            arcpy.management.Delete(path)


def _add_field_if_missing(feature_class, name, field_type, **kwargs):
    existing = {field.name.upper() for field in arcpy.ListFields(feature_class)}
    if name.upper() not in existing:
        arcpy.management.AddField(feature_class, name, field_type, **kwargs)


def _maximum_configured_buffer():
    distances = []
    for rules in (
        FLOWLINE_BUFFER_RULES,
        WATERBODY_BUFFER_RULES,
        AREA_BUFFER_RULES,
    ):
        distances.extend(float(value[1]) for value in rules.values())
    return max(distances, default=0.0)


def _classify_hydrography(
    clipped_features,
    output_features,
    source_name,
    rules,
):
    """
    Retain supported FCodes and attach the editable riparian-buffer rule.

    Features assigned a zero-foot distance are removed from the buffering
    workflow. Their rules remain visible in the configuration dictionaries so
    they can be enabled later by changing only the distance.
    """
    _delete_if_exists(output_features)
    arcpy.management.CopyFeatures(clipped_features, output_features)

    fcode_field = _field_name(output_features, "FCode")
    _add_field_if_missing(output_features, "HYDRO_SRC", "TEXT", field_length=20)
    _add_field_if_missing(output_features, "HYDRO_CLASS", "TEXT", field_length=80)
    _add_field_if_missing(output_features, "BUFFER_FT", "DOUBLE")
    _add_field_if_missing(output_features, "BUFFER_TXT", "TEXT", field_length=30)
    _add_field_if_missing(output_features, "RULE_SOURCE", "TEXT", field_length=80)

    kept = 0
    removed = 0
    with arcpy.da.UpdateCursor(
        output_features,
        [
            fcode_field,
            "HYDRO_SRC",
            "HYDRO_CLASS",
            "BUFFER_FT",
            "BUFFER_TXT",
            "RULE_SOURCE",
        ],
    ) as cursor:
        for row in cursor:
            try:
                fcode = int(row[0])
            except (TypeError, ValueError):
                cursor.deleteRow()
                removed += 1
                continue

            rule = rules.get(fcode)
            if rule is None or float(rule[1]) <= 0:
                cursor.deleteRow()
                removed += 1
                continue

            hydro_class, distance_ft = rule
            row[1] = source_name
            row[2] = hydro_class
            row[3] = float(distance_ft)
            row[4] = f"{float(distance_ft):g} Feet"
            row[5] = "Editable NHD FCode assumption"
            cursor.updateRow(row)
            kept += 1

    message(
        f"{source_name}: retained {kept:,} bufferable features; "
        f"excluded {removed:,} unsupported or zero-buffer features"
    )
    return kept


def _classify_original_stands_as_wetland(stand_aoi, scratch):
    """
    Calculate NWI wetland coverage for each original stand.

    Every polygon from every feature class listed in NWI_WETLAND_INPUTS is
    treated as wetland. Inputs are clipped to the stand footprint, merged, and
    dissolved before intersection so overlapping state files or polygons are
    never counted twice. Wetland status is assigned to the original TreeMap
    stand using a majority threshold, and all later riparian-split pieces
    inherit that classification through SEG_ID.
    """
    for field_name, field_type in (
        ("WETLAND_AC", "DOUBLE"),
        ("WETLAND_PCT", "DOUBLE"),
        ("WETLAND", "SHORT"),
    ):
        _add_field_if_missing(segment_polygons, field_name, field_type)

    # Reset fields so rerunning the script cannot retain stale classifications.
    with arcpy.da.UpdateCursor(
        segment_polygons,
        ["WETLAND_AC", "WETLAND_PCT", "WETLAND"],
    ) as cursor:
        for row in cursor:
            row[0] = 0.0
            row[1] = 0.0
            row[2] = 0
            cursor.updateRow(row)

    if not CLASSIFY_WETLANDS:
        message("Wetland classification is disabled.")
        _delete_if_exists(wetland_footprint)
        return

    if not NWI_WETLAND_INPUTS:
        raise ValueError(
            "CLASSIFY_WETLANDS is True, but NWI_WETLAND_INPUTS is empty. "
            "Add at least one NWI polygon feature class or shapefile."
        )

    missing_inputs = [
        path for path in NWI_WETLAND_INPUTS if not arcpy.Exists(path)
    ]
    if missing_inputs:
        raise FileNotFoundError(
            "The following NWI wetland inputs were not found:\n"
            + "\n".join(missing_inputs)
        )

    wetland_merged = scratch["wetland_merged"]
    wetland_intersection = scratch["wetland_intersection"]
    clip_outputs = []

    _delete_if_exists(wetland_merged, wetland_intersection, wetland_footprint)

    total_source_features = 0
    for index, wetland_input in enumerate(NWI_WETLAND_INPUTS, start=1):
        clip_output = os.path.join(output_gdb, f"_nwi_clip_{index:03d}")
        _delete_if_exists(clip_output)
        arcpy.analysis.PairwiseClip(wetland_input, stand_aoi, clip_output)

        clipped_count = int(arcpy.management.GetCount(clip_output)[0])
        total_source_features += clipped_count
        message(
            f"NWI input {index}: retained {clipped_count:,} wetland features "
            f"within the stand footprint"
        )
        if clipped_count > 0:
            clip_outputs.append(clip_output)
        else:
            _delete_if_exists(clip_output)

    if not clip_outputs:
        message("No NWI wetland polygons occurred within the stand area.")
        return

    if len(clip_outputs) == 1:
        arcpy.management.CopyFeatures(clip_outputs[0], wetland_merged)
    else:
        arcpy.management.Merge(clip_outputs, wetland_merged)

    # Dissolving removes overlap within and between state NWI inputs.
    arcpy.analysis.PairwiseDissolve(wetland_merged, wetland_footprint)
    arcpy.analysis.PairwiseIntersect(
        [segment_polygons, wetland_footprint],
        wetland_intersection,
        join_attributes="ALL",
    )

    seg_id_field = _field_name(wetland_intersection, "SEG_ID")
    wetland_acres_by_segment = defaultdict(float)
    with arcpy.da.SearchCursor(
        wetland_intersection,
        [seg_id_field, "SHAPE@"],
    ) as cursor:
        for segment_id, geometry in cursor:
            if segment_id is None or geometry is None:
                continue
            wetland_acres_by_segment[int(segment_id)] += float(
                geometry.getArea("GEODESIC", "ACRES")
            )

    original_seg_id = _field_name(segment_polygons, "SEG_ID")
    wetland_stands = 0
    wetland_stand_acres = 0.0
    total_wetland_overlap = 0.0

    with arcpy.da.UpdateCursor(
        segment_polygons,
        [
            original_seg_id,
            "SHAPE@",
            "WETLAND_AC",
            "WETLAND_PCT",
            "WETLAND",
        ],
    ) as cursor:
        for row in cursor:
            segment_id = int(row[0])
            stand_acres = float(row[1].getArea("GEODESIC", "ACRES"))
            wetland_acres = min(
                wetland_acres_by_segment.get(segment_id, 0.0),
                stand_acres,
            )
            wetland_pct = (
                100.0 * wetland_acres / stand_acres
                if stand_acres > 0
                else 0.0
            )
            is_wetland = int(wetland_pct >= WETLAND_MAJORITY_THRESHOLD)

            row[2] = wetland_acres
            row[3] = wetland_pct
            row[4] = is_wetland
            cursor.updateRow(row)

            total_wetland_overlap += wetland_acres
            if is_wetland:
                wetland_stands += 1
                wetland_stand_acres += stand_acres

    message(
        f"NWI wetland features used after clipping: "
        f"{total_source_features:,}; overlap with stands: "
        f"{total_wetland_overlap:,.2f} acres"
    )
    message(
        f"Original stands classified as wetland at "
        f">= {WETLAND_MAJORITY_THRESHOLD:g}% coverage: "
        f"{wetland_stands:,} ({wetland_stand_acres:,.2f} stand acres)"
    )

    # Intermediate state clips are not needed after the dissolved QA layer is
    # created. The dissolved NWI_Wetlands_Dissolved output is retained.
    _delete_if_exists(*clip_outputs, wetland_merged, wetland_intersection)


def _assign_management_classes(feature_class):
    """Assign final management class with Riparian > Wetland > Upland."""
    _add_field_if_missing(feature_class, "MGMT_CLASS", "SHORT")
    _add_field_if_missing(feature_class, "MGMT_TYPE", "TEXT", field_length=20)

    fields = ["RIPARIAN", "WETLAND", "MGMT_CLASS", "MGMT_TYPE"]
    counts = Counter()
    acres = defaultdict(float)

    with arcpy.da.UpdateCursor(feature_class, fields + ["Acres"]) as cursor:
        for row in cursor:
            riparian = int(row[0] or 0)
            wetland = int(row[1] or 0)

            if riparian == 1:
                code = RIPARIAN_MANAGEMENT_CODE
                label = "Riparian"
            elif wetland == 1:
                code = WETLAND_MANAGEMENT_CODE
                label = "Wetland"
            else:
                code = UPLAND_MANAGEMENT_CODE
                label = "Upland"

            row[2] = code
            row[3] = label
            cursor.updateRow(row)
            counts[label] += 1
            acres[label] += float(row[4] or 0.0)

    message("Final management classes:")
    for label in ("Upland", "Riparian", "Wetland"):
        message(
            f"  {label}: {counts[label]:,} pieces "
            f"({acres[label]:,.2f} acres)"
        )


def _calculate_split_fields(feature_class):
    """Assign stable parent labels and recalculate acreage from geometry."""
    seg_id_field = _field_name(feature_class, "SEG_ID")
    _add_field_if_missing(feature_class, "FINAL_ID", "LONG")
    _add_field_if_missing(feature_class, "MU_ID", "LONG")
    _add_field_if_missing(feature_class, "ParentID", "TEXT", field_length=40)
    _add_field_if_missing(feature_class, "Acres", "DOUBLE")

    with arcpy.da.UpdateCursor(
        feature_class,
        ["OID@", seg_id_field, "FINAL_ID", "MU_ID", "ParentID", "Acres", "SHAPE@"],
    ) as cursor:
        for row in cursor:
            row[2] = int(row[0])
            row[3] = int(row[0])
            row[4] = f"SEG_{int(row[1])}"
            row[5] = float(row[6].getArea("GEODESIC", "ACRES"))
            cursor.updateRow(row)



def _classify_owner_value(raw_owner):
    """Map Forest-* and Other-* labels to the same ownership class."""
    text = "" if raw_owner is None else str(raw_owner).strip()
    suffix = text.casefold().split("-", 1)[-1] if "-" in text else text.casefold()
    return OWNERSHIP_SUFFIX_GROUPS.get(
        suffix, (UNKNOWN_OWNERSHIP_CODE, UNKNOWN_OWNERSHIP_TYPE)
    )


def prepare_ownership_boundaries():
    """Dissolve parcel ownership classes and rasterize them to the model grid."""
    if not USE_OWNERSHIP_HARD_BOUNDARIES:
        message("Ownership hard boundaries are disabled; using Unknown ownership.")
        return np.zeros((int(arcpy.Raster(age_raster).height), int(arcpy.Raster(age_raster).width)), dtype=np.int16)
    if not arcpy.Exists(PARCEL_OWNERSHIP_LAYER):
        raise FileNotFoundError(
            "Parcel ownership layer was not found. Update PARCEL_OWNERSHIP_LAYER:\n"
            + PARCEL_OWNERSHIP_LAYER
        )

    owner_field = _field_name(PARCEL_OWNERSHIP_LAYER, PARCEL_OWNER_FIELD)
    scratch_copy = os.path.join(output_gdb, "_ownership_classified")
    _delete_if_exists(scratch_copy, ownership_dissolved, ownership_raster)
    message("Classifying parcel ownership values...")
    arcpy.management.CopyFeatures(PARCEL_OWNERSHIP_LAYER, scratch_copy)
    _add_field_if_missing(scratch_copy, "OWN_CODE", "SHORT")
    _add_field_if_missing(scratch_copy, "OWN_TYPE", "TEXT", field_length=30)
    with arcpy.da.UpdateCursor(scratch_copy, [owner_field, "OWN_CODE", "OWN_TYPE"]) as cursor:
        for raw_owner, _, _ in cursor:
            code, owner_type = _classify_owner_value(raw_owner)
            cursor.updateRow((raw_owner, code, owner_type))

    message("Dissolving parcels by standardized ownership class...")
    arcpy.analysis.PairwiseDissolve(
        scratch_copy,
        ownership_dissolved,
        ["OWN_CODE", "OWN_TYPE"],
        multi_part="MULTI_PART",
    )
    message("Rasterizing ownership classes as hard stand boundaries...")
    arcpy.conversion.PolygonToRaster(
        ownership_dissolved,
        "OWN_CODE",
        ownership_raster,
        cell_assignment="MAXIMUM_AREA",
        priority_field="OWN_CODE",
        cellsize=float(arcpy.Raster(age_raster).meanCellWidth),
    )
    arcpy.management.DefineProjection(
        ownership_raster, arcpy.Describe(age_raster).spatialReference
    )
    ownership = arcpy.RasterToNumPyArray(
        ownership_raster, nodata_to_value=UNKNOWN_OWNERSHIP_CODE
    ).astype(np.int16)
    _delete_if_exists(scratch_copy)
    message(f"Dissolved ownership boundaries written:\n{ownership_dissolved}")
    message(f"Aligned ownership raster written:\n{ownership_raster}")
    return ownership


def finalize_database_groups(feature_class):
    """Route riparian units separately while retaining ownership attributes."""
    _add_field_if_missing(feature_class, "DB_GROUP", "TEXT", field_length=30)
    owner_type_field = _field_name(feature_class, "OWN_TYPE")
    with arcpy.da.UpdateCursor(
        feature_class, ["MGMT_CLASS", owner_type_field, "DB_GROUP"]
    ) as cursor:
        for mgmt_class, owner_type, _ in cursor:
            group = owner_type or UNKNOWN_OWNERSHIP_TYPE
            if int(mgmt_class or 0) == RIPARIAN_MANAGEMENT_CODE:
                group = "Riparian"
            cursor.updateRow((mgmt_class, owner_type, group))

def create_riparian_stands():
    """
    Create a dissolved NHD-based riparian footprint and split finished stands.

    Polygon water features are buffered from their shoreline. Flowlines are
    buffered from their centerline. Overlapping source buffers are dissolved
    into one binary footprint before the stands are split, preventing duplicate
    stand acreage where multiple NHD layers represent the same water feature.
    """
    required = [segment_polygons, nhd_flowline, nhd_area, nhd_waterbody]
    missing = [path for path in required if not arcpy.Exists(path)]
    if missing:
        raise FileNotFoundError(
            "The following riparian inputs were not found:\n" + "\n".join(missing)
        )

    scratch = {
        "aoi_core": os.path.join(output_gdb, "_rip_aoi_core"),
        "aoi": os.path.join(output_gdb, "_rip_aoi"),
        "flow_clip": os.path.join(output_gdb, "_rip_flow_clip"),
        "area_clip": os.path.join(output_gdb, "_rip_area_clip"),
        "water_clip": os.path.join(output_gdb, "_rip_water_clip"),
        "flow_class": os.path.join(output_gdb, "NHD_Flowline_Buffer_Rules"),
        "area_class": os.path.join(output_gdb, "NHD_Area_Buffer_Rules"),
        "water_class": os.path.join(output_gdb, "NHD_Waterbody_Buffer_Rules"),
        "flow_buffer": os.path.join(output_gdb, "_rip_flow_buffer"),
        "area_buffer": os.path.join(output_gdb, "_rip_area_buffer"),
        "water_buffer": os.path.join(output_gdb, "_rip_water_buffer"),
        "wetland_merged": os.path.join(output_gdb, "_nwi_merged"),
        "wetland_intersection": os.path.join(output_gdb, "_nwi_intersection"),
        "merged_buffer": os.path.join(output_gdb, "_rip_merged_buffer"),
        "riparian_piece": os.path.join(output_gdb, "_riparian_piece"),
        "nonriparian_piece": os.path.join(output_gdb, "_nonriparian_piece"),
        "merged_stands": os.path.join(output_gdb, "_riparian_merged_stands"),
    }

    _delete_if_exists(
        *scratch.values(),
        riparian_buffer,
        riparian_split_stands,
        wetland_footprint,
    )

    maximum_buffer = _maximum_configured_buffer()
    if maximum_buffer <= 0:
        raise ValueError("All configured riparian buffer distances are zero.")

    message("Creating a hydrography search area around the stand footprint...")
    arcpy.analysis.PairwiseDissolve(segment_polygons, scratch["aoi_core"])
    arcpy.analysis.PairwiseBuffer(
        scratch["aoi_core"],
        scratch["aoi"],
        f"{maximum_buffer:g} Feet",
        dissolve_option="ALL",
        method="PLANAR",
    )

    # Clip the national layers first so all later operations are local.
    message("Clipping NHD layers to the buffered stand footprint...")
    arcpy.analysis.PairwiseClip(nhd_flowline, scratch["aoi"], scratch["flow_clip"])
    arcpy.analysis.PairwiseClip(nhd_area, scratch["aoi"], scratch["area_clip"])
    arcpy.analysis.PairwiseClip(nhd_waterbody, scratch["aoi"], scratch["water_clip"])

    message("Classifying original stands from USFWS NWI wetland coverage...")
    _classify_original_stands_as_wetland(scratch["aoi_core"], scratch)

    counts = {
        "flow": _classify_hydrography(
            scratch["flow_clip"],
            scratch["flow_class"],
            "Flowline",
            FLOWLINE_BUFFER_RULES,
        ),
        "area": _classify_hydrography(
            scratch["area_clip"],
            scratch["area_class"],
            "Area",
            AREA_BUFFER_RULES,
        ),
        "water": _classify_hydrography(
            scratch["water_clip"],
            scratch["water_class"],
            "Waterbody",
            WATERBODY_BUFFER_RULES,
        ),
    }

    buffer_outputs = []
    for key, classified_key, buffer_key in (
        ("flow", "flow_class", "flow_buffer"),
        ("area", "area_class", "area_buffer"),
        ("water", "water_class", "water_buffer"),
    ):
        if counts[key] == 0:
            continue
        message(f"Buffering classified NHD {key} features...")
        arcpy.analysis.PairwiseBuffer(
            scratch[classified_key],
            scratch[buffer_key],
            "BUFFER_TXT",
            dissolve_option="NONE",
            method="PLANAR",
        )
        buffer_outputs.append(scratch[buffer_key])

    if not buffer_outputs:
        raise RuntimeError(
            "No supported NHD features occurred near the stand polygons. "
            "Review the input extent and FCode rules."
        )

    arcpy.management.Merge(buffer_outputs, scratch["merged_buffer"])
    arcpy.analysis.PairwiseDissolve(scratch["merged_buffer"], riparian_buffer)

    message("Splitting stands into riparian and non-riparian pieces...")
    arcpy.analysis.PairwiseIntersect(
        [segment_polygons, riparian_buffer],
        scratch["riparian_piece"],
        join_attributes="ALL",
    )
    arcpy.analysis.PairwiseErase(
        segment_polygons,
        riparian_buffer,
        scratch["nonriparian_piece"],
    )

    for feature_class, value, label in (
        (scratch["riparian_piece"], 1, "Riparian Buffer"),
        (scratch["nonriparian_piece"], 0, "Non-Riparian"),
    ):
        _add_field_if_missing(feature_class, "RIPARIAN", "SHORT")
        _add_field_if_missing(feature_class, "RIP_CLASS", "TEXT", field_length=40)
        arcpy.management.CalculateField(
            feature_class,
            "RIPARIAN",
            str(value),
            "PYTHON3",
        )
        arcpy.management.CalculateField(
            feature_class,
            "RIP_CLASS",
            repr(label),
            "PYTHON3",
        )

    arcpy.management.Merge(
        [scratch["riparian_piece"], scratch["nonriparian_piece"]],
        scratch["merged_stands"],
    )
    arcpy.management.MultipartToSinglepart(
        scratch["merged_stands"],
        riparian_split_stands,
    )
    _calculate_split_fields(riparian_split_stands)
    _assign_management_classes(riparian_split_stands)
    finalize_database_groups(riparian_split_stands)

    riparian_count = 0
    nonriparian_count = 0
    riparian_acres = 0.0
    nonriparian_acres = 0.0
    with arcpy.da.SearchCursor(
        riparian_split_stands,
        ["RIPARIAN", "Acres"],
    ) as cursor:
        for is_riparian, acres in cursor:
            if int(is_riparian) == 1:
                riparian_count += 1
                riparian_acres += float(acres or 0.0)
            else:
                nonriparian_count += 1
                nonriparian_acres += float(acres or 0.0)

    message(f"Riparian pieces: {riparian_count:,} ({riparian_acres:,.2f} acres)")
    message(
        f"Non-riparian pieces: {nonriparian_count:,} "
        f"({nonriparian_acres:,.2f} acres)"
    )
    message(f"Combined riparian buffer written:\n{riparian_buffer}")
    if CLASSIFY_WETLANDS and arcpy.Exists(wetland_footprint):
        message(f"Dissolved NWI wetland footprint written:\n{wetland_footprint}")
    message(f"Final management-classified stands written:\n{riparian_split_stands}")

    # Keep the three classified source layers for QA and future rule changes,
    # but remove temporary clips, buffers, and intermediate split layers.
    keep = {
        scratch["flow_class"],
        scratch["area_class"],
        scratch["water_class"],
    }
    for path in scratch.values():
        if path not in keep and arcpy.Exists(path):
            arcpy.management.Delete(path)


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
        [forest_type_raster, basal_area_raster, qmd_raster, tpa_raster],
    )

    arcpy.env.snapRaster = age_raster
    arcpy.env.cellSize = age_raster
    arcpy.env.extent = age_raster
    arcpy.env.outputCoordinateSystem = age_raster

    cell_area_map_units = geometry["cell_width"] * geometry["cell_height"]
    cell_acres = cell_area_map_units / 4046.8564224
    initial_seed_cells = max(1, int(round(INITIAL_SEED_ACRES / cell_acres)))
    minimum_stand_cells = max(1, int(math.ceil(MIN_STAND_ACRES / cell_acres)))
    maximum_stand_cells = max(
        minimum_stand_cells,
        int(math.floor(MAX_STAND_ACRES / cell_acres)),
    )

    message(f"Cell size: {geometry['cell_width']:.3f} × {geometry['cell_height']:.3f} map units")
    message(f"Cell area: {cell_acres:.6f} acres")
    message(f"Initial seed size: approximately {initial_seed_cells:,} cells")
    message(f"Minimum stand size: {minimum_stand_cells:,} cells")
    message(f"Maximum stand size: approximately {maximum_stand_cells:,} cells")

    message("=" * 60)
    message("READING SEGMENTATION RASTERS")
    message("=" * 60)

    age = read_raster(age_raster)
    forest_type = read_raster(forest_type_raster)
    basal_area = read_raster(basal_area_raster, scale=10.0)
    qmd = read_raster(qmd_raster, scale=10.0)
    tpa = read_raster(tpa_raster, scale=10.0)
    ownership = prepare_ownership_boundaries()
    if ownership.shape != age.shape:
        raise RuntimeError(
            f"Ownership raster shape {ownership.shape} does not match segmentation rasters {age.shape}."
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
        raise RuntimeError("No cells contained complete data across all five rasters.")
    message(f"Cells with complete segmentation data: {valid_cells:,}")

    message("=" * 60)
    message("STANDARDIZING CONTINUOUS VARIABLES")
    message("=" * 60)

    age_z, age_center, age_scale = robust_standardize(age, valid_mask)
    ba_z, ba_center, ba_scale = robust_standardize(basal_area, valid_mask)
    qmd_z, qmd_center, qmd_scale = robust_standardize(qmd, valid_mask)
    tpa_z, tpa_center, tpa_scale = robust_standardize(tpa, valid_mask)

    message(f"STDAGE center={age_center:.3f}, scale={age_scale:.3f}")
    message(f"BALIVE center={ba_center:.3f}, scale={ba_scale:.3f}")
    message(f"QMD center={qmd_center:.3f}, scale={qmd_scale:.3f}")
    message(f"TPA center={tpa_center:.3f}, scale={tpa_scale:.3f}")

    feature_arrays = {"STDAGE": age_z, "BALIVE": ba_z, "QMD": qmd_z, "TPA": tpa_z}
    raw_arrays = {"STDAGE": age, "BALIVE": basal_area, "QMD": qmd, "TPA": tpa}

    message("=" * 60)
    message("CREATING INITIAL SEEDS")
    message("=" * 60)
    start = time.perf_counter()
    labels = create_initial_segments(valid_mask, initial_seed_cells, ownership)
    initial_count = int(np.count_nonzero(np.bincount(labels.ravel())[1:]))
    message(f"Initial connected seed regions: {initial_count:,}")
    timed_message("Initial seed creation", start)

    message("=" * 60)
    message("ITERATIVE CELL REASSIGNMENT")
    message("=" * 60)
    labels = iterative_reassignment(labels, valid_mask, feature_arrays, forest_type, ownership)

    message("=" * 60)
    message("SEPARATING DISCONNECTED PIECES")
    message("=" * 60)
    start = time.perf_counter()
    labels = split_disconnected_segments(labels)
    pre_merge_count = int(np.count_nonzero(np.bincount(labels.ravel())[1:]))
    message(f"Connected segments before maximum-size splitting: {pre_merge_count:,}")
    timed_message("Disconnected-piece separation", start)

    message("=" * 60)
    message("SPLITTING OVERSIZED SEGMENTS")
    message("=" * 60)
    labels = split_oversized_segments(labels, maximum_stand_cells)
    labels = split_disconnected_segments(labels)
    count_after_split = int(np.count_nonzero(np.bincount(labels.ravel())[1:]))
    message(f"Connected segments after maximum-size splitting: {count_after_split:,}")

    message("=" * 60)
    message("MERGING UNDERSIZED SEGMENTS")
    message("=" * 60)
    labels = merge_small_segments(
        labels,
        feature_arrays,
        forest_type,
        ownership,
        minimum_stand_cells,
        maximum_cells=maximum_stand_cells,
    )
    labels = split_disconnected_segments(labels)

    counts = np.bincount(labels.ravel())
    remaining_small_before_delete = int(
        np.sum((counts[1:] > 0) & (counts[1:] < minimum_stand_cells))
    )
    message(
        f"Undersized segments remaining after merge: "
        f"{remaining_small_before_delete:,}"
    )

    diagnose_remaining_small_segments(
        labels,
        minimum_stand_cells,
        maximum_stand_cells,
    )

    if MERGE_SIMILAR_STANDS:
        message("=" * 60)
        message("MERGING SIMILAR ADJACENT STANDS")
        message("=" * 60)
        message(
            f"Criteria: age difference <= "
            f"{SIMILAR_MERGE_MAX_AGE_DIFFERENCE:.1f} years; "
            f"same dominant FORTYPCD="
            f"{SIMILAR_MERGE_REQUIRE_SAME_FORTYPCD}; "
            f"combined size <= {MAX_STAND_ACRES:.1f} acres"
        )
        labels = merge_similar_adjacent_stands(
            labels,
            age,
            forest_type,
            ownership,
            maximum_stand_cells,
            SIMILAR_MERGE_MAX_AGE_DIFFERENCE,
            require_same_forest_type=(
                SIMILAR_MERGE_REQUIRE_SAME_FORTYPCD
            ),
            minimum_shared_edges=SIMILAR_MERGE_MIN_SHARED_EDGES,
            maximum_passes=SIMILAR_MERGE_MAX_PASSES,
        )
        labels = split_disconnected_segments(labels)

        counts = np.bincount(labels.ravel())
        remaining_small_before_delete = int(
            np.sum(
                (counts[1:] > 0)
                & (counts[1:] < minimum_stand_cells)
            )
        )
        message(
            f"Undersized segments after similar-stand merging: "
            f"{remaining_small_before_delete:,}"
        )
        diagnose_remaining_small_segments(
            labels,
            minimum_stand_cells,
            maximum_stand_cells,
        )

    message("=" * 60)
    message("DELETING REMAINING UNDERSIZED SEGMENTS")
    message("=" * 60)
    labels, deleted_segments, deleted_cells = delete_small_segments(labels, minimum_stand_cells)
    labels = split_disconnected_segments(labels)
    labels = renumber_segments(labels)

    message(f"Deleted undersized segments: {deleted_segments:,}")
    message(f"Deleted forest cells: {deleted_cells:,} ({deleted_cells * cell_acres:,.2f} acres)")
    final_count = int(labels.max())
    message(f"Final segment count: {final_count:,}")

    message("=" * 60)
    message("CREATING SEGMENT SUMMARY")
    message("=" * 60)
    homogeneity_scales = {
        "STDAGE": float(np.nanstd(age[valid_mask])),
        "BALIVE": float(np.nanstd(basal_area[valid_mask])),
        "QMD": float(np.nanstd(qmd[valid_mask])),
        "TPA": float(np.nanstd(tpa[valid_mask])),
    }
    message("Landscape standard deviations used for homogeneity:")
    for variable, landscape_sd in homogeneity_scales.items():
        message(f"  {variable}: {landscape_sd:.3f}")

    start = time.perf_counter()
    summary = create_summary(labels, raw_arrays, forest_type, ownership, cell_acres, homogeneity_scales)
    timed_message("Segment summary", start)
    summary.to_csv(segment_summary_csv, index=False)
    message(f"Segment summary written:\n{segment_summary_csv}")

    message(f"Mean stand acreage: {summary['Acres'].mean():.2f}")
    message(f"Median stand acreage: {summary['Acres'].median():.2f}")
    message(f"Minimum stand acreage: {summary['Acres'].min():.2f}")
    message(f"Maximum stand acreage: {summary['Acres'].max():.2f}")

    mean_homogeneity = float(summary["HOMOGENEITY"].mean())
    area_weighted_homogeneity = float(np.average(summary["HOMOGENEITY"], weights=summary["Acres"]))
    message(f"Mean stand homogeneity score: {mean_homogeneity:.2f} / 100")
    message(f"Area-weighted landscape homogeneity score: {area_weighted_homogeneity:.2f} / 100")

    message("Mean within-stand standard deviations:")
    for variable in ["STDAGE", "BALIVE", "QMD", "TPA"]:
        sd_field = f"{variable}_SD"
        unweighted_mean_sd = float(summary[sd_field].mean())
        area_weighted_mean_sd = float(np.average(summary[sd_field], weights=summary["Acres"]))
        message(f"  {variable}: mean={unweighted_mean_sd:.3f}, area-weighted={area_weighted_mean_sd:.3f}")

    message("=" * 60)
    message("WRITING SEGMENT RASTER")
    message("=" * 60)
    write_segment_raster(labels, geometry)
    message(f"Segment raster written:\n{segment_raster}")

    message("=" * 60)
    message("CONVERTING SEGMENTS TO POLYGONS")
    message("=" * 60)
    write_segment_polygons(summary)
    message(f"Segment polygons written:\n{segment_polygons}")

    if CREATE_RIPARIAN_STANDS:
        message("=" * 60)
        message("CREATING RIPARIAN BUFFERS AND SPLIT STANDS")
        message("=" * 60)
        create_riparian_stands()

    message("=" * 60)
    message("SEGMENTATION COMPLETE")
    message("=" * 60)


if __name__ == "__main__":
    main()