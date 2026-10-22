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

from . import gis_backend as gis
import numpy as np
import pandas as pd
from scipy import ndimage

from .small_stands import plan_proximity_merges


# ============================================================
# USER INPUTS
# ============================================================

input_folder = (
    r"D:\Misc\Maine_5_co\ME_Output_stands"
)

age_raster = os.path.join(input_folder, "TreeMap_STDAGE.tif")
forest_type_raster = os.path.join(input_folder, "TreeMap_FORTYPCD.tif")
basal_area_raster = os.path.join(input_folder, "TreeMap_BALIVE.tif")
qmd_raster = os.path.join(input_folder, "TreeMap_QMD.tif")
tpa_raster = os.path.join(input_folder, "TreeMap_TPA.tif")

output_folder = os.path.join(
    r"D:\Misc\Maine_5_co\ME_Output_stands"
)

output_gdb = os.path.join(output_folder, "TreeMap_Segmentation.workspace")
segment_raster = os.path.join(output_folder, "TreeMap_Segments_Initial.tif")
segment_polygons = os.path.join(output_gdb, "TreeMap_Segments_Initial")
segment_summary_csv = os.path.join(output_folder, "TreeMap_Segment_Summary.csv")

# NHDPlus High Resolution national hydrography inputs.
nhd_flowline = (
    r"D:\Misc\Stream_buffer_data\NHDPlus_H_National_Release_2_GDB"
    r"\NHDPlus_H_National_Release_2_GDB.gdb\NHDPlus_H_National"
    r"\NetworkNHDFlowline"
)
nhd_area = (
    r"D:\Misc\Stream_buffer_data\NHDPlus_H_National_Release_2_GDB"
    r"\NHDPlus_H_National_Release_2_GDB.gdb\NHDPlus_H_National"
    r"\NHDArea"
)
nhd_waterbody = (
    r"D:\Misc\Stream_buffer_data\NHDPlus_H_National_Release_2_GDB"
    r"\NHDPlus_H_National_Release_2_GDB.gdb\NHDPlus_H_National"
    r"\NHDWaterbody"
)

# USFWS National Wetlands Inventory polygon inputs.
# Add one or more state feature classes or shapefiles to this list. The script
# clips, merges, and dissolves them automatically, so state boundaries do not
# affect wetland-percentage calculations.
NWI_WETLAND_INPUTS = []

# Riparian-processing outputs.
riparian_buffer = os.path.join(output_gdb, "Riparian_Buffer")
management_units_final = os.path.join(
    output_gdb,
    "ManagementUnits_Final",
)
# Backward-compatible internal alias used by the riparian functions.
riparian_split_stands = management_units_final
ownership_aligned_raster = os.path.join(
    output_folder,
    "Ownership_Classes_Aligned.tif",
)
wetland_footprint = os.path.join(output_gdb, "NWI_Wetlands_Dissolved")


# ============================================================
# SEGMENTATION SETTINGS
# ============================================================

# Populated from the selected regional configuration by pipeline.py.
INITIAL_SEED_ACRES = None
MIN_STAND_ACRES = None
MAX_STAND_ACRES = None
MAX_ITERATIONS = None
CONVERGENCE_THRESHOLD = None
MINIMUM_SCORE_IMPROVEMENT = None
USE_EIGHT_NEIGHBORS = None
WEIGHTS = {}
SHARED_EDGE_BONUS = None
Z_CLIP = None
READ_NODATA = None
REMAINING_SMALL_STAND_ACTION = None
PROXIMITY_MAXIMUM_DISTANCE_METERS = None
PROXIMITY_REQUIRE_SAME_OWNERSHIP = None
PROXIMITY_REQUIRE_SAME_FOREST_TYPE = None
PROXIMITY_MAXIMUM_COMBINED_ACRES = None
PROXIMITY_RETAIN_WHEN_NO_CANDIDATE = None
PROXIMITY_QA_CSV = None
PROXIMITY_MERGE_QA = {}
REPORT_TIMING = None
MERGE_SIMILAR_STANDS = None
SIMILAR_MERGE_MAX_AGE_DIFFERENCE = None
SIMILAR_MERGE_REQUIRE_SAME_FORTYPCD = None
SIMILAR_MERGE_MIN_SHARED_EDGES = None
SIMILAR_MERGE_MAX_PASSES = None


# ============================================================
# RIPARIAN BUFFER SETTINGS
# ============================================================

# Set to False to run only the original stand segmentation.
CREATE_RIPARIAN_STANDS = None

# These are editable first-pass assumptions. They are intentionally stored in
# one place so other states or later Florida BMP refinements can use different
# distances without changing the geometry-processing functions.
# Populated from the selected regional configuration by pipeline.py.
FLOWLINE_BUFFER_RULES = {}
WATERBODY_BUFFER_RULES = {}
AREA_BUFFER_RULES = {}

# No minimum-area deletion is performed after the riparian split. Narrow
# riparian pieces are intentional management units rather than segmentation
# errors.


# ============================================================
# WETLAND AND MANAGEMENT-CLASS SETTINGS
# ============================================================

# Classify each original TreeMap stand from all polygons in the listed USFWS
# National Wetlands Inventory inputs. No attribute filtering is performed.
# The final riparian-split pieces inherit the original-stand classification.
CLASSIFY_WETLANDS = None
WETLAND_MAJORITY_THRESHOLD = None

# Final management codes. Riparian always has first priority.
UPLAND_MANAGEMENT_CODE = None
RIPARIAN_MANAGEMENT_CODE = None
WETLAND_MANAGEMENT_CODE = None




# ============================================================
# OWNERSHIP RASTER SETTINGS
# ============================================================

USE_OWNERSHIP_HARD_BOUNDARIES = None

# National Forest Service ownership raster from RDS-2025-0045.
# Update this path to the downloaded raster on your computer.
OWNERSHIP_SOURCE_RASTER = (
    r"D:\forest_condition_2026\Stand_delineation\RDS-2025-0045 (1)\Data\US_forest_ownership.tif"
)

OWNERSHIP_TYPES = {}
UNKNOWN_OWNERSHIP_CODE = None
UNKNOWN_OWNERSHIP_TYPE = None

# Non-forest and water are not ownership classes. Converting them to
# Unknown prevents them from creating extra hard boundaries.
COLLAPSE_NONFOREST_AND_WATER = None
NONFOREST_AND_WATER_CODES = ()

ownership_projected_temp = os.path.join(
    output_gdb,
    "_ownership_projected",
)


# ============================================================
# GENERAL SETTINGS
# ============================================================

gis.env.overwriteOutput = True
# Workspace folders and geodatabases are initialized by run_leto.py.


# ============================================================
# HELPERS
# ============================================================

def message(text):
    print(text, flush=True)
    try:
        gis.AddMessage(text)
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

    if USE_OWNERSHIP_HARD_BOUNDARIES:
        paths.append(OWNERSHIP_SOURCE_RASTER)
    missing = [path for path in paths if not gis.Exists(path)]
    if missing:
        raise FileNotFoundError(
            "The following input rasters were not found:\n" + "\n".join(missing)
        )


def raster_geometry(raster_path):
    desc = gis.Describe(raster_path)
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
    array = gis.RasterToNumPyArray(path, nodata_to_value=READ_NODATA).astype(np.float32)
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
    data["OWN_TYPE"] = [
        OWNERSHIP_TYPES.get(
            int(code),
            UNKNOWN_OWNERSHIP_TYPE,
        )
        for code in data["OWN_CODE"]
    ]
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
    lower_left = gis.Point(geometry["xmin"], geometry["ymin"])
    raster = gis.NumPyArrayToRaster(
        export_array,
        lower_left,
        geometry["cell_width"],
        geometry["cell_height"],
        value_to_nodata=output_nodata,
    )
    raster.save(segment_raster)
    gis.management.DefineProjection(segment_raster, geometry["spatial_reference"])
    gis.management.CalculateStatistics(segment_raster)
    try:
        gis.management.BuildRasterAttributeTable(segment_raster, "Overwrite")
    except gis.ExecuteError:
        message("Raster attribute table was not built, but the segment raster was created.")


def write_segment_polygons(summary):
    if gis.Exists(segment_polygons):
        gis.management.Delete(segment_polygons)

    gis.conversion.RasterToPolygon(
        in_raster=segment_raster,
        out_polygon_features=segment_polygons,
        simplify="NO_SIMPLIFY",
        raster_field="Value",
        create_multipart_features="SINGLE_OUTER_PART",
    )

    existing_fields = {field.name.upper() for field in gis.ListFields(segment_polygons)}
    if "SEG_ID" not in existing_fields:
        gis.management.AddField(segment_polygons, "SEG_ID", "LONG")

    value_field = None
    for field in gis.ListFields(segment_polygons):
        if field.name.upper() in {"GRIDCODE", "GRID_CODE"}:
            value_field = field.name
            break
    if value_field is None:
        raise RuntimeError("Could not identify the raster-value field in the polygon output.")

    gis.management.CalculateField(segment_polygons, "SEG_ID", f"!{value_field}!", "PYTHON3")

    summary_table = os.path.join(output_gdb, "TreeMap_Segment_Summary")
    if gis.Exists(summary_table):
        gis.management.Delete(summary_table)

    gis.conversion.TableToTable(segment_summary_csv, output_gdb, os.path.basename(summary_table))
    join_fields = [field for field in summary.columns if field != "SEG_ID"]
    gis.management.JoinField(segment_polygons, "SEG_ID", summary_table, "SEG_ID", join_fields)


# ============================================================
# RIPARIAN BUFFERING AND STAND SPLITTING
# ============================================================

def _field_name(feature_class, requested_name):
    """Return the actual field name using a case-insensitive match."""
    requested_upper = requested_name.upper()
    for field in gis.ListFields(feature_class):
        if field.name.upper() == requested_upper:
            return field.name
    raise RuntimeError(
        f"Field '{requested_name}' was not found in:\n{feature_class}"
    )


def _delete_if_exists(*paths):
    for path in paths:
        if gis.Exists(path):
            gis.management.Delete(path)


def _add_field_if_missing(feature_class, name, field_type, **kwargs):
    existing = {field.name.upper() for field in gis.ListFields(feature_class)}
    if name.upper() not in existing:
        gis.management.AddField(feature_class, name, field_type, **kwargs)


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
    gis.management.CopyFeatures(clipped_features, output_features)

    fcode_field = _field_name(output_features, "FCode")
    _add_field_if_missing(output_features, "HYDRO_SRC", "TEXT", field_length=20)
    _add_field_if_missing(output_features, "HYDRO_CLASS", "TEXT", field_length=80)
    _add_field_if_missing(output_features, "BUFFER_FT", "DOUBLE")
    _add_field_if_missing(output_features, "BUFFER_TXT", "TEXT", field_length=30)
    _add_field_if_missing(output_features, "RULE_SOURCE", "TEXT", field_length=80)

    kept = 0
    removed = 0
    with gis.da.UpdateCursor(
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
    with gis.da.UpdateCursor(
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
        path for path in NWI_WETLAND_INPUTS if not gis.Exists(path)
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
        gis.analysis.PairwiseClip(wetland_input, stand_aoi, clip_output)

        clipped_count = int(gis.management.GetCount(clip_output)[0])
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
        gis.management.CopyFeatures(clip_outputs[0], wetland_merged)
    else:
        gis.management.Merge(clip_outputs, wetland_merged)

    # Dissolving removes overlap within and between state NWI inputs.
    gis.analysis.PairwiseDissolve(wetland_merged, wetland_footprint)
    gis.analysis.PairwiseIntersect(
        [segment_polygons, wetland_footprint],
        wetland_intersection,
        join_attributes="ALL",
    )

    seg_id_field = _field_name(wetland_intersection, "SEG_ID")
    wetland_acres_by_segment = defaultdict(float)
    with gis.da.SearchCursor(
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

    with gis.da.UpdateCursor(
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

    with gis.da.UpdateCursor(feature_class, fields + ["Acres"]) as cursor:
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

    with gis.da.UpdateCursor(
        feature_class,
        ["OID@", seg_id_field, "FINAL_ID", "MU_ID", "ParentID", "Acres", "SHAPE@"],
    ) as cursor:
        for row in cursor:
            row[2] = int(row[0])
            if row[3] is None:
                row[3] = int(row[0])
            row[4] = f"SEG_{int(row[1])}"
            row[5] = float(row[6].getArea("GEODESIC", "ACRES"))
            cursor.updateRow(row)



def _same_spatial_reference(first, second):
    """Return True when two ArcGIS spatial references are equivalent."""
    if first is None or second is None:
        return False

    first_code = int(getattr(first, "factoryCode", 0) or 0)
    second_code = int(getattr(second, "factoryCode", 0) or 0)

    if first_code > 0 and second_code > 0:
        return first_code == second_code

    return first.exportToString() == second.exportToString()


def prepare_ownership_raster(reference_raster):
    """
    Project, clip, and align the national ownership raster to the model grid.

    Nearest-neighbor resampling is used because ownership is categorical.
    """
    reference = gis.Raster(reference_raster)
    expected_shape = (int(reference.height), int(reference.width))

    if not USE_OWNERSHIP_HARD_BOUNDARIES:
        message("Ownership hard boundaries are disabled.")
        return np.full(
            expected_shape,
            UNKNOWN_OWNERSHIP_CODE,
            dtype=np.int16,
        )

    if not gis.Exists(OWNERSHIP_SOURCE_RASTER):
        raise FileNotFoundError(
            "Ownership source raster was not found. Update "
            "OWNERSHIP_SOURCE_RASTER:\n"
            f"{OWNERSHIP_SOURCE_RASTER}"
        )

    _delete_if_exists(
        ownership_projected_temp,
        ownership_aligned_raster,
    )

    reference_desc = gis.Describe(reference_raster)
    source_desc = gis.Describe(OWNERSHIP_SOURCE_RASTER)
    processing_source = OWNERSHIP_SOURCE_RASTER

    with gis.EnvManager(
        snapRaster=reference_raster,
        cellSize=reference_raster,
        extent=reference_raster,
        outputCoordinateSystem=reference_desc.spatialReference,
        resamplingMethod="NEAREST",
    ):
        if not _same_spatial_reference(
            source_desc.spatialReference,
            reference_desc.spatialReference,
        ):
            message(
                "Projecting ownership raster with nearest-neighbor "
                "resampling..."
            )
            gis.management.ProjectRaster(
                OWNERSHIP_SOURCE_RASTER,
                ownership_projected_temp,
                reference_desc.spatialReference,
                "NEAREST",
                f"{float(reference.meanCellWidth)} "
                f"{float(reference.meanCellHeight)}",
            )
            processing_source = ownership_projected_temp

        message(
            "Clipping and aligning ownership raster to the "
            "segmentation grid..."
        )
        aligned = gis.sa.ExtractByMask(
            processing_source,
            reference_raster,
        )
        aligned.save(ownership_aligned_raster)

    validate_alignment(
        reference_raster,
        [ownership_aligned_raster],
    )

    ownership = gis.RasterToNumPyArray(
        ownership_aligned_raster,
        nodata_to_value=UNKNOWN_OWNERSHIP_CODE,
    ).astype(np.int16)

    if ownership.shape != expected_shape:
        raise RuntimeError(
            "Aligned ownership raster shape does not match the "
            "segmentation raster.\n"
            f"Ownership: {ownership.shape}\n"
            f"Expected: {expected_shape}"
        )

    observed_codes = {int(value) for value in np.unique(ownership)}
    unexpected_codes = sorted(
        observed_codes - set(OWNERSHIP_TYPES)
    )
    if unexpected_codes:
        raise RuntimeError(
            "Ownership raster contains undocumented values: "
            f"{unexpected_codes}"
        )

    if COLLAPSE_NONFOREST_AND_WATER:
        collapse_mask = np.isin(
            ownership,
            NONFOREST_AND_WATER_CODES,
        )
        collapsed_cells = int(collapse_mask.sum())
        ownership[collapse_mask] = UNKNOWN_OWNERSHIP_CODE
        message(
            "NonForest and Water cells converted to Unknown: "
            f"{collapsed_cells:,}"
        )

    _delete_if_exists(ownership_projected_temp)

    message(
        "Aligned ownership raster written:\n"
        f"{ownership_aligned_raster}"
    )
    return ownership


def report_ownership_cells(ownership, valid_mask):
    """Print ownership-cell counts within the valid TreeMap footprint."""
    message("Ownership cells within the valid TreeMap footprint:")
    values, counts = np.unique(
        ownership[valid_mask],
        return_counts=True,
    )
    count_lookup = {
        int(value): int(count)
        for value, count in zip(values, counts)
    }

    for code in sorted(OWNERSHIP_TYPES):
        message(
            f"  {code} - {OWNERSHIP_TYPES[code]}: "
            f"{count_lookup.get(code, 0):,} cells"
        )


def finalize_database_groups(feature_class):
    """Route riparian units separately while retaining ownership attributes."""
    _add_field_if_missing(feature_class, "DB_GROUP", "TEXT", field_length=30)
    owner_type_field = _field_name(feature_class, "OWN_TYPE")
    with gis.da.UpdateCursor(
        feature_class, ["MGMT_CLASS", owner_type_field, "DB_GROUP"]
    ) as cursor:
        for mgmt_class, owner_type, _ in cursor:
            group = owner_type or UNKNOWN_OWNERSHIP_TYPE
            if int(mgmt_class or 0) == RIPARIAN_MANAGEMENT_CODE:
                group = "Riparian"
            cursor.updateRow((mgmt_class, owner_type, group))

def apply_proximity_merge():
    """Assign shared MU_ID values to eligible nearby segment polygons."""
    global PROXIMITY_MERGE_QA

    required_fields = ["OID@", "SEG_ID", "Acres", "OWN_CODE", "FORTYPCD_DOM", "STDAGE_MEAN"]
    units = []
    oid_to_segment = {}
    with gis.da.SearchCursor(segment_polygons, required_fields) as cursor:
        for oid, seg_id, acres, owner, forest_type, age in cursor:
            seg_id = int(seg_id)
            oid_to_segment[int(oid)] = seg_id
            units.append({
                "SEG_ID": seg_id,
                "Acres": acres,
                "OWN_CODE": owner,
                "FORTYPCD_DOM": forest_type,
                "STDAGE_MEAN": age,
            })

    near_table = os.path.join(output_gdb, "_proximity_merge_near")
    _delete_if_exists(near_table)
    gis.analysis.GenerateNearTable(
        in_features=segment_polygons,
        near_features=segment_polygons,
        out_table=near_table,
        search_radius=f"{float(PROXIMITY_MAXIMUM_DISTANCE_METERS):g} Meters",
        location="NO_LOCATION",
        angle="NO_ANGLE",
        closest="ALL",
        method="GEODESIC",
    )
    candidate_rows = []
    with gis.da.SearchCursor(
        near_table, ["IN_FID", "NEAR_FID", "NEAR_DIST"]
    ) as cursor:
        for source_oid, target_oid, distance in cursor:
            source_id = oid_to_segment.get(int(source_oid))
            target_id = oid_to_segment.get(int(target_oid))
            if source_id is None or target_id is None:
                continue
            candidate_rows.append({
                "SOURCE_SEG_ID": source_id,
                "TARGET_SEG_ID": target_id,
                "DISTANCE_METERS": float(distance),
            })

    mapping, merges, qa = plan_proximity_merges(
        pd.DataFrame(units),
        pd.DataFrame(
            candidate_rows,
            columns=["SOURCE_SEG_ID", "TARGET_SEG_ID", "DISTANCE_METERS"],
        ),
        minimum_acres=float(MIN_STAND_ACRES),
        maximum_distance_meters=float(PROXIMITY_MAXIMUM_DISTANCE_METERS),
        maximum_combined_acres=float(PROXIMITY_MAXIMUM_COMBINED_ACRES),
        require_same_ownership=bool(PROXIMITY_REQUIRE_SAME_OWNERSHIP),
        require_same_forest_type=bool(PROXIMITY_REQUIRE_SAME_FOREST_TYPE),
        retain_when_no_candidate=bool(PROXIMITY_RETAIN_WHEN_NO_CANDIDATE),
    )
    lookup = mapping.set_index("SEG_ID").to_dict(orient="index")
    _add_field_if_missing(segment_polygons, "MU_ID", "LONG")
    _add_field_if_missing(segment_polygons, "MU_PART", "LONG")
    _add_field_if_missing(
        segment_polygons, "MU_SEG_IDS", "TEXT", field_length=4000
    )
    with gis.da.UpdateCursor(
        segment_polygons, ["SEG_ID", "MU_ID", "MU_PART", "MU_SEG_IDS"]
    ) as cursor:
        for row in cursor:
            assignment = lookup[int(row[0])]
            row[1] = int(assignment["MU_ID"])
            row[2] = int(assignment["MU_PART"])
            row[3] = assignment["MU_SEG_IDS"]
            cursor.updateRow(row)

    if PROXIMITY_QA_CSV:
        merges.to_csv(PROXIMITY_QA_CSV, index=False)
    PROXIMITY_MERGE_QA = qa
    message(f"Proximity merges: {qa['proximity_merges']:,}")
    message(
        "Proximity merge distance (average / maximum): "
        f"{qa['average_merge_distance_meters']:.2f} / "
        f"{qa['maximum_merge_distance_meters']:.2f} meters"
    )
    message(
        f"Multipart management units: {qa['multipart_management_units']:,}; "
        "average segments per multipart unit: "
        f"{qa['average_segments_per_multipart_unit']:.2f}"
    )
    message(f"Small stands remaining: {qa['small_stands_remaining']:,}")
    message(
        "Proximity acreage before / after: "
        f"{qa['total_acres_before']:.6f} / {qa['total_acres_after']:.6f}"
    )
    _delete_if_exists(near_table)
    return qa


def create_riparian_stands():
    """
    Create a dissolved NHD-based riparian footprint and split finished stands.

    Polygon water features are buffered from their shoreline. Flowlines are
    buffered from their centerline. Overlapping source buffers are dissolved
    into one binary footprint before the stands are split, preventing duplicate
    stand acreage where multiple NHD layers represent the same water feature.
    """
    required = [segment_polygons, nhd_flowline, nhd_area, nhd_waterbody]
    missing = [path for path in required if not gis.Exists(path)]
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
    gis.analysis.PairwiseDissolve(segment_polygons, scratch["aoi_core"])
    gis.analysis.PairwiseBuffer(
        scratch["aoi_core"],
        scratch["aoi"],
        f"{maximum_buffer:g} Feet",
        dissolve_option="ALL",
        method="PLANAR",
    )

    # Clip the national layers first so all later operations are local.
    message("Clipping NHD layers to the buffered stand footprint...")
    gis.analysis.PairwiseClip(nhd_flowline, scratch["aoi"], scratch["flow_clip"])
    gis.analysis.PairwiseClip(nhd_area, scratch["aoi"], scratch["area_clip"])
    gis.analysis.PairwiseClip(nhd_waterbody, scratch["aoi"], scratch["water_clip"])

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
        gis.analysis.PairwiseBuffer(
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

    gis.management.Merge(buffer_outputs, scratch["merged_buffer"])
    gis.analysis.PairwiseDissolve(scratch["merged_buffer"], riparian_buffer)

    message("Splitting stands into riparian and non-riparian pieces...")
    gis.analysis.PairwiseIntersect(
        [segment_polygons, riparian_buffer],
        scratch["riparian_piece"],
        join_attributes="ALL",
    )
    gis.analysis.PairwiseErase(
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
        gis.management.CalculateField(
            feature_class,
            "RIPARIAN",
            str(value),
            "PYTHON3",
        )
        gis.management.CalculateField(
            feature_class,
            "RIP_CLASS",
            repr(label),
            "PYTHON3",
        )

    gis.management.Merge(
        [scratch["riparian_piece"], scratch["nonriparian_piece"]],
        scratch["merged_stands"],
    )
    gis.management.MultipartToSinglepart(
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
    with gis.da.SearchCursor(
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
    if CLASSIFY_WETLANDS and gis.Exists(wetland_footprint):
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
        if path not in keep and gis.Exists(path):
            gis.management.Delete(path)


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

    gis.env.snapRaster = age_raster
    gis.env.cellSize = age_raster
    gis.env.extent = age_raster
    gis.env.outputCoordinateSystem = age_raster

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
    ownership = prepare_ownership_raster(age_raster)
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
    report_ownership_cells(ownership, valid_mask)

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
    message("HANDLING REMAINING UNDERSIZED SEGMENTS")
    message("=" * 60)
    deleted_segments = 0
    deleted_cells = 0
    action = str(REMAINING_SMALL_STAND_ACTION).lower().strip()
    if action == "force_merge":
        labels = merge_small_segments(
            labels, feature_arrays, forest_type, ownership,
            minimum_stand_cells, maximum_cells=None,
        )
        labels = split_disconnected_segments(labels)
        labels = renumber_segments(labels)
        message("Forced additional ownership-compatible merges without a maximum-size limit.")
    elif action == "retain":
        labels = split_disconnected_segments(labels)
        labels = renumber_segments(labels)
        message("Remaining undersized segments were retained.")
    elif action == "delete":
        labels, deleted_segments, deleted_cells = delete_small_segments(labels, minimum_stand_cells)
        labels = split_disconnected_segments(labels)
        labels = renumber_segments(labels)
        message(f"Deleted undersized segments: {deleted_segments:,}")
        message(f"Deleted forest cells: {deleted_cells:,} ({deleted_cells * cell_acres:,.2f} acres)")
    elif action == "proximity_merge":
        message(
            "Residual small segments retained for non-geometric proximity "
            "grouping after polygon creation."
        )
    else:
        raise ValueError(
            "REMAINING_SMALL_STAND_ACTION must be 'delete', 'retain', "
            "'force_merge', or 'proximity_merge'."
        )
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

    if str(REMAINING_SMALL_STAND_ACTION).lower().strip() == "proximity_merge":
        message("=" * 60)
        message("GROUPING NEARBY UNDERSIZED MANAGEMENT UNITS")
        message("=" * 60)
        apply_proximity_merge()

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
