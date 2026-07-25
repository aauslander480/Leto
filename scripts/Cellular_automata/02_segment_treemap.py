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

output_gdb = os.path.join(output_folder, "TreeMap_Segmentation.gdb")
segment_raster = os.path.join(output_folder, "TreeMap_Segments_Initial.tif")
segment_polygons = os.path.join(output_gdb, "TreeMap_Segments_Initial")
segment_summary_csv = os.path.join(output_folder, "TreeMap_Segment_Summary.csv")


# ============================================================
# SEGMENTATION SETTINGS
# ============================================================

INITIAL_SEED_ACRES = 40.0
MIN_STAND_ACRES = 5.0
MAX_STAND_ACRES = 200.0

MAX_ITERATIONS = 40
CONVERGENCE_THRESHOLD = 0.001
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


def create_initial_segments(valid_mask, cells_per_seed):
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
            local_labels, component_count = ndimage.label(block_mask, structure=structure)
            target = labels[row_start:row_end, col_start:col_end]
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


def iterative_reassignment(labels, valid_mask, feature_arrays, forest_type):
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


def merge_small_segments(labels, feature_arrays, forest_type, minimum_cells, maximum_cells=None):
    """Merge many small segments per pass and remap the raster only once per pass."""
    pass_number = 0

    while True:
        pass_number += 1
        start = time.perf_counter()
        stats = segment_statistics_arrays(labels, feature_arrays, forest_type)
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

def create_summary(labels, raw_arrays, forest_type, cell_acres, homogeneity_scales):
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
    labels = create_initial_segments(valid_mask, initial_seed_cells)
    initial_count = int(np.count_nonzero(np.bincount(labels.ravel())[1:]))
    message(f"Initial connected seed regions: {initial_count:,}")
    timed_message("Initial seed creation", start)

    message("=" * 60)
    message("ITERATIVE CELL REASSIGNMENT")
    message("=" * 60)
    labels = iterative_reassignment(labels, valid_mask, feature_arrays, forest_type)

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
    summary = create_summary(labels, raw_arrays, forest_type, cell_acres, homogeneity_scales)
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

    message("=" * 60)
    message("SEGMENTATION COMPLETE")
    message("=" * 60)


if __name__ == "__main__":
    main()