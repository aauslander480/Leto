"""Vectorized ownership cleanup for TreeMap-covered NonForest patches."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage


@dataclass(frozen=True)
class OwnershipCleanupResult:
    """Cleaned ownership values, change mask, and scalar QA metrics."""

    ownership: np.ndarray
    changed_mask: np.ndarray
    qa: dict[str, object]


def reassign_nonforest_patches(
    ownership: np.ndarray,
    treemap_valid_mask: np.ndarray,
    *,
    nonforest_code: int,
    eligible_codes: tuple[int, ...] | list[int],
    cell_area_acres: float,
) -> OwnershipCleanupResult:
    """Reassign each targeted NonForest patch by longest shared cell edge.

    Target patches use 8-neighbor connectivity, matching cells that meet at a
    corner into one patch. Candidate boundary length deliberately uses only
    the four cardinal directions because diagonal contact shares no cell edge.
    Exact ties select the lowest eligible ownership code, making the result
    independent of raster traversal order. Patches with no eligible shared
    edge remain NonForest.

    Parameters
    ----------
    ownership:
        Two-dimensional aligned ownership-code array.
    treemap_valid_mask:
        Boolean array identifying valid TreeMap forest coverage.
    nonforest_code:
        Existing ownership code whose TreeMap-covered patches are targeted.
    eligible_codes:
        Existing ownership codes allowed as destinations. Unknown, NonForest,
        and Water should be excluded by configuration.
    cell_area_acres:
        Actual raster cell area converted to acres by the caller.
    """
    ownership = np.asarray(ownership)
    valid = np.asarray(treemap_valid_mask, dtype=bool)
    if ownership.ndim != 2 or valid.ndim != 2:
        raise ValueError("Ownership and TreeMap masks must be two-dimensional.")
    if ownership.shape != valid.shape:
        raise ValueError("Ownership and TreeMap masks must have identical shapes.")
    if not np.isfinite(cell_area_acres) or cell_area_acres <= 0:
        raise ValueError("cell_area_acres must be a positive finite number.")

    eligible = tuple(sorted({int(code) for code in eligible_codes}))
    if int(nonforest_code) in eligible:
        raise ValueError("NonForest cannot be an eligible replacement code.")

    target = valid & (ownership == int(nonforest_code))
    initial_cells = int(target.sum())
    structure = np.ones((3, 3), dtype=np.uint8)
    patch_labels, patch_count = ndimage.label(target, structure=structure)
    changed = np.zeros(ownership.shape, dtype=bool)
    cleaned = ownership.copy()

    if patch_count == 0 or not eligible:
        qa = _qa_summary(
            initial_cells, patch_count, changed, patch_labels, cleaned,
            eligible, cell_area_acres,
        )
        return OwnershipCleanupResult(cleaned, changed, qa)

    # Rows are patch IDs; columns are destination codes. Accumulate only
    # cardinal adjacency so each increment represents one shared cell edge.
    shared_edges = np.zeros((patch_count + 1, len(eligible)), dtype=np.int64)
    pairs = (
        (patch_labels[:-1, :], ownership[1:, :]),
        (patch_labels[1:, :], ownership[:-1, :]),
        (patch_labels[:, :-1], ownership[:, 1:]),
        (patch_labels[:, 1:], ownership[:, :-1]),
    )
    for column, code in enumerate(eligible):
        for patch_side, owner_side in pairs:
            touching = (patch_side > 0) & (owner_side == code)
            if np.any(touching):
                shared_edges[:, column] += np.bincount(
                    patch_side[touching], minlength=patch_count + 1
                )

    maxima = shared_edges.max(axis=1)
    # eligible is sorted, so argmax resolves exact ties to the lowest code.
    destinations = np.asarray(eligible, dtype=ownership.dtype)[
        np.argmax(shared_edges, axis=1)
    ]
    resolved_patch = maxima > 0
    resolved_patch[0] = False
    changeable = resolved_patch[patch_labels]
    cleaned[changeable] = destinations[patch_labels[changeable]]
    changed[changeable] = True

    qa = _qa_summary(
        initial_cells, patch_count, changed, patch_labels, cleaned,
        eligible, cell_area_acres,
    )
    return OwnershipCleanupResult(cleaned, changed, qa)


def _qa_summary(
    initial_cells: int,
    patch_count: int,
    changed: np.ndarray,
    patch_labels: np.ndarray,
    cleaned: np.ndarray,
    eligible: tuple[int, ...],
    cell_area_acres: float,
) -> dict[str, object]:
    changed_cells = int(changed.sum())
    reassigned_patch_ids = np.unique(patch_labels[changed])
    reassigned_patches = int(np.count_nonzero(reassigned_patch_ids))
    by_destination = {
        int(code): int(np.count_nonzero(changed & (cleaned == code)))
        for code in eligible
        if np.any(changed & (cleaned == code))
    }
    return {
        "initial_nonforest_cells": initial_cells,
        "initial_nonforest_acres": initial_cells * float(cell_area_acres),
        "patches_evaluated": int(patch_count),
        "patches_reassigned": reassigned_patches,
        "cells_reassigned": changed_cells,
        "acres_reassigned": changed_cells * float(cell_area_acres),
        "unresolved_patches": int(patch_count) - reassigned_patches,
        "cells_by_destination": by_destination,
        "acres_by_destination": {
            code: count * float(cell_area_acres)
            for code, count in by_destination.items()
        },
    }

