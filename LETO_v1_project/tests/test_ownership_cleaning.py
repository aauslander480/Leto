import numpy as np
import pytest

from leto.ownership_cleaning import reassign_nonforest_patches


NONFOREST = 1
WATER = 2
FAMILY = 3
CORPORATE = 4
FEDERAL = 6
ELIGIBLE = (FAMILY, CORPORATE, FEDERAL)


def clean(ownership, valid=None, cell_area=0.25):
    ownership = np.asarray(ownership, dtype=np.int16)
    if valid is None:
        valid = np.ones(ownership.shape, dtype=bool)
    return reassign_nonforest_patches(
        ownership,
        np.asarray(valid, dtype=bool),
        nonforest_code=NONFOREST,
        eligible_codes=ELIGIBLE,
        cell_area_acres=cell_area,
    )


def test_patch_surrounded_primarily_by_family_becomes_family():
    result = clean([
        [3, 3, 3, 3],
        [3, 1, 1, 4],
        [3, 1, 1, 4],
        [3, 3, 3, 4],
    ])
    assert np.all(result.ownership[1:3, 1:3] == FAMILY)
    assert result.qa["patches_reassigned"] == 1
    assert result.qa["cells_reassigned"] == 4


def test_longest_shared_edge_selects_corporate_for_entire_patch():
    result = clean([
        [4, 4, 4, 4, 4],
        [3, 1, 1, 1, 4],
        [3, 1, 1, 1, 4],
        [3, 3, 4, 4, 4],
    ])
    assert np.all(result.ownership[1:3, 1:4] == CORPORATE)
    assert len(np.unique(result.ownership[1:3, 1:4])) == 1


def test_diagonal_contact_does_not_count_as_shared_border():
    result = clean([
        [4, 2, 2],
        [2, 1, 3],
        [2, 3, 3],
    ])
    assert result.ownership[1, 1] == FAMILY


def test_diagonally_connected_nonforest_cells_are_one_patch():
    result = clean([
        [3, 1, 3],
        [3, 3, 1],
        [3, 3, 4],
    ])
    assert result.qa["patches_evaluated"] == 1
    assert result.ownership[0, 1] == FAMILY
    assert result.ownership[1, 2] == FAMILY


def test_nonforest_outside_treemap_coverage_is_unchanged():
    ownership = np.array([[3, 1, 3], [3, 1, 3]], dtype=np.int16)
    valid = np.array([[1, 1, 1], [1, 0, 1]], dtype=bool)
    result = clean(ownership, valid)
    assert result.ownership[0, 1] == FAMILY
    assert result.ownership[1, 1] == NONFOREST
    assert not result.changed_mask[1, 1]


def test_tie_uses_lowest_eligible_code_deterministically():
    result = clean([[2, 3, 2], [2, 1, 2], [2, 4, 2]])
    assert result.ownership[1, 1] == FAMILY


def test_patch_without_eligible_neighbor_remains_nonforest():
    result = clean([[2, 2, 2], [2, 1, 2], [2, 2, 2]])
    assert result.ownership[1, 1] == NONFOREST
    assert result.qa["unresolved_patches"] == 1
    assert not result.changed_mask.any()


def test_shape_dtype_untargeted_values_and_acreage_are_preserved():
    ownership = np.array([[0, 3, 4], [2, 1, 6]], dtype=np.int16)
    original = ownership.copy()
    result = clean(ownership, cell_area=0.37)
    assert result.ownership.shape == original.shape
    assert result.ownership.dtype == original.dtype
    assert np.array_equal(result.ownership[original != NONFOREST], original[original != NONFOREST])
    assert result.qa["initial_nonforest_acres"] == pytest.approx(0.37)
    assert result.qa["acres_reassigned"] == pytest.approx(0.37)


def test_mismatched_dimensions_fail_early():
    with pytest.raises(ValueError, match="identical shapes"):
        clean([[1, 3]], valid=[[1], [1]])

