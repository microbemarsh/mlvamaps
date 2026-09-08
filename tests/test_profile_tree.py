import numpy as np
import pytest

from mlvamaps.profile_tree import neighbor_joining_tree_from_matrix
from newick_helpers import _parse_newick, _tip_patristic_distances

def test_all_tip_distances_use_tree_branch_lengths():
    root = _parse_newick("((A:1,B:2):3,C:4);")
    assert _tip_patristic_distances(root) == {
        ("A", "B"): pytest.approx(3.0),
        ("A", "C"): pytest.approx(8.0),
        ("B", "C"): pytest.approx(9.0),
    }


def test_vectorized_neighbor_joining_is_deterministic():
    labels = ["A", "B", "C", "D"]
    distances = {
        ("A", "B"): 5.0,
        ("A", "C"): 9.0,
        ("A", "D"): 9.0,
        ("B", "C"): 10.0,
        ("B", "D"): 10.0,
        ("C", "D"): 8.0,
    }
    matrix = np.zeros((len(labels), len(labels)))
    for (left, right), value in distances.items():
        i, j = labels.index(left), labels.index(right)
        matrix[i, j] = matrix[j, i] = value
    first = neighbor_joining_tree_from_matrix(labels, matrix)
    assert neighbor_joining_tree_from_matrix(labels, matrix) == first
    assert set(_tip_patristic_distances(_parse_newick(first))) == set(distances)
