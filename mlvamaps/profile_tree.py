"""Deterministic neighbor joining for repeat-profile similarity exports."""
from __future__ import annotations

import numpy as np


def _newick_label(label: str) -> str:
    # Newick escapes a quote inside a quoted label by doubling it. Preserve
    # sample IDs exactly instead of silently rewriting meaningful characters.
    return "'" + str(label).replace("'", "''") + "'"


def neighbor_joining_tree_from_matrix(
    labels: list[str], distances: np.ndarray
) -> str:
    """Run neighbor joining on a dense matrix through compiled NumPy kernels."""
    if not labels:
        return ";\n"
    if len(labels) == 1:
        return f"({_newick_label(labels[0])}:0.00000000);\n"
    active = list(labels)
    clusters = {label: _newick_label(label) for label in labels}
    matrix = np.asarray(distances, dtype=np.float64).copy()
    if matrix.shape != (len(labels), len(labels)):
        raise ValueError("Neighbor-joining distance matrix shape does not match labels")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Neighbor-joining distances must be finite")
    np.maximum(matrix, 0.0, out=matrix)
    np.fill_diagonal(matrix, 0.0)

    node_number = 0
    while len(active) > 2:
        size = len(active)
        row_sums = matrix.sum(axis=1)
        q_matrix = (size - 2) * matrix - row_sums[:, None] - row_sums[None, :]
        q_matrix[np.tril_indices(size)] = np.inf
        left_index, right_index = np.unravel_index(
            int(np.argmin(q_matrix)), q_matrix.shape
        )
        left = active[left_index]
        right = active[right_index]
        pair_distance = float(matrix[left_index, right_index])
        left_length = 0.5 * pair_distance + (
            float(row_sums[left_index]) - float(row_sums[right_index])
        ) / (2 * (size - 2))
        right_length = pair_distance - left_length
        left_length = max(0.0, left_length)
        right_length = max(0.0, right_length)
        node_number += 1
        joined = f"__NJ_{node_number}"
        clusters[joined] = (
            f"({clusters[left]}:{left_length:.8f},"
            f"{clusters[right]}:{right_length:.8f})"
        )
        retained = [
            index for index in range(size) if index not in {left_index, right_index}
        ]
        joined_distances = np.maximum(
            0.0,
            (
                matrix[left_index, retained]
                + matrix[right_index, retained]
                - pair_distance
            )
            / 2,
        )
        next_matrix = np.zeros((size - 1, size - 1), dtype=float)
        next_matrix[:-1, :-1] = matrix[np.ix_(retained, retained)]
        next_matrix[-1, :-1] = joined_distances
        next_matrix[:-1, -1] = joined_distances
        matrix = next_matrix
        active = [active[index] for index in retained]
        active.append(joined)

    left, right = active
    final_distance = float(matrix[0, 1]) / 2
    return (
        f"({clusters[left]}:{max(0.0, final_distance):.8f},"
        f"{clusters[right]}:{max(0.0, final_distance):.8f});\n"
    )
