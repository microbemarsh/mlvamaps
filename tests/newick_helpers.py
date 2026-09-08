"""Read generated Newick trees for independent export assertions."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class _Node:
    name: str | None
    children: list[tuple["_Node", float]]


def _parse_newick(text: str) -> _Node:
    """Parse the branch lengths and leaf names needed from a Newick tree."""
    text = text.strip()
    position = 0

    def skip_space() -> None:
        nonlocal position
        while position < len(text) and text[position].isspace():
            position += 1

    def label() -> str:
        nonlocal position
        skip_space()
        if position < len(text) and text[position] in {"'", '"'}:
            quote = text[position]
            position += 1
            parts = []
            while position < len(text):
                if text[position] == quote:
                    if position + 1 < len(text) and text[position + 1] == quote:
                        parts.append(quote)
                        position += 2
                        continue
                    position += 1
                    return "".join(parts)
                parts.append(text[position])
                position += 1
            raise ValueError("Unterminated quoted name in Newick tree")
        start = position
        while position < len(text) and text[position] not in ",():;":
            position += 1
        return text[start:position].strip()

    def branch_length() -> float:
        nonlocal position
        skip_space()
        if position >= len(text) or text[position] != ":":
            return 0.0
        position += 1
        start = position
        while position < len(text) and text[position] not in ",);":
            position += 1
        value = text[start:position].strip()
        try:
            length = float(value)
        except ValueError as exc:
            raise ValueError(f"Invalid branch length {value!r} in tree") from exc
        return length

    def subtree() -> tuple[_Node, float]:
        nonlocal position
        skip_space()
        if position < len(text) and text[position] == "(":
            position += 1
            children: list[tuple[_Node, float]] = []
            while True:
                children.append(subtree())
                skip_space()
                if position < len(text) and text[position] == ",":
                    position += 1
                    continue
                if position >= len(text) or text[position] != ")":
                    raise ValueError("Malformed Newick tree")
                position += 1
                break
            # A tree may put a support value or internal-node label here. It is
            # irrelevant for patristic distance but must still be consumed.
            internal_name = label()
            node = _Node(internal_name or None, children)
            length = branch_length()
            return node, length
        name = label()
        if not name:
            raise ValueError("Blank leaf name in Newick tree")
        node = _Node(name, [])
        length = branch_length()
        return node, length

    root, _root_length = subtree()
    skip_space()
    if position < len(text) and text[position] == ";":
        position += 1
    skip_space()
    if position != len(text):
        raise ValueError("Unexpected trailing content in Newick tree")
    return root


def _tip_patristic_distance_matrix(root: _Node) -> tuple[list[str], np.ndarray]:
    """Build a dense tip-distance matrix using NumPy's compiled matrix kernels."""
    edge_lengths: list[float] = []
    names: list[str] = []
    tip_paths: list[list[int]] = []

    def collect(node: _Node, path: list[int]) -> None:
        if not node.children and node.name is not None:
            names.append(str(node.name))
            tip_paths.append(path)
        for child, length in node.children:
            edge_index = len(edge_lengths)
            edge_lengths.append(float(length))
            collect(child, [*path, edge_index])

    collect(root, [])
    if len(names) != len(set(names)):
        raise ValueError("Reference tree contains duplicate tip names")
    incidence = np.zeros((len(names), len(edge_lengths)), dtype=np.float64)
    for tip_index, path in enumerate(tip_paths):
        incidence[tip_index, path] = 1.0
    weighted_paths = incidence * np.asarray(edge_lengths, dtype=np.float64)
    root_distances = weighted_paths.sum(axis=1)
    shared_distances = weighted_paths @ incidence.T
    matrix = (
        root_distances[:, None]
        + root_distances[None, :]
        - 2.0 * shared_distances
    )
    np.maximum(matrix, 0.0, out=matrix)
    return names, matrix


def _tip_patristic_distances(root: _Node) -> dict[tuple[str, str], float]:
    """Return the public pair-keyed view of the compiled dense calculation."""
    names, matrix = _tip_patristic_distance_matrix(root)
    return {
        tuple(sorted((left, right))): float(matrix[left_index, right_index])
        for left_index, left in enumerate(names)
        for right_index, right in enumerate(names[left_index + 1 :], left_index + 1)
    }
