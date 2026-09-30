from __future__ import annotations

import parasail


_DNA_MATRIX = parasail.matrix_create("ACGTN", 0, -1)


def alignment_metrics(query: str, target: str) -> dict[str, int | str]:
    """Return exact edit metrics from a SIMD Parasail global traceback."""
    if not query or not target:
        aligned_query = query or ("-" * len(target))
        aligned_target = target or ("-" * len(query))
        matches = 0
    else:
        result = parasail.nw_trace_scan_32(query, target, 1, 1, _DNA_MATRIX)
        traceback = result.traceback
        aligned_query = traceback.query
        aligned_target = traceback.ref
        matches = traceback.comp.count("|")

    insertions = aligned_target.count("-")
    deletions = aligned_query.count("-")
    substitutions = len(aligned_query) - insertions - deletions - matches
    return {
        "aligned_repeat_sequence": aligned_query,
        "aligned_representative_sequence": aligned_target,
        "insertions_vs_representative": insertions,
        "deletions_vs_representative": deletions,
        "substitutions_vs_representative": substitutions,
        "edit_distance_to_representative": insertions + deletions + substitutions,
    }


def alignment_metrics_pair(pair: tuple[str, str]) -> dict[str, int | str]:
    return alignment_metrics(*pair)