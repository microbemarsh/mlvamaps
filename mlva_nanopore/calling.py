from __future__ import annotations

from .models import Locus


def repeat_unit_length(locus: Locus) -> int:
    if locus.repeat_unit_length_bp:
        return locus.repeat_unit_length_bp
    if locus.repeat_motif:
        return len(locus.repeat_motif)
    return 0


def expected_nonrepeat_bp(locus: Locus) -> int | None:
    repeat_bp = repeat_unit_length(locus)
    if not repeat_bp or not locus.expected_product_size_bp or not locus.nominal_repeat_units:
        return None
    nonrepeat = locus.expected_product_size_bp - (locus.nominal_repeat_units * repeat_bp)
    return max(nonrepeat, len(locus.forward_primer) + len(locus.reverse_primer))


def estimate_repeat_count_from_product_length(locus: Locus, product_size_bp: int) -> float | None:
    repeat_bp = repeat_unit_length(locus)
    if not repeat_bp:
        return None
    nonrepeat_bp = expected_nonrepeat_bp(locus)
    if nonrepeat_bp is None:
        nonrepeat_bp = len(locus.forward_primer) + len(locus.reverse_primer)
    repeat_region_bp = max(0, product_size_bp - nonrepeat_bp)
    return repeat_region_bp / repeat_bp


def estimate_repeat_count_from_inner_length(locus: Locus, inner_size_bp: int) -> float | None:
    repeat_bp = repeat_unit_length(locus)
    if not repeat_bp:
        return None
    nonrepeat_bp = expected_nonrepeat_bp(locus)
    if nonrepeat_bp is None:
        return inner_size_bp / repeat_bp
    inner_nonrepeat_bp = max(0, nonrepeat_bp - len(locus.forward_primer) - len(locus.reverse_primer))
    return max(0, inner_size_bp - inner_nonrepeat_bp) / repeat_bp
