from __future__ import annotations

import math

from .models import Locus


def normalize_allele(value: float, precision: int = 6) -> int | float:
    """Return integer-valued alleles as ints and retain true partial alleles."""
    rounded = round(float(value), precision)
    return int(rounded) if rounded.is_integer() else rounded


def allele_grid(
    locus: Locus, step: float = 0.5, padding: float = 1.0
) -> list[int | float]:
    """Create the explicit allele states considered by probabilistic callers."""
    if step <= 0 or step > 1:
        raise ValueError("allele grid step must be greater than 0 and at most 1")
    if padding < 0:
        raise ValueError("allele grid padding cannot be negative")
    lower = max(0.0, float(locus.expected_min_repeats) - padding)
    upper = float(locus.expected_max_repeats) + padding
    count = int(math.floor((upper - lower) / step + 1e-9))
    values = [normalize_allele(lower + index * step) for index in range(count + 1)]
    if not values or float(values[-1]) < upper - 1e-9:
        values.append(normalize_allele(upper))
    return values


def gaussian_allele_probabilities(
    value: float,
    candidates: list[int | float],
    sigma: float,
) -> list[float]:
    """Evaluate and normalize a Gaussian measurement model on an allele grid."""
    if sigma <= 0:
        raise ValueError("allele measurement sigma must be positive")
    weights = [
        math.exp(-((float(value) - float(candidate)) ** 2) / (2 * sigma * sigma))
        for candidate in candidates
    ]
    total = sum(weights)
    if total <= 0:
        nearest = min(range(len(candidates)), key=lambda index: abs(float(candidates[index]) - value))
        return [1.0 if index == nearest else 0.0 for index in range(len(candidates))]
    return [weight / total for weight in weights]


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
    # MLVA_finder encoded these three values in names such as
    # ``vrrA_12bp_314bp_10U`` and used this exact calculation.  Keep that
    # convention when the metadata are available so assembly calls remain
    # directly comparable with historical profiles.
    if locus.expected_product_size_bp and locus.nominal_repeat_units:
        return abs(
            locus.nominal_repeat_units
            - ((locus.expected_product_size_bp - product_size_bp) / repeat_bp)
        )
    nonrepeat_bp = expected_nonrepeat_bp(locus)
    if nonrepeat_bp is None:
        nonrepeat_bp = len(locus.forward_primer) + len(locus.reverse_primer)
    repeat_region_bp = max(0, product_size_bp - nonrepeat_bp)
    return repeat_region_bp / repeat_bp


def legacy_round_repeat_count(value: float, tolerance: float = 0.25) -> int | float:
    """Round an assembly allele with the historical MLVA_finder convention.

    Values within ``tolerance`` of an integer become that integer; all other
    values become the intervening half allele.  The old default tolerance was
    0.25.
    """
    if not 0 <= tolerance <= 0.5:
        raise ValueError("repeat-count rounding tolerance must be between 0 and 0.5")
    lower = math.floor(value)
    upper = math.ceil(value)
    if value < lower + tolerance:
        return lower
    if value > upper - tolerance:
        return upper
    return lower + 0.5


def estimate_repeat_count_from_inner_length(locus: Locus, inner_size_bp: int) -> float | None:
    repeat_bp = repeat_unit_length(locus)
    if not repeat_bp:
        return None
    nonrepeat_bp = expected_nonrepeat_bp(locus)
    if nonrepeat_bp is None:
        return inner_size_bp / repeat_bp
    inner_nonrepeat_bp = max(0, nonrepeat_bp - len(locus.forward_primer) - len(locus.reverse_primer))
    return max(0, inner_size_bp - inner_nonrepeat_bp) / repeat_bp
