"""Canonical MLVA product-length and repeat-count calibration.

Every input technology uses these functions when it directly observes a
primer-bounded product.  The compatibility imports in :mod:`mlvamaps.calling`
remain available, but new callers should import this module explicitly.
"""

from __future__ import annotations

from .calling import (
    assembly_equivalent_product_allele,
    estimate_repeat_count_from_inner_length,
    estimate_repeat_count_from_product_length,
    estimate_repeat_count_from_spanning_read,
    expected_nonrepeat_bp,
    legacy_round_repeat_count,
    normalize_allele,
    repeat_unit_length,
)

__all__ = [
    "assembly_equivalent_product_allele",
    "estimate_repeat_count_from_inner_length",
    "estimate_repeat_count_from_product_length",
    "estimate_repeat_count_from_spanning_read",
    "expected_nonrepeat_bp",
    "legacy_round_repeat_count",
    "normalize_allele",
    "repeat_unit_length",
]