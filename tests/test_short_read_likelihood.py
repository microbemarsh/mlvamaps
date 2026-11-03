"""Conditional fragment geometry improves estimates without creating measurements."""
from dataclasses import replace
import math
from statistics import NormalDist, mean

import pytest

from mlvamaps.alignment_evidence import CandidateEvidence
from mlvamaps.allele_inference import InferenceThresholds, infer_alleles
from mlvamaps.models import Locus
from mlvamaps.short_read_likelihood import _log_normal_interval, infer_short_read_length

LOCUS = Locus("L", repeat_unit_length_bp=4)


def selected_pairs(count=40, molecules=120):
    """Quantiles from an independently calibrated, selected fragment library."""
    library = NormalDist(260, 30)
    low, high = count * 4 + 90, count * 4 + 150
    first, last = library.cdf(low), library.cdf(high)
    return [CandidateEvidence("L", None, str(index), metadata={
        "pair_repeat_count": count + (library.mean - length) / 4,
        "pair_repeat_sigma": library.stdev / 4, "pair_alignment_score": 100,
        "pair_selection_low_delta": (low - length) / 4,
        "pair_selection_high_delta": (high - length) / 4,
    }) for index in range(molecules)
        for length in [library.inv_cdf(first + (last - first) * (index + .5) / molecules)]]


def union_fragments(count=40, molecules=2000):
    """Uniform normal quantiles, retaining exactly the documented class union."""
    library, anchor, mate = NormalDist(245, 30), 45, 100
    low, high, repeat = 12, 80, count * 4
    rows = []
    for index in range(molecules):
        size = round(library.inv_cdf((index + .5) / molecules))
        if anchor + mate <= size <= anchor + repeat:
            kind = "repeat_mate"
        elif anchor + repeat + low <= size <= anchor + repeat + high:
            kind = "boundary_pair"
        else:
            continue
        metadata = dict(fragment_length_scheme="anchor_union_v1", fragment_length_class=kind,
            fragment_anchor_bp=anchor, fragment_mate_length_bp=mate, fragment_boundary_min_bp=low,
            fragment_boundary_max_bp=high, fragment_insert_mean_bp=library.mean,
            fragment_insert_sd_bp=library.stdev, fragment_allele_offset=0,
            fragment_repeat_unit_bp=4, pair_alignment_score=100, repeat_lower_bound=13.75)
        if kind == "boundary_pair":
            metadata["fragment_opposite_anchor_bp"] = size - anchor - repeat
        rows.append(CandidateEvidence("L", None, str(index), metadata=metadata))
    return rows


@pytest.mark.parametrize("count", [40, 40.5])
def test_selection_conditioning_recovers_novel_lengths_without_measured_votes(count):
    rows = selected_pairs(count)
    assert abs(mean(row.metadata["pair_repeat_count"] for row in rows) - count) > 3
    calls, molecules = infer_alleles(rows, [LOCUS], [], "s", "illumina", require_measurement=True,
                                   sample_mode="isolate")
    result = calls[0]
    assert result["repeat_count_raw"] == result["repeat_count"] == count
    assert result["estimation_method"] == "joint_short_read_length"
    assert result["status"] == "estimated" and result["confidence"] == "low"
    assert result["finite_length_support"] and not result["upper_unbounded"]
    assert result["repeat_interval_low"] <= count <= result["repeat_interval_high"]
    assert result["estimate_support"] == 120 and result["molecule_support"] == 0 and molecules == {}
    assert result["dominant_repeat"] == count and result["allele_fractions"] == f"{count}:1.00000000"


def test_joint_boundary_and_repeat_mate_classes_recover_known_length():
    rows = union_fragments(molecules=600)
    result = infer_short_read_length(rows, LOCUS, InferenceThresholds())
    assert result["repeat_count"] == 40
    assert result["length_null_fraction"] == 0
    assert result["repeat_rich_mate_support"] > 0 and result["paired_span_support"] > 0
    assert result["estimate_support"] == len(rows)
    assert result == infer_short_read_length(rows * 4, LOCUS, InferenceThresholds())
    pure = [row for row in rows if row.metadata["fragment_length_class"] == "repeat_mate"]
    only = infer_short_read_length(pure, LOCUS, InferenceThresholds())
    assert only["upper_unbounded"] and not only["finite_length_support"]
    assert "repeat_count" not in only and only["repeat_lower_bound"] == 13.75
    outlier = replace(pure[0], molecule_id="unrelated", metadata={**pure[0].metadata,
                      "fragment_mate_length_bp": 500})
    protected = infer_short_read_length(rows + [outlier], LOCUS, InferenceThresholds())
    assert protected["repeat_count"] == 40 and protected["repeat_lower_bound"] == 13.75
    calls, molecules = infer_alleles(rows, [LOCUS], [], "s", "illumina", require_measurement=True,
                                   sample_mode="isolate")
    assert calls[0]["repeat_count"] == 40 and calls[0]["molecule_support"] == 0 and not molecules


def test_reference_copies_mates_and_repeated_partial_bounds_do_not_multiply_votes():
    rows = selected_pairs()
    bound = CandidateEvidence("L", None, "bound", metadata={"repeat_lower_bound": 41.25})
    original = infer_short_read_length(rows + [bound], LOCUS, InferenceThresholds())
    duplicate = infer_short_read_length(list(reversed(rows * 4)) + [bound] * 30, LOCUS, InferenceThresholds())
    assert original == duplicate
    assert original["repeat_count_raw"] == original["repeat_lower_bound"] == 41.25
    assert original["repeat_count"] == 41.5
    assert original["estimate_support"] == 120 and original["partial_bound_support"] == 1


def test_ambiguous_locus_molecules_and_direct_measurements_keep_precedence():
    rows = selected_pairs()
    ambiguous = CandidateEvidence("L", None, "ambiguous", quality_flags=("ambiguous_locus",),
                                  metadata={"pair_repeat_count": math.nan, "pair_repeat_sigma": 0,
                                            "repeat_lower_bound": 10000})
    expected = infer_short_read_length(rows, LOCUS, InferenceThresholds())
    assert infer_short_read_length(rows + [ambiguous], LOCUS, InferenceThresholds()) == expected
    direct = CandidateEvidence("L", 5, "direct", measured_repeat_count=5, direct_product_measurement=True)
    assert infer_short_read_length(rows + [direct], LOCUS, InferenceThresholds()) == {}


def test_bounds_only_and_saturated_pairs_do_not_claim_a_finite_upper_length():
    bound = CandidateEvidence("L", None, "m", metadata={"repeat_lower_bound": 17.5})
    result = infer_short_read_length([bound], LOCUS, InferenceThresholds())
    assert result["length_model_status"] == "lower_bound_only" and result["upper_unbounded"]
    assert "repeat_count" not in result and not result["finite_length_support"]
    tail = CandidateEvidence("L", None, "m", metadata={"pair_repeat_count": 0, "pair_repeat_sigma": 2,
                             "pair_selection_low_delta": 0, "pair_selection_high_delta": 5})
    result = infer_short_read_length([tail], LOCUS, InferenceThresholds())
    assert result["upper_unbounded"] and not result["finite_length_support"]
    assert result["repeat_interval_high"] == ""


def test_legacy_pairs_and_metagenomes_do_not_gain_an_unjustified_single_allele_fit():
    rows = selected_pairs()
    legacy = [replace(row, metadata={key: value for key, value in row.metadata.items()
                                    if not key.startswith("pair_selection")}) for row in rows]
    assert infer_short_read_length(legacy, LOCUS, InferenceThresholds()) == {}
    assert infer_short_read_length(rows, LOCUS, InferenceThresholds(), sample_mode="metagenome") == {
        "length_model_status": "unsupported_sample_mode"}


def test_repeat_mate_tail_beyond_boundary_grid_is_not_a_finite_length():
    metadata = dict(fragment_length_scheme="anchor_union_v1", fragment_length_class="repeat_mate",
        fragment_anchor_bp=45, fragment_mate_length_bp=150, fragment_boundary_min_bp=12,
        fragment_boundary_max_bp=146, fragment_insert_mean_bp=245, fragment_insert_sd_bp=4,
        fragment_allele_offset=0, fragment_repeat_unit_bp=4, repeat_lower_bound=26.25)
    rows = [CandidateEvidence("L", None, str(index), metadata=metadata) for index in range(10)]
    rows.append(CandidateEvidence("L", None, "boundary", metadata={**metadata,
        "fragment_length_class": "boundary_pair", "fragment_opposite_anchor_bp": 100}))
    result = infer_short_read_length(rows, LOCUS, InferenceThresholds())
    # Boundary centers alone stop near 34 repeats, before the pure-mate event
    # becomes feasible. The long-array plateau must remain upper-unbounded.
    assert result["candidate_grid_high"] >= 58
    assert result["repeat_count_raw"] > 50 and result["length_null_fraction"] == .01
    assert result["upper_unbounded"] and not result["finite_length_support"]
    assert result["repeat_interval_high"] == ""


def test_candidate_and_cell_limits_are_explicit_and_fallback_respects_bounds():
    rows = selected_pairs()
    limited = infer_short_read_length(rows, LOCUS, InferenceThresholds(maximum_candidates=8))
    assert limited["computational_fallback"] == "short_read_length_grid_limit"
    assert len(limited["candidate_distribution"].split(";")) <= 8 and limited["estimate_support"] == 120
    bound = CandidateEvidence("L", None, "bound", metadata={"repeat_lower_bound": 50})
    skipped = infer_short_read_length(rows + [bound], LOCUS, InferenceThresholds(maximum_likelihood_cells=2))
    assert skipped["length_model_status"] == "insufficient_budget"
    assert skipped["upper_unbounded"] and skipped["repeat_count"] == 50
    assert skipped["best_probability"] == 0 and skipped["repeat_interval_high"] == ""


@pytest.mark.parametrize("updates", [{"pair_repeat_sigma": 0}, {"pair_repeat_count": math.nan},
    {"pair_selection_low_delta": 1, "pair_selection_high_delta": 2},
    {"pair_selection_low_delta": 0, "pair_selection_high_delta": 0},
    {"pair_selection_low_delta": None}])
def test_invalid_pair_geometry_is_rejected(updates):
    row = selected_pairs(molecules=1)[0]
    with pytest.raises(ValueError, match="paired-fragment"):
        infer_short_read_length([replace(row, metadata={**row.metadata, **updates})], LOCUS, InferenceThresholds())


def test_conditional_normal_probability_is_stable_in_far_library_tails():
    assert _log_normal_interval(-2, 3) == pytest.approx(math.log(NormalDist().cdf(3) - NormalDist().cdf(-2)))
    value = _log_normal_interval(40, 40.5)
    assert math.isfinite(value) and value < -800
    assert _log_normal_interval(-40.5, -40) == pytest.approx(value)
