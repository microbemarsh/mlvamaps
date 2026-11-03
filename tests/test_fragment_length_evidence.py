"""Calibrated paired length observations retain boundary and sampling geometry."""
from dataclasses import replace
import random

import pytest

from mlvamaps.alignment_evidence import CandidateAlignment
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.sequence import revcomp
from mlvamaps.short_read_evidence import _partial_pair_geometry, extract_short_read_evidence
from mlvamaps.short_read_fragments import RecruitedFragment


def partial_pair(left_bases=30, right_bases=30):
    rng = random.Random(731)
    forward, left, right, reverse = ["".join(rng.choices("ACGT", k=n)) for n in (20, 50, 50, 20)]
    locus = Locus("L", forward_primer=forward, reverse_primer=reverse,
        left_flank_sequence=left, right_flank_sequence=right, repeat_motif="ATGC",
        repeat_unit_length_bp=4, expected_product_size_bp=156, nominal_repeat_units=4)
    product = forward + left + "ATGC" * 40 + right + revcomp(reverse)
    reference = forward + left + "ATGC" * 4 + right + revcomp(reverse)
    context = CandidateContext("ref", "L", reference, 4, 4, 70, 86, left, right)
    first = product[70 - left_bases:170 - left_bases]
    end = 230 + right_bases
    second = revcomp(product[end - 100:end])
    pair = ReadPair("m", ReadRecord("m/1", first, "I" * len(first)),
                    ReadRecord("m/2", second, "I" * len(second)))
    return locus, context, pair, product


@pytest.mark.parametrize("left,right", [(30, 30), (40, 15), (15, 40)])
def test_partial_boundaries_retain_calibrated_length_and_selection_window(left, right):
    locus, context, pair, _ = partial_pair(left, right)
    result = _partial_pair_geometry(pair, locus, [context], (140,), 160 + left + right, 12)
    assert result["pair_repeat_count"] == 40 and result["pair_repeat_sigma"] == 3
    assert result["pair_nonrepeat_bases"] == left + right
    assert result["pair_geometry_source"] == "partial_repeat_boundaries"
    assert result["fragment_length_scheme"] == "anchor_union_v1"
    assert result["fragment_length_class"] == "boundary_pair"
    observed = right if left >= 21 else left
    assert result["pair_selection_condition"] == ("left_anchor" if left >= 21 else "right_anchor")
    assert result["pair_selection_low_delta"] == (12 - observed) / 4
    assert result["pair_selection_high_delta"] == ((70 if left >= 21 else 20) - observed) / 4
    assert result["pair_repeat_lower_bound"] == max(100 - left, 100 - right) / 4
    swapped = ReadPair("m", pair.read2, pair.read1)
    assert _partial_pair_geometry(swapped, locus, [context], (140,), 160 + left + right, 12) == result


def test_partial_pair_geometry_needs_strong_unique_anchor_and_unambiguous_calibration():
    locus, context, pair, _ = partial_pair(15, 15)
    assert not _partial_pair_geometry(pair, locus, [context], (140,), 190, 12)
    locus, context, pair, _ = partial_pair()
    assert not _partial_pair_geometry(pair, locus, [context], (140, 144), 220, 12)
    assert not _partial_pair_geometry(pair, locus, [context], (140,), None, None)
    assert not _partial_pair_geometry(replace(pair, read2=None), locus, [context], (140,), 220, 12)
    pure = replace(pair, read2=ReadRecord("m/2", "ATGC" * 25, "I" * 100))
    censored = _partial_pair_geometry(pure, locus, [context], (140,), 220, 12)
    assert censored["fragment_length_class"] == "repeat_mate"
    assert censored["fragment_anchor_bp"] == 30 and censored["fragment_mate_length_bp"] == 100
    assert "pair_repeat_count" not in censored and "pair_repeat_lower_bound" not in censored
    damaged = replace(pair, read1=replace(pair.read1, sequence=pair.read1.sequence[:-1] + "N"))
    assert not _partial_pair_geometry(damaged, locus, [context], (140,), 220, 12)
    assert not _partial_pair_geometry(pair, locus, [context, replace(context, repeat_start=69)], (140,), 220, 12)


def test_partial_pair_extraction_adds_an_estimate_without_direct_measurement():
    locus, context, pair, _ = partial_pair(40, 15)
    fragments = {"L": {"m": RecruitedFragment(pair, True, True, False)}}
    rows = extract_short_read_evidence([], [context], [locus], recruited_fragments=fragments,
                                      insert_mean=215, insert_sd=12)
    assert len(rows) == 1
    assert rows[0].metadata["pair_repeat_count"] == 40
    assert rows[0].measured_repeat_count is None and not rows[0].direct_product_measurement
    assert rows[0].molecule_id == "m"


def test_complete_products_keep_priority_over_paired_length_estimates():
    locus, context, pair, product = partial_pair()
    full = ReadPair("m", ReadRecord("m/1", product, "I" * len(product)),
                    ReadRecord("m/2", revcomp(product), "I" * len(product)))
    rows = extract_short_read_evidence([], [context], [locus],
        recruited_fragments={"L": {"m": RecruitedFragment(full, True, True, False)}},
        insert_mean=215, insert_sd=12)
    assert rows[0].measured_repeat_count == 40 and rows[0].direct_product_measurement
    assert "pair_repeat_count" not in rows[0].metadata


def test_repeat_rich_mate_is_statistical_evidence_without_an_observed_upper_boundary():
    locus, context, pair, _ = partial_pair(30, 30)
    pair = replace(pair, read2=ReadRecord("m/2", "ATGC" * 25, "I" * 100))
    rows = extract_short_read_evidence([], [context], [locus],
        recruited_fragments={"L": {"m": RecruitedFragment(pair, True, False, True)}},
        insert_mean=220, insert_sd=12)
    metadata = rows[0].metadata
    assert metadata["repeat_rich_mate"] and metadata["fragment_length_class"] == "repeat_mate"
    assert metadata["repeat_lower_bound"] == 17.5  # Only 70 anchored repeat bases.
    assert "pair_repeat_count" not in metadata and rows[0].measured_repeat_count is None


def test_whole_flank_pairs_record_the_accepted_opposite_flank_window():
    locus = Locus("L", repeat_unit_length_bp=4, expected_product_size_bp=220, nominal_repeat_units=5)
    context = CandidateContext("ref", "L", "A" * 220, 5, 4, 100, 120, "", "")
    first = CandidateAlignment("m", "m/1", 1, "L", "ref", 5, "R", 90, 60, 1, 1, .2,
        0, 40, 30, 70, "40M", "", True, False, False, query_sequence="A" * 40, query_quality="I" * 40)
    second = replace(first, read_id="m/2", mate=2, reference_start=150, reference_end=190, reverse=True)
    rows = extract_short_read_evidence([first, second], [context], [locus], insert_mean=168, insert_sd=8)
    metadata = rows[0].metadata
    assert metadata["pair_repeat_count"] == 7
    assert metadata["pair_selection_low_delta"] == -7.5
    assert metadata["pair_selection_high_delta"] == 7.5
    assert metadata["pair_geometry_source"] == "opposite_flanks"
