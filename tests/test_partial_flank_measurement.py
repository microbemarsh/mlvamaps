from dataclasses import replace
import random

import pytest

from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.locus_measurement import measure_locus_product
from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.sequence import revcomp
from mlvamaps.short_read_evidence import extract_short_read_evidence
from mlvamaps.short_read_fragments import RecruitedFragment


def _locus():
    rng = random.Random(9821)
    forward, left, right, reverse = [
        "".join(rng.choices("ACGT", k=size)) for size in (20, 50, 50, 20)
    ]
    return Locus("partial", forward_primer=forward, reverse_primer=reverse,
                 left_flank_sequence=left, right_flank_sequence=right,
                 repeat_motif="GATC", repeat_unit_length_bp=4,
                 nominal_repeat_units=4, expected_product_size_bp=156,
                 expected_min_repeats=2, expected_max_repeats=8)


def _substitute(sequence, index):
    base = "A" if sequence[index] != "A" else "C"
    return sequence[:index] + base + sequence[index + 1:]


@pytest.mark.parametrize("count", [0, 7, 19])
@pytest.mark.parametrize("truncation", [(19, 0), (0, 19), (19, 19)])
def test_truncated_flanks_measure_novel_and_zero_counts_despite_anchor_errors(count, truncation):
    locus = _locus()
    left = _substitute(locus.left_flank_sequence, 40)
    right = _substitute(locus.right_flank_sequence, 10)
    sequence = left[truncation[0]:] + locus.repeat_motif * count + right[:50 - truncation[1]]
    measurement = measure_locus_product(sequence, locus, "?" * len(sequence),
                                        reference_nonrepeat_bp=(140,))
    assert measurement.status == "REPEAT_INFORMATIVE"
    assert measurement.raw_repeat_count == measurement.called_allele == count
    assert measurement.repeat_sequence == locus.repeat_motif * count
    assert measurement.product_sequence == ""
    assert not measurement.metadata["boundary_ambiguous"]


@pytest.mark.parametrize("edit", ["insert", "delete"])
def test_partial_anchor_indels_do_not_change_the_repeat_interval(edit):
    locus = _locus()
    left = locus.left_flank_sequence[-31:]
    left = left[:20] + ("A" + left[20:] if edit == "insert" else left[21:])
    sequence = left + locus.repeat_motif * 19 + locus.right_flank_sequence[:31]
    measurement = measure_locus_product(sequence, locus, reference_nonrepeat_bp=(140,))
    assert measurement.status == "REPEAT_INFORMATIVE"
    assert measurement.raw_repeat_count == 19


@pytest.mark.parametrize("calibration", [None, (), (140, 144)])
def test_partial_flanks_require_unique_reference_calibration(calibration):
    locus = _locus()
    sequence = locus.left_flank_sequence[-31:] + locus.repeat_motif * 7 + locus.right_flank_sequence[:31]
    measurement = measure_locus_product(sequence, locus, reference_nonrepeat_bp=calibration)
    assert measurement.raw_repeat_count is None


@pytest.mark.parametrize("primer", ["forward", "reverse"])
def test_surviving_primer_distinguishes_a_flank_deletion_from_read_truncation(primer):
    locus = _locus()
    sequence = (locus.left_flank_sequence[-21:] + locus.repeat_motif * 19
                + locus.right_flank_sequence[:21])
    sequence = locus.forward_primer + sequence if primer == "forward" else sequence + revcomp(locus.reverse_primer)
    measurement = measure_locus_product(sequence, locus, reference_nonrepeat_bp=(140,))
    assert measurement.raw_repeat_count is None


@pytest.mark.parametrize("problem", ["duplicate", "low_complexity", "interior_failure", "reversed"])
def test_partial_flanks_reject_ambiguous_weak_or_nonterminal_anchors(problem):
    locus = _locus()
    left, right = locus.left_flank_sequence[-21:], locus.right_flank_sequence[:21]
    if problem == "duplicate":
        left = left + "C" + left
    elif problem == "low_complexity":
        locus = replace(locus, left_flank_sequence="A" * 50)
        left = "A" * 21
    elif problem == "interior_failure":
        left = locus.left_flank_sequence
        for index in (2, 6, 10, 14, 18):
            left = _substitute(left, index)
    sequence = left + locus.repeat_motif * 7 + right
    if problem == "reversed":
        sequence = right + locus.repeat_motif * 7 + left
    measurement = measure_locus_product(sequence, locus, reference_nonrepeat_bp=(140,))
    assert measurement.raw_repeat_count is None


@pytest.mark.parametrize("reverse", [False, True])
def test_short_read_precheck_retains_partial_flank_measurements(reverse):
    locus = _locus()
    reference = (locus.forward_primer + locus.left_flank_sequence + locus.repeat_motif * 4
                 + locus.right_flank_sequence + revcomp(locus.reverse_primer))
    context = CandidateContext("reference", locus.locus_id, reference, 4, 4, 70, 86,
                               locus.left_flank_sequence, locus.right_flank_sequence)
    sequence = (_substitute(locus.left_flank_sequence[-31:], 20) + locus.repeat_motif * 19
                + _substitute(locus.right_flank_sequence[:31], 10))
    if reverse:
        sequence = revcomp(sequence)
    read = ReadRecord("molecule", sequence, "?" * len(sequence))
    fragments = {locus.locus_id: {"molecule": RecruitedFragment(
        ReadPair("molecule", read), True, True, False)}}
    evidence = extract_short_read_evidence([], [context], [locus], recruited_fragments=fragments)
    assert len(evidence) == 1
    assert evidence[0].molecule_id == "molecule"
    assert evidence[0].raw_repeat_count == evidence[0].measured_repeat_count == 19
    assert evidence[0].full_repeat_span
    assert not evidence[0].direct_product_measurement
