"""Independent molecule identity, ambiguity and coordinate regressions."""
from dataclasses import replace

import pysam
import pytest

from mlvamaps.alignment_evidence import CandidateAlignment, CandidateEvidence, evidence_row
from mlvamaps.allele_inference import infer_alleles, estimate_unresolved_calls
from mlvamaps.calling import gaussian_allele_probabilities
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.locus_measurement import reference_interval_to_query_interval, extract_reference_interval_from_original_read
from mlvamaps.mapped_reconstruction import recruit_mapping_fragments
from mlvamaps.minimap_mapping import parse_candidate_alignments
from mlvamaps.models import Locus
from mlvamaps.sequence import revcomp


def test_supplementary_and_secondary_restore_original_sequence_and_coordinates(tmp_path):
    context = CandidateContext("c", "L", "A" * 10, 1, 2, 4, 6, "", "")
    original = "ACGTTGCAAC"
    oriented = revcomp(original)
    path = tmp_path / "reads.sam"
    path.write_text("@SQ\tSN:c\tLN:10\n"
        f"m\t2064\tc\t7\t20\t6H4M\t*\t0\t0\t{oriented[6:]}\tIIII\tAS:i:4\n"
        "m\t272\tc\t1\t0\t10M\t*\t0\t0\t*\t*\tAS:i:5\n"
        f"m\t16\tc\t1\t60\t6M4S\t*\t0\t0\t{oriented}\tABCDEFGHIJ\tAS:i:6\n")
    rows = parse_candidate_alignments(path, [context])
    assert len(rows) == 3 and rows[0].supplementary and not rows[0].primary
    assert (rows[0].query_start, rows[0].query_end, rows[0].cigar) == (6, 10, "6S4M")
    assert all(row.query_sequence == original and row.query_quality == "JIHGFEDCBA" for row in rows)
    assert rows[0].query_sequence is rows[-1].query_sequence


def test_cross_locus_ties_survive_as_ambiguity_without_an_estimated_count():
    contexts = [CandidateContext(name, name, "A" * 80, 5, 4, 30, 50, "", "") for name in ("a", "b")]
    row = CandidateAlignment("m", "m", None, "a", "a", 5, "ref", 80, 0, 1, 1, 1,
                             0, 80, 0, 80, "80M", "", True, False, False, "A" * 80)
    ambiguous = {}
    assert recruit_mapping_fragments([row, replace(row, locus_id="b", candidate_id="b")],
                                    contexts, ambiguous_assignments=ambiguous) == {}
    assert ambiguous == {"m": ("a", "b")}
    evidence = [CandidateEvidence(locus, None, "m", candidate_loci=("a", "b"),
                                  quality_flags=("ambiguous_locus",)) for locus in ("a", "b")]
    calls, molecules = infer_alleles(evidence, [Locus("a"), Locus("b")], contexts, "s", "ont", require_measurement=True)
    estimate_unresolved_calls(calls, evidence, contexts)
    assert not molecules
    assert all(call["presence_support"] == 1 and call["repeat_count"] == "" for call in calls)
    assert evidence_row("s", evidence[0])["candidate_loci"] == "a;b"


def test_conflicting_duplicate_measurements_never_vote_by_reference_copies():
    rows = [CandidateEvidence("L", 3, "m", measured_repeat_count=3, direct_product_measurement=True)] * 100
    rows += [CandidateEvidence("L", 7, "m", measured_repeat_count=7, direct_product_measurement=True)]
    for require_measurement in (False, True):
        calls, molecules = infer_alleles(rows, [Locus("L")], [], "s", "hifi", require_measurement=require_measurement)
        assert not molecules and calls[0]["repeat_count"] == ""


def test_coordinate_projection_is_bounded_by_cigar_operations():
    alignment = pysam.AlignedSegment()
    alignment.query_sequence = "ACGT"
    alignment.reference_start = 0
    alignment.cigarstring = "2M100000000D2M"
    assert reference_interval_to_query_interval(alignment, 0, 100000004) == (0, 4)
    assert reference_interval_to_query_interval(alignment, 2, 5) is None
    assert reference_interval_to_query_interval(alignment, -1, 1) is None
    assert extract_reference_interval_from_original_read(alignment, "AC", "II", 0, 100000004) is None


def test_likelihood_underflow_preserves_tied_candidates():
    assert gaussian_allele_probabilities(1000, [999, 1001], .0001) == [.5, .5]
    assert gaussian_allele_probabilities(1000, [998, 999], .0001)[1] == pytest.approx(1)


def test_phase_shifted_expansion_recruits_without_counting_flank_indels_as_a_repeat():
    from mlvamaps.mapped_reconstruction import recruitment_identity
    context = CandidateContext("c", "L", "A" * 106, 4, 4, 45, 61, "", "")
    row = CandidateAlignment("m", "m", None, "L", "c", 4, "r", 60, 60, 106 / 126,
        1, 1, 0, 126, 0, 106, "44M20I62M", "", True, False, False, "A" * 126)
    assert recruitment_identity(row, context) == pytest.approx(1)
    assert "m" in recruit_mapping_fragments([row], [context])["L"]


def test_one_base_of_the_opposite_flank_is_not_spanning_evidence():
    from mlvamaps.short_read_evidence import extract_short_read_evidence
    context = CandidateContext("c", "L", "A" * 106, 4, 4, 45, 61, "", "")
    row = CandidateAlignment("m", "m", None, "L", "c", 4, "r", 62, 60, 1,
        1, .6, 0, 62, 0, 62, "62M", "", True, False, False, "A" * 62)
    evidence = extract_short_read_evidence([row], [context], [Locus("L")])
    assert not evidence[0].full_repeat_span
    assert evidence[0].metadata["spanning_repeat_indel_delta"] is None


def test_long_read_parser_retains_only_local_sequence_with_original_coordinate_offset(tmp_path):
    from mlvamaps.io import ReadSpool
    from mlvamaps.models import ReadRecord
    context = CandidateContext("c", "L", "ACGT" * 25, 5, 4, 40, 60, "", "")
    original = "C" * 10000 + context.sequence + "G" * 10000
    path = tmp_path / "reads.sam"
    path.write_text("@SQ\tSN:c\tLN:100\n" +
        f"m\t0\tc\t1\t60\t10000S100M10000S\t*\t0\t0\t{original}\t{'I' * len(original)}\n")
    with ReadSpool() as reads:
        reads.append(ReadRecord("m", original, "I" * len(original)))
        row = parse_candidate_alignments(path, [context], query_reads=reads)[0]
    assert row.query_offset == 9976 and row.query_start == 24
    assert len(row.query_sequence) == 148 and row.cigar == "24S100M24S"
