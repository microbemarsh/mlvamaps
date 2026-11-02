from dataclasses import replace

from mlvamaps.alignment_evidence import CandidateAlignment, CandidateEvidence
from mlvamaps.allele_inference import InferenceThresholds, infer_alleles, estimate_unresolved_calls
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.io import ReadSpool
from mlvamaps.long_read_evidence import flank_error_model, refine_repeat_sequence
from mlvamaps.models import Locus, ReadRecord
from mlvamaps.qc import filter_reads


def test_spooled_qc_uses_one_stream_and_preserves_quality_and_random_access():
    records = (ReadRecord(str(i), "ACGT" * 20, "I" * 80) for i in range(100))
    with ReadSpool() as spool:
        kept, qc = filter_reads(records, retained=spool)
        assert kept is spool and len(spool) == 100 and qc[0]["value"] == 100
        assert next(iter(spool)) == spool.get("0")
        assert spool.get("99").quality == "I" * 80
        assert all(isinstance(offset, int) for offset in spool.offsets.values())
    assert spool.handle.closed


def test_resource_limit_retains_evidence_without_a_false_point_estimate():
    evidence = [CandidateEvidence("L", 4, f"m{i}", measured_repeat_count=4,
                                  raw_repeat_count=4, repeat_sigma=100) for i in range(3)]
    calls, molecules = infer_alleles(evidence, [Locus("L", repeat_unit_length_bp=4)], [], "s", "ont",
                                   InferenceThresholds(maximum_candidates=20), require_measurement=True)
    estimate_unresolved_calls(calls, evidence, [])
    assert calls[0]["repeat_count"] == "" and len(molecules) == 3
    assert calls[0]["computational_fallback"] == "candidate_limit"
    assert calls[0]["repeat_interval_high"] >= 304


def test_isolate_bimodality_is_retained_without_calling_an_unobserved_middle():
    evidence = [CandidateEvidence("L", 4, f"m{i}", measured_repeat_count=allele,
                                  raw_repeat_count=allele, repeat_sigma=.1)
                for i, allele in enumerate([3] * 12 + [9] * 8)]
    calls, _ = infer_alleles(evidence, [Locus("L", repeat_unit_length_bp=4)], [], "s", "ont",
                           require_measurement=True, sample_mode="isolate")
    estimate_unresolved_calls(calls, evidence, [])
    assert calls[0]["repeat_count"] == "" and calls[0]["possible_secondary"]
    assert calls[0]["secondary_repeat"] == 9


def test_native_motif_refinement_is_bounded_and_reports_limits():
    limits = InferenceThresholds()
    scores, reason = refine_repeat_sequence("ATGC" * 8, "ATGC", limits)
    assert not reason and max(scores, key=scores.get) == 8
    scores, reason = refine_repeat_sequence("ATGC" * 3 + "A" + "ATGC" * 5, "ATGC", limits)
    assert not reason and max(scores, key=scores.get) == 8
    assert refine_repeat_sequence("A" * 10001, "A", limits)[1] == "alignment_resource_limit"
    assert refine_repeat_sequence("ATGC" * 8, "ATGC", replace(limits, maximum_candidates=2))[1] == "candidate_limit"


def test_error_model_counts_independent_flanks_and_excludes_repeat_indels():
    context = CandidateContext("c", "L", "A" * 240, 10, 4, 100, 140, "", "")
    prototype = CandidateAlignment("m", "m", None, "L", "c", 10, "ref", 240, 60, .99,
        1, 1, 0, 240, 0, 240, "100M4I140M", ":100+aaaa:140", True, False, False)
    rows = [replace(prototype, molecule_id=str(i)) for i in range(20)]
    model = flank_error_model(rows, [context])
    assert model["source"] == "pooled_flanks" and model["flank_bases"] == 4000
    assert model["insertion"] == .001
    assert flank_error_model(rows * 3, [context]) == model
    assert flank_error_model(rows[:2], [context])["source"] == "conservative_prior"


def test_rescued_products_share_ont_uncertainty_and_refinement_limits():
    from mlvamaps.long_read_evidence import long_read_measurement_metadata
    from mlvamaps.locus_measurement import measure_locus_product
    locus = Locus("L", left_flank_sequence="GACTGACCTGATCGACGTTAGCTA",
                  right_flank_sequence="CCGTAACCGGTTAGGCTACTGAC", repeat_motif="ATGC", repeat_unit_length_bp=4)
    sequence = locus.left_flank_sequence + "ATGC" * 100 + locus.right_flank_sequence
    measurement = measure_locus_product(sequence, locus, "I" * len(sequence))
    model = flank_error_model([], [])
    metadata = long_read_measurement_metadata(measurement, locus, model, InferenceThresholds())
    assert metadata["error_model"] == model
    assert metadata["repeat_sigma"] > measurement.metadata["repeat_sigma"]
    noisy = replace(measurement, repeat_motif_identity=.5)
    limited = long_read_measurement_metadata(noisy, locus, model, InferenceThresholds(), refine=False)
    assert limited["refinement_fallback"] == "refinement_read_limit"


def test_legacy_consensus_and_variant_mixtures_also_obey_resource_limits():
    from mlvamaps.local_assembly import sequence_consensus
    from mlvamaps.mixture import estimate_variant_mixtures
    limits = InferenceThresholds(maximum_alignment_cells=20)
    assert sequence_consensus(["A" * 10000] * 3, limits) == "A" * 10000
    assert sequence_consensus(["ACGT" * 20, "ACGA" * 20], limits) == ""
    rows = [dict(sample_id="s", locus_id="L", variant_id=str(i), support_reads=4,
                 repeat_count=i, representative_sequence=sequence)
            for i, sequence in enumerate(("ACGT" * 20, "ACGA" * 20))]
    mixtures = estimate_variant_mixtures(rows, limits=limits)
    assert all(row["computational_fallback"] == "observed_counts_only" for row in mixtures)
    assert all(row["meaningful"] == "no" for row in mixtures)


def test_bounded_anchor_search_matches_searching_the_same_window():
    from mlvamaps.locus_measurement import find_anchor
    pattern = "ACGTACGACAGTCGATGCAA"
    for index in range(len(pattern)):
        window = "TT" + pattern[:index] + "T" + pattern[index + 1:] + "GG"
        sequence = "CGATCG" + window + "AACCGG"
        expected = find_anchor(pattern, window, 2)
        bounded = find_anchor(pattern, sequence, 2, 6, len(sequence) - 6)
        assert expected is not None and bounded is not None
        assert (bounded.start - 6, bounded.end - 6, bounded.edit_distance) == (expected.start, expected.end, expected.edit_distance)
