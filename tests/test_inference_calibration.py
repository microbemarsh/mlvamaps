"""Assembly units survive direct evidence, fallbacks, and final rounding."""
from dataclasses import replace

import pytest

from mlvamaps.alignment_evidence import CandidateEvidence
from mlvamaps.allele_inference import InferenceThresholds, infer_alleles, estimate_unresolved_calls, finalize_repeat_count_calls
from mlvamaps.calling import legacy_round_repeat_count
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.models import Locus


def test_complete_product_precedes_a_conflicting_partial_measurement_of_same_molecule():
    locus = Locus("L", repeat_unit_length_bp=4)
    complete = CandidateEvidence("L", 7, "m", measured_repeat_count=7, raw_repeat_count=7,
                                 direct_product_measurement=True, repeat_sigma=.08)
    partial = replace(complete, measured_repeat_count=5, raw_repeat_count=5, direct_product_measurement=False)
    calls, molecules = infer_alleles([complete, partial], [locus], [], "s", "illumina",
                                    require_measurement=True, sample_mode="isolate")
    assert calls[0]["repeat_count"] == 7 and calls[0]["molecule_support"] == 1
    assert calls[0]["presence_support"] == 1 and molecules == {("L", "m"): 7}


def test_primary_support_uses_current_rounding_of_raw_measurements():
    locus = Locus("L", repeat_unit_length_bp=4)
    evidence = [CandidateEvidence("L", 5.5, str(i), measured_repeat_count=5.5, raw_repeat_count=5.25,
                 repeat_sigma=.08, direct_product_measurement=True) for i in range(3)]
    calls, molecules = infer_alleles(evidence, [locus], [], "s", "hifi", require_measurement=True,
                                    sample_mode="isolate", round_tolerance=.4)
    finalize_repeat_count_calls(calls, evidence, [locus], [], round_tolerance=.4)
    assert calls[0]["repeat_count"] == 5 and calls[0]["primary_molecule_support"] == 3
    assert set(molecules.values()) == {5}


@pytest.mark.parametrize("tolerance", [0, .25, .4])
def test_spanning_alignment_adjusts_raw_reference_before_assembly_rounding(tolerance):
    locus = Locus("L", repeat_unit_length_bp=4, expected_product_size_bp=100, nominal_repeat_units=5)
    reference = CandidateContext("c", "L", "A" * 101, legacy_round_repeat_count(5.25, tolerance),
                                 4, 40, 61, "", "")
    evidence = [CandidateEvidence("L", reference.repeat_count, "m", candidate_id="c", alignment_score=100,
                 alignment_identity=1, query_coverage=1,
                 metadata={"aligned_flank_bases": 40, "spanning_repeat_indel_delta": .25})]
    calls, _ = infer_alleles(evidence, [locus], [reference], "s", "illumina", require_measurement=True)
    estimate_unresolved_calls(calls, evidence, [reference], loci=[locus], round_tolerance=tolerance)
    assert calls[0]["repeat_count_raw"] == 5.5
    assert calls[0]["repeat_count"] == legacy_round_repeat_count(5.5, tolerance)
    assert calls[0]["estimate_support"] == 1


@pytest.mark.parametrize("tolerance", [0, .25, .4])
def test_conflicting_and_partial_fallbacks_keep_raw_counts_and_use_final_rounding(tolerance):
    locus = Locus("L", repeat_unit_length_bp=4, expected_product_size_bp=100, nominal_repeat_units=5)
    conflicting = [CandidateEvidence("L", legacy_round_repeat_count(raw, tolerance), molecule,
                   measured_repeat_count=legacy_round_repeat_count(raw, tolerance), raw_repeat_count=raw,
                   repeat_sigma=.08, direct_product_measurement=True)
                   for molecule, values in (("a", (5.25, 8)), ("b", (5.25, 9))) for raw in values]
    partial = [CandidateEvidence("L", None, "m", metadata={"repeat_lower_bound": 7.8889})]
    for evidence, raw in ((conflicting, 5.25), (partial, 7.8889)):
        calls, _ = infer_alleles(evidence, [locus], [], "s", "illumina", require_measurement=True)
        estimate_unresolved_calls(calls, evidence, [], loci=[locus], round_tolerance=tolerance)
        finalize_repeat_count_calls(calls, evidence, [locus], [], round_tolerance=tolerance)
        assert calls[0]["repeat_count_raw"] == raw
        assert calls[0]["repeat_count"] == legacy_round_repeat_count(raw, tolerance)
        assert not calls[0]["prior_only"]


def test_mixture_at_minimum_observed_support_keeps_both_alleles_and_actual_depths():
    locus = Locus("L", repeat_unit_length_bp=4)
    evidence = [CandidateEvidence("L", allele, str(index), measured_repeat_count=allele,
                 raw_repeat_count=allele, repeat_sigma=.175, direct_product_measurement=True)
                for index, allele in enumerate([5] * 6 + [8] * 3)]
    evidence += [CandidateEvidence("L", None, f"partial{index}") for index in range(5)]
    calls, _ = infer_alleles(evidence, [locus], [], "s", "illumina", require_measurement=True,
                            sample_mode="metagenome")
    finalize_repeat_count_calls(calls, evidence, [locus], [])
    assert calls[0]["status"] == "mixed" and calls[0]["repeat_count"] == 5
    assert calls[0]["secondary_repeat"] == 8
    assert calls[0]["presence_support"] == 14 and calls[0]["molecule_support"] == 9
    assert calls[0]["primary_molecule_support"] == 6 and calls[0]["secondary_molecule_support"] == 3


def test_supported_consensus_keeps_low_likelihood_and_can_export_matching_sample_sequence():
    from mlvamaps.mapped_reconstruction import apply_reconstructed_calls
    from mlvamaps.unified_fastq import taxonomic_query_sequences
    locus = Locus("L", repeat_unit_length_bp=4)
    evidence = [CandidateEvidence("L", raw, str(i), measured_repeat_count=raw, raw_repeat_count=raw,
                 repeat_sigma=.5, direct_product_measurement=True) for i, raw in enumerate([5, 5.5] * 2)]
    limits = InferenceThresholds(minimum_probability=1)
    calls, _ = infer_alleles(evidence, [locus], [], "s", "hifi", limits, require_measurement=True,
        sample_mode="isolate", consensus_measurements={"L": {"raw_repeat_count": 5.25, "supporting_reads": 4}})
    before = dict(calls[0])
    assert before["consensus_status"] == "resolved_ambiguity" and before["confidence"] == "low"
    apply_reconstructed_calls(calls, [dict(locus_id="L", status="RECONSTRUCTED", method="mapped_flank_consensus",
        repeat_count=5.5, raw_repeat_count=5.25, reconstruction_molecules=4, reconstruction_fraction=1.)], limits)
    assert calls[0] == before
    assert taxonomic_query_sequences(calls, {"L": ("ACGT", 5.5)}) == {"L": "ACGT"}
    assert taxonomic_query_sequences(calls, {"L": ("ACGT", 5)}) == {}


def test_rejected_consensus_preserves_molecule_likelihood_and_count():
    from mlvamaps.mapped_reconstruction import apply_reconstructed_calls
    from mlvamaps.unified_fastq import taxonomic_query_sequences
    locus = Locus("L", repeat_unit_length_bp=4)
    evidence = [CandidateEvidence("L", raw, str(i), measured_repeat_count=raw, raw_repeat_count=raw,
                 repeat_sigma=.5, direct_product_measurement=True) for i, raw in enumerate([5, 5.5] * 2)]
    limits = InferenceThresholds(minimum_probability=1)
    calls, _ = infer_alleles(evidence, [locus], [], "s", "hifi", limits, require_measurement=True,
        sample_mode="isolate", consensus_measurements={"L": {"raw_repeat_count": 8, "supporting_reads": 4}})
    before = dict(calls[0])
    assert before["consensus_status"] == "disagrees" and before["status"] == "ambiguous"
    assert before["repeat_count"] == 5.5 and before["molecule_support"] == 4
    apply_reconstructed_calls(calls, [dict(locus_id="L", status="RECONSTRUCTED", method="mapped_flank_consensus",
        repeat_count=8, raw_repeat_count=8, reconstruction_molecules=4, reconstruction_fraction=1.)], limits)
    assert calls[0] == before
    assert taxonomic_query_sequences(calls, {"L": ("ACGT", 8)}) == {}
