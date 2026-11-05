"""Primary FASTQ calls abstain when sample evidence does not resolve an allele."""

import csv
from copy import deepcopy
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from mlvamaps.alignment_evidence import CandidateEvidence
from mlvamaps.allele_inference import (
    CONFIDENCE_FIELDS,
    InferenceThresholds,
    finalize_repeat_count_calls,
    infer_alleles,
)
from mlvamaps.evidence_reporting import write_evidence_reports
from mlvamaps.mapped_reconstruction import apply_reconstructed_calls
from mlvamaps.models import Locus
from mlvamaps.myoga_export import build_call_matrix
from mlvamaps.pipeline import simple_call_rows_from_alleles
from mlvamaps.profile_matching import build_fingerprint
from mlvamaps.unified_fastq import common_calls_to_compatibility


LOCUS = Locus("L", repeat_unit_length_bp=4, nominal_repeat_units=4,
              expected_product_size_bp=200)
LIMITS = InferenceThresholds()
NOISY_MEASUREMENTS = [(15.75, 15.5, .235), (15.5, 15.5, .235), (16, 16, .295),
                      (15.5, 15.5, .235), (15.5, 15.5, .475), (16, 16, .235)]


def measured(depth, count=6):
    return [CandidateEvidence("L", count, str(index), measured_repeat_count=count,
                              raw_repeat_count=count, repeat_sigma=.08,
                              direct_product_measurement=True, full_repeat_span=True)
            for index in range(depth)]


def infer(evidence, locus=LOCUS, technology="ont"):
    calls, _ = infer_alleles(evidence, [locus], [], "sample", technology, LIMITS,
                            require_measurement=True, sample_mode="isolate")
    return calls


def finish(evidence, *, locus=LOCUS, technology="ont", updates=None):
    calls = infer(evidence, locus, technology)
    calls[0].update(updates or {})
    finalize_repeat_count_calls(calls, evidence, [locus], [], LIMITS)
    return calls[0]


@pytest.mark.parametrize("technology", ["illumina", "ont"])
@pytest.mark.parametrize("depth", [0, 2, 3, 200])
def test_primary_calls_require_independent_support_without_erasing_predictions(technology, depth):
    call = finish(measured(depth), technology=technology)
    assert call["call_accepted"] is (depth >= LIMITS.minimum_molecules)
    assert call["repeat_count"] == (6 if depth >= LIMITS.minimum_molecules else "")
    assert call["predicted_repeat_count"] == (6 if depth else 4)
    assert call["predicted_repeat_count_raw"] == (6 if depth else 4)
    assert call["presence_support"] == call["independent_support"] == depth
    assert bool(call["call_rejection_reason"]) is (depth < LIMITS.minimum_molecules)


@pytest.mark.parametrize("count", [0, 4.5, 400])
def test_supported_zero_fractional_and_novel_counts_remain_accepted(count):
    call = finish(measured(3, count))
    assert call["call_accepted"]
    assert call["repeat_count"] == call["predicted_repeat_count"] == count


def test_missing_data_without_an_eligible_prior_is_a_valid_unresolved_call():
    call = finish([], locus=Locus("L", repeat_unit_length_bp=4))
    assert not call["call_accepted"] and call["call_rejection_reason"]
    assert call["repeat_count"] == call["predicted_repeat_count"] == ""
    assert call["presence_support"] == call["independent_support"] == 0


def test_deep_partial_evidence_does_not_make_an_unbounded_count_identifiable():
    evidence = [CandidateEvidence("L", None, str(index),
                                  metadata={"repeat_lower_bound": 15})
                for index in range(300)]
    call = finish(evidence)
    assert call["independent_support"] == 300
    assert call["upper_unbounded"] and not call["finite_length_support"]
    assert not call["call_accepted"] and call["call_rejection_reason"]
    assert call["repeat_count"] == "" and call["predicted_repeat_count"] >= 15


def test_candidate_copies_cannot_turn_two_molecules_into_an_accepted_call():
    evidence = measured(2)
    duplicated = [replace(row, candidate_id=f"copy-{index}")
                  for row in evidence for index in range(100)]
    call = finish(duplicated)
    assert call == finish(evidence)
    assert call["independent_support"] == 2 and not call["call_accepted"]


@pytest.mark.parametrize("evidence,locus", [
    (measured(2), LOCUS),
    (measured(3), LOCUS),
    ([], LOCUS),
    ([], Locus("L")),
    ([CandidateEvidence("L", None, "m", metadata={"repeat_lower_bound": 15})], LOCUS),
])
def test_finalizing_again_preserves_predictions_and_acceptance(evidence, locus):
    call = finish(evidence, locus=locus)
    original = deepcopy(call)
    finalize_repeat_count_calls([call], evidence, [locus], [], LIMITS)
    assert call == original


@pytest.mark.parametrize("updates", [
    {"consensus_status": "disagrees"},
    {"consensus_status": "resolved_ambiguity"},
    {"computational_fallback": "candidate_limit"},
    {"possible_secondary": True},
    {"upper_unbounded": True, "finite_length_support": False},
    {"status": "estimated", "finite_length_support": True},
    {"status": "mixed"},
    {"best_probability": .5},
    {"margin": .1},
])
def test_depth_and_a_numeric_prediction_cannot_override_unresolved_evidence(updates):
    call = finish(measured(100), updates=updates)
    assert not call["call_accepted"] and call["call_rejection_reason"]
    assert call["repeat_count"] == "" and call["predicted_repeat_count"] == 6


def test_supported_reconstruction_is_accepted_without_a_likelihood_margin():
    calls = infer([])
    audit = [dict(locus_id="L", status="RECONSTRUCTED", method="mapped_flank_consensus",
                  repeat_count=6, raw_repeat_count=6, reconstruction_molecules=4,
                  reconstruction_fraction=.8)]
    apply_reconstructed_calls(calls, audit, LIMITS)
    assert calls[0]["probability_kind"] == "reconstruction_support_fraction"
    assert calls[0]["margin"] == 0
    finalize_repeat_count_calls(calls, [], [LOCUS], [], LIMITS)
    assert calls[0]["call_accepted"] and calls[0]["repeat_count"] == 6
    assert calls[0]["independent_support"] == 4


@pytest.mark.parametrize("method,support,fraction,allele,accepted", [
    ("mapped_flank_consensus", 4, 1., 6, True),
    ("mapped_flank_consensus", 2, 1., 6, False),
    ("mapped_flank_consensus", 4, .5, 6, False),
    ("mapped_flank_consensus", 4, 1., 8, False),
    ("complete_molecules", 4, 1., 6, False),
    ("mapped_local_assembly", 4, 1., 6, False),
])
def test_supported_consensus_can_resolve_low_depth_measurements(method, support, fraction, allele, accepted):
    evidence = measured(2) + [CandidateEvidence("L", None, f"unmeasured-{i}") for i in range(2)]
    calls = infer(evidence)
    assert calls[0]["status"] == "low_coverage"
    apply_reconstructed_calls(calls, [dict(
        locus_id="L", status="RECONSTRUCTED", method=method, repeat_count=allele,
        raw_repeat_count=allele, reconstruction_molecules=support,
        reconstruction_fraction=fraction)], LIMITS)
    finalize_repeat_count_calls(calls, evidence, [LOCUS], [], LIMITS)
    call = calls[0]
    assert call["call_accepted"] is accepted
    assert call["repeat_count"] == (6 if accepted else "")
    assert call["presence_support"] == 4
    if accepted:
        assert call["independent_support"] == call["molecule_support"] == 4
        assert call["primary_molecule_support"] == 2


@pytest.mark.parametrize("copies", [1, 10])
@pytest.mark.parametrize("reconstructed", [False, True])
def test_high_model_support_needs_molecule_agreement_or_a_supported_consensus(copies, reconstructed):
    # Numeric measurements from the retained noisy synthetic challenge: the
    # old model strongly selected 15.5 although only four of six molecules
    # rounded to it. Increasing coverage preserves that disagreement.
    evidence = [CandidateEvidence("L", allele, str(index), measured_repeat_count=allele,
                                  raw_repeat_count=raw, repeat_sigma=sigma, full_repeat_span=True)
                for index, (raw, allele, sigma) in enumerate(NOISY_MEASUREMENTS * copies)]
    calls = infer(evidence, technology="illumina")
    assert calls[0]["status"] == "called" and calls[0]["best_probability"] > .97
    assert calls[0]["repeat_count"] == 15.5
    if reconstructed:
        apply_reconstructed_calls(calls, [dict(
            locus_id="L", status="RECONSTRUCTED", method="mapped_flank_consensus",
            repeat_count=15.5, raw_repeat_count=15.75, reconstruction_molecules=6 * copies,
            reconstruction_fraction=1.0)], LIMITS)
        assert calls[0]["consensus_status"] == "agrees"
    finalize_repeat_count_calls(calls, evidence, [LOCUS], [], LIMITS)
    call = calls[0]
    assert call["predicted_repeat_count"] == 15.5
    assert call["primary_molecule_support"] == 4 * copies
    assert call["molecule_support"] == 6 * copies
    assert call["call_accepted"] is reconstructed
    assert call["repeat_count"] == (15.5 if reconstructed else "")


@pytest.mark.parametrize("reconstructed_allele", [15.5, 16])
def test_incidental_microassembly_support_cannot_bypass_a_conflicting_allele(reconstructed_allele):
    evidence = [CandidateEvidence("L", allele, str(index), measured_repeat_count=allele,
                                  raw_repeat_count=raw, repeat_sigma=sigma, full_repeat_span=True)
                for index, (raw, allele, sigma) in enumerate(NOISY_MEASUREMENTS * 5)]
    # A graph can use already measured molecules too. Only its one additional
    # supporting molecule appears as derived evidence in the shared caller.
    evidence.append(CandidateEvidence("L", reconstructed_allele, "derived",
                    measured_repeat_count=reconstructed_allele, raw_repeat_count=reconstructed_allele,
                    repeat_sigma=.08, evidence_tier="microassembly"))
    calls = infer(evidence, technology="illumina")
    assert calls[0]["status"] == "called" and calls[0]["best_probability"] > .999
    assert calls[0]["repeat_count"] == 15.5
    assert calls[0]["molecule_allele_counts"] == {15.5: 20, 16: 10}
    assert calls[0]["derived_molecule_support"] == 1
    finalize_repeat_count_calls(calls, evidence, [LOCUS], [], LIMITS)
    call = calls[0]
    assert call["predicted_repeat_count"] == 15.5
    assert call["call_accepted"] is (reconstructed_allele == 15.5)
    assert call["repeat_count"] == (15.5 if reconstructed_allele == 15.5 else "")


def test_supported_microassembly_without_original_measurements_remains_callable():
    evidence = [replace(row, evidence_tier="microassembly", direct_product_measurement=False,
                        full_repeat_span=False) for row in measured(3)]
    call = finish(evidence, technology="illumina")
    assert call["molecule_allele_counts"] == {} and call["derived_molecule_support"] == 3
    assert call["call_accepted"] and call["repeat_count"] == 6


@pytest.mark.parametrize("method", ["mapped_local_assembly", "complete_molecules"])
@pytest.mark.parametrize("reconstructed_allele", [6, 8])
def test_reconstructed_products_check_agreement_without_derived_only_molecules(method, reconstructed_allele):
    evidence = measured(3)
    calls = infer(evidence)
    assert calls[0]["derived_molecule_support"] == 0
    # Graph products retain their prior support validation, but their audit
    # has zero new molecules when all supporting molecules were already sized.
    audit = [dict(locus_id="L", status="RECONSTRUCTED", method=method,
                  repeat_count=reconstructed_allele, raw_repeat_count=reconstructed_allele,
                  reconstruction_molecules=0 if method == "mapped_local_assembly" else 3,
                  reconstruction_fraction=0.0 if method == "mapped_local_assembly" else 1.0)]
    apply_reconstructed_calls(calls, audit, LIMITS)
    assert calls[0]["consensus_status"] == ("agrees" if reconstructed_allele == 6 else "disagrees")
    finalize_repeat_count_calls(calls, evidence, [LOCUS], [], LIMITS)
    call = calls[0]
    assert call["predicted_repeat_count"] == 6
    assert call["call_accepted"] is (reconstructed_allele == 6)
    assert call["repeat_count"] == (6 if reconstructed_allele == 6 else "")


def test_rejected_predictions_stay_out_of_primary_exports_and_distance_inputs(tmp_path):
    evidence = measured(2)
    call = finish(evidence)
    write_evidence_reports([call], evidence, [LOCUS], tmp_path, LIMITS)
    with (tmp_path / "genotype_evidence.tsv").open() as handle:
        report = next(csv.DictReader(handle, delimiter="\t"))
    assert report["repeat_count"] == "" and float(report["predicted_repeat_count"]) == 6
    assert report["call_accepted"] == "False" and report["call_rejection_reason"]
    compact = common_calls_to_compatibility([call])[0]
    allele = dict(locus_id="L", called_repeat_count=compact["repeat_count"],
                  call_status=compact["status"], read_depth=compact["read_depth"],
                  dominant_vntr_asv=6,
                  **{field: compact[field] for field in CONFIDENCE_FIELDS})
    simple = simple_call_rows_from_alleles("sample", [allele])[0]
    fingerprint, probabilistic = build_fingerprint("sample", [allele], [LOCUS])
    matrix, _ = build_call_matrix([SimpleNamespace(sample_id="sample", rows=[compact])], ["L"])
    assert compact["repeat_count"] == simple["repeat_count"] == fingerprint[0]["L"] == ""
    assert probabilistic[0]["repeat_count"] == "" and np.isnan(matrix[0, 0])
    for row in (compact, simple, probabilistic[0]):
        assert row["predicted_repeat_count"] == 6
        assert not row["call_accepted"] and row["call_rejection_reason"]
    assert compact["read_depth"] == simple["read_depth"] == 2
