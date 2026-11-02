"""Complete numeric predictions retain independent evidence semantics."""
from mlvamaps.pipeline import simple_call_rows_from_alleles
from mlvamaps.profile_matching import build_fingerprint
from mlvamaps.models import Locus
from dataclasses import replace
import pytest
from mlvamaps.alignment_evidence import CandidateEvidence
from mlvamaps.allele_inference import (InferenceThresholds, infer_alleles, estimate_unresolved_calls,
                                      finalize_repeat_count_calls, locus_repeat_prior)
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.evidence_reporting import write_evidence_reports
from mlvamaps.profile_matching import match_profiles, profile_match_locus_rows


def complete(evidence=(), locus=None, contexts=(), limits=None, mode="isolate"):
    locus = locus or Locus("L", repeat_unit_length_bp=4, nominal_repeat_units=4, expected_product_size_bp=200)
    limits = limits or InferenceThresholds()
    calls, molecules = infer_alleles(list(evidence), [locus], list(contexts), "s", "ont", limits,
                                   require_measurement=True, sample_mode=mode)
    estimate_unresolved_calls(calls, list(evidence), list(contexts))
    finalize_repeat_count_calls(calls, list(evidence), [locus], list(contexts), limits)
    return calls[0], molecules


def test_zero_depth_exports_keep_selected_count_and_absence():
    alleles = [dict(locus_id="L", called_repeat_count=4, read_depth=0,
                    call_status="LOCUS_DROPOUT")]
    row = simple_call_rows_from_alleles("s", alleles)[0]
    fingerprint, probabilistic = build_fingerprint("s", alleles, [Locus("L")])
    assert row["repeat_count"] == fingerprint[0]["L"] == probabilistic[0]["repeat_count"] == 4
    assert row["present"] == "no" and row["read_depth"] == 0


def test_prior_is_explicit_and_reference_copy_invariant():
    locus = Locus("L")
    context = CandidateContext("c", "L", "ACGT", 4, 4, 0, 4, "", "", observed_or_synthetic="observed")
    contexts = [context, replace(context, repeat_count=8)]
    assert locus_repeat_prior(locus, contexts)[0] == {4: .5, 8: .5}
    assert locus_repeat_prior(locus, contexts + [context] * 100) == locus_repeat_prior(locus, contexts)
    row, molecules = complete(locus=locus, contexts=contexts)
    assert row["repeat_count"] == 4 and row["prior_only"] and not molecules
    assert row["model_support"] == .5 and row["data_evidence_confidence"] == row["independent_support"] == 0
    assert row["presence_support"] == 0 and row["status"] == "not_found"
    assert row["second_best_repeat_count"] == 8
    with pytest.raises(ValueError, match="eligible locus-specific prior"):
        complete(locus=locus)
    assert complete(locus=replace(locus, expected_min_repeats=2, expected_max_repeats=8))[0]["repeat_count"] == 2


@pytest.mark.parametrize("mode", ["isolate", "metagenome"])
@pytest.mark.parametrize("limit", [dict(maximum_candidates=2), dict(maximum_likelihood_cells=1)])
def test_bounded_likelihood_keeps_large_noisy_measurements(mode, limit):
    rows = [CandidateEvidence("L", 4, str(i), measured_repeat_count=400,
            raw_repeat_count=400, repeat_sigma=100, direct_product_measurement=True) for i in range(3)]
    call, molecules = complete(rows, limits=InferenceThresholds(**limit), mode=mode)
    assert call["repeat_count"] == 400 and not call["prior_only"] and len(molecules) == 3
    assert call["computational_fallback"] and call["data_evidence_confidence"] < .01
    assert call["repeat_interval_high"] >= 700
    assert complete(rows * 10, limits=InferenceThresholds(**limit), mode=mode)[0] == call


def test_ambiguous_locus_is_never_sample_support_for_a_prior(tmp_path):
    locus = Locus("L", expected_min_repeats=4, expected_max_repeats=4)
    evidence = [CandidateEvidence("L", 9, "m", measured_repeat_count=9,
                quality_flags=("ambiguous_locus",), metadata={"pair_repeat_count": 9})]
    call, molecules = complete(evidence, locus)
    write_evidence_reports([call], evidence, [locus], tmp_path, InferenceThresholds())
    assert call["repeat_count"] == 4 and call["prior_only"] and not molecules
    assert call["data_evidence_confidence"] == call["independent_support"] == 0
    assert call["genotype_status"] == "AMBIGUOUS_LOCUS" and call["presence_support"] == 1


def test_prior_counts_are_exported_but_not_compared():
    row, _ = complete()
    alleles = [dict(locus_id="L", called_repeat_count=row["repeat_count"], **row)]
    fingerprint, _ = build_fingerprint("s", alleles, [Locus("L")])
    profiles = [dict(profile_id="reference", L=4)]
    matches = match_profiles("s", fingerprint[0], profiles, alleles)
    assert matches[0]["compared_loci"] == 0
    assert profile_match_locus_rows("s", fingerprint[0], profiles, matches, alleles)[0]["match_status"] == "NOT_COMPARED"
    from mlvamaps.myoga_export import build_call_matrix
    from types import SimpleNamespace
    import numpy as np
    matrix, rows = build_call_matrix([SimpleNamespace(sample_id="s", rows=[dict(locus_id="L", **row)])], ["L"])
    assert np.isnan(matrix[0, 0]) and rows["s"]["L"]["repeat_count"] == 4
