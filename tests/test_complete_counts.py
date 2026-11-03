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
    estimate_unresolved_calls(calls, list(evidence), list(contexts), thresholds=limits)
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
    with pytest.raises(ValueError, match="eligible locus-specific prior"):
        complete(locus=replace(locus, expected_product_size_bp=200, repeat_unit_length_bp=4))
    assert complete(locus=replace(locus, expected_min_repeats=2, expected_max_repeats=8))[0]["repeat_count"] == 2
    zero = complete(locus=replace(locus, expected_max_repeats=0))[0]
    assert zero["repeat_count"] == 0 and zero["prior_only"]
    bounded = complete(locus=replace(locus, expected_max_repeats=1000000), limits=InferenceThresholds(maximum_candidates=1))[0]
    assert bounded["repeat_interval_high"] == 1000000 and bounded["computational_fallback"]


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


def test_partial_bounds_condition_prior_without_counting_mates_or_copies(tmp_path):
    locus = Locus("L", expected_product_size_bp=156, nominal_repeat_units=4, repeat_unit_length_bp=4)
    rows = [CandidateEvidence("L", 4, "m", metadata={"repeat_lower_bound": 15}),
            CandidateEvidence("L", 4, "m", metadata={"repeat_lower_bound": 12})]
    call, molecules = complete(rows, locus)
    assert call["repeat_count"] == 15 and not molecules and not call["prior_only"]
    assert call["independent_support"] == 1 and call["data_evidence_confidence"] < .1
    assert call["repeat_interval_high"] == "" and call["upper_unbounded"]
    assert call["candidate_grid_high"] > 15 and call["candidate_grid_low"] == 15
    assert complete(rows * 20, locus)[0] == call
    write_evidence_reports([call], rows, [locus], tmp_path, InferenceThresholds())
    assert call["genotype_status"] == "LOWER_BOUND"
    assert complete(rows, replace(locus, expected_product_size_bp=0))[0]["repeat_count"] == 15


def test_conflicting_candidates_use_joint_likelihood_and_never_reference_copy_votes():
    rows = [CandidateEvidence("L", value, molecule, measured_repeat_count=value,
                raw_repeat_count=value, repeat_sigma=.2, direct_product_measurement=True)
            for molecule, values in (("a", (4, 8)), ("b", (6, 8))) for value in values]
    call, molecules = complete(rows)
    assert call["repeat_count"] == 8 and call["estimate_source"] == "conflicting_measurements"
    assert call["independent_support"] == 2 and not molecules
    assert complete(rows + [rows[0]] * 100)[0] == call
    limited, _ = complete(rows, limits=InferenceThresholds(maximum_candidates=1))
    assert limited["repeat_count"] == 8 and limited["computational_fallback"]


def test_spanning_alignment_supersedes_partial_only_uncertainty():
    locus = Locus("L", repeat_unit_length_bp=4, nominal_repeat_units=4, expected_product_size_bp=156)
    context = CandidateContext("c", "L", "A" * 70 + "ATGC" * 4 + "C" * 70, 4, 4, 70, 86, "", "")
    evidence = [CandidateEvidence("L", 4, "partial", metadata={"repeat_lower_bound": 3}),
                CandidateEvidence("L", 4, "span", candidate_id="c", alignment_score=100,
                    alignment_identity=1, query_coverage=1,
                    metadata={"aligned_flank_bases": 30, "spanning_repeat_indel_delta": 2})]
    calls, molecules = infer_alleles(evidence, [locus], [context], "s", "illumina",
                                   require_measurement=True, sample_mode="isolate")
    assert calls[0]["upper_unbounded"]
    estimate_unresolved_calls(calls, evidence, [context])
    finalize_repeat_count_calls(calls, evidence, [locus], [context])
    call, = calls
    assert call["repeat_count"] == 6 and call["repeat_interval_high"] == 6
    assert not call["upper_unbounded"] and call["finite_length_support"]
    assert call["length_model_status"] == "superseded_by_spanning_alignment" and not molecules


def test_competing_pair_lengths_precede_partial_or_prior_predictions():
    rows = [CandidateEvidence("L", None, molecule, metadata={"pair_repeat_count": value,
                "pair_repeat_sigma": .2, "pair_alignment_score": 100, "repeat_lower_bound": 2})
            for molecule, values in (("a", (4, 8)), ("b", (6, 8))) for value in values]
    call, molecules = complete(rows)
    assert call["repeat_count"] == 8 and call["estimate_source"] == "paired_fragment_length"
    assert call["independent_support"] == 2 and not molecules and call["molecule_support"] == 0
    assert complete(rows + [rows[0]] * 100)[0] == call


def test_partial_bound_uses_full_flank_for_product_calibration():
    from mlvamaps.short_read_evidence import partial_repeat_lower_bound
    locus = Locus("L", forward_primer="A" * 20, reverse_primer="T" * 20,
                  left_flank_sequence="A" * 29 + "GACTGACCTGATCGACGTTAG",
                  right_flank_sequence="G" * 50, repeat_motif="ATGC", repeat_unit_length_bp=4,
                  expected_product_size_bp=156, nominal_repeat_units=4)
    locus = replace(locus, expected_product_size_bp=20 + len(locus.left_flank_sequence) + 50 + 20 + 16)
    assert partial_repeat_lower_bound(locus.left_flank_sequence + "ATGC" * 12, locus) == 12


def test_direct_sizing_takes_priority_over_partial_and_reference_evidence():
    direct = [CandidateEvidence("L", 400, str(i), measured_repeat_count=400,
              raw_repeat_count=400, repeat_sigma=.08, direct_product_measurement=True) for i in range(3)]
    partial = CandidateEvidence("L", 4, "other", metadata={"repeat_lower_bound": 12})
    assert complete(direct + [partial])[0]["repeat_count"] == 400
    noisy = [replace(row, repeat_sigma=.8) for row in direct]
    assert complete(noisy)[0]["data_evidence_confidence"] < complete(direct)[0]["data_evidence_confidence"]
    assert complete(noisy)[0]["repeat_interval_high"] > complete(direct)[0]["repeat_interval_high"]


@pytest.mark.parametrize("layout", ["paired", "long"])
@pytest.mark.parametrize("mode", ["isolate", "metagenome"])
def test_challenge_completion_export_agreement_and_determinism(tmp_path, layout, mode):
    import csv
    import math
    import shutil
    import sys
    from pathlib import Path
    if not shutil.which("minimap2") or not shutil.which("sassy"):
        pytest.skip("native FASTQ tools required")
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    from fastq_challenge_fixture import prepare_challenge
    from mlvamaps.cli import main
    argv, _, truth = prepare_challenge(tmp_path / "inputs", layout)
    argv[argv.index("--sample-mode") + 1] = mode
    assert main(argv) == 0
    root = Path(argv[argv.index("-o") + 1])
    def read(name, delimiter="\t"):
        with (root / name).open() as handle:
            return list(csv.DictReader(handle, delimiter=delimiter))
    calls = {row["locus"]: row for row in read("common_locus_calls.tsv")}
    assert set(calls) == set(truth) and len(read("common_locus_calls.tsv")) == len(truth)
    for call in calls.values():
        assert math.isfinite(float(call["repeat_count"])) and float(call["repeat_count"]) >= 0
    for name in ("calls.tsv", "locus_repeat_counts.tsv", "mlva_fingerprint_probabilistic.tsv", "genotype_evidence.tsv", "myoga_loci.csv"):
        rows = read(name, "," if name.endswith("csv") else "\t")
        assert len(rows) == len(truth)
        for row in rows:
            locus = row.get("locus") or row["locus_id"]
            assert float(row["repeat_count"]) == float(calls[locus]["repeat_count"])
            assert row["prior_only"] == calls[locus]["prior_only"]
    if layout == "long":
        for row in read("allele_calls.tsv"):
            assert float(row["called_repeat_count"]) == float(calls[row["locus_id"]]["repeat_count"])
        assert float(calls["long_array"]["repeat_count"]) == 400
    fingerprint = read("mlva_fingerprint.tsv")[0]
    assert all(float(fingerprint[locus]) == float(call["repeat_count"]) for locus, call in calls.items())
    for locus in ("absent", "zero_coverage"):
        assert calls[locus]["prior_only"] == "True"
        assert float(calls[locus]["data_evidence_confidence"]) == float(calls[locus]["independent_support"]) == 0
        assert calls[locus]["presence_support"] == "0"
    snapshots = {name: (root / name).read_bytes() for name in
                 ("calls.tsv", "common_locus_calls.tsv", "genotype_evidence.tsv", "mlva_fingerprint.tsv")}
    argv[argv.index("-o") + 1] = str(tmp_path / "repeat")
    assert main(argv) == 0
    assert all((tmp_path / "repeat" / name).read_bytes() == data for name, data in snapshots.items())
