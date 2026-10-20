"""Raw-length inference, mixture identifiability, and derived-evidence guards."""

from dataclasses import replace
import math

import pytest

from mlvamaps.alignment_evidence import CandidateAlignment, CandidateEvidence
from mlvamaps.allele_inference import InferenceThresholds, infer_alleles
from mlvamaps.calling import legacy_round_repeat_count
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.models import Locus
from mlvamaps.short_read_evidence import extract_short_read_evidence


LOCUS = Locus("L", repeat_motif="ATGC", repeat_unit_length_bp=4,
              expected_product_size_bp=220, nominal_repeat_units=5,
              expected_min_repeats=2, expected_max_repeats=8)


def rows(values, sigma=.2):
    return [CandidateEvidence(
        "L", 4, f"m{index}", direct_product_measurement=True,
        measured_repeat_count=legacy_round_repeat_count(value),
        raw_repeat_count=value, repeat_sigma=sigma,
    ) for index, value in enumerate(values)]


def call(evidence, **kwargs):
    return infer_alleles(evidence, [LOCUS], [], "sample", "hifi",
                        require_measurement=True, **kwargs)


def test_raw_lengths_are_combined_before_rounding_and_references_add_no_votes():
    evidence = rows([4.75, 4.75, 5, 5, 5.25, 5.25] * 2, sigma=.3)
    result, molecules = call(evidence, sample_mode="isolate")
    duplicated, other_molecules = call(evidence * 30, sample_mode="isolate")
    assert result == duplicated and molecules == other_molecules
    assert result[0]["repeat_count_raw"] == result[0]["repeat_count"] == 5
    assert result[0]["status"] == "called"
    assert result[0]["molecule_support"] == 12
    assert result[0]["em_iterations"] == 0


@pytest.mark.parametrize("value", [0, 5.25, 5.75, 15.5])
@pytest.mark.parametrize("tolerance", [0, .25, .4])
def test_novel_zero_and_fractional_counts_keep_assembly_rounding(value, tolerance):
    result, _ = call(rows([value] * 5, .08), sample_mode="isolate",
                     round_tolerance=tolerance)
    assert result[0]["repeat_count_raw"] == value
    assert result[0]["repeat_count"] == legacy_round_repeat_count(value, tolerance)


def test_em_recovers_allele_fractions_without_returning_the_middle_of_a_mixture():
    evidence = rows([3] * 12 + [7] * 8, sigma=.12)
    result, molecules = call(evidence, sample_mode="metagenome")
    row = result[0]
    assert row["status"] == "mixed" and row["repeat_count"] == 3
    assert row["dominant_fraction"] == pytest.approx(.6, abs=.01)
    assert row["secondary_repeat"] == 7
    assert row["secondary_fraction"] == pytest.approx(.4, abs=.01)
    assert row["em_converged"] and 0 < row["em_iterations"] <= 200
    assert row["probability_kind"] == "single_allele_model_support"
    assert set(molecules.values()) == {3, 7}


def test_weak_lengths_and_flank_mappings_do_not_manufacture_a_mixture():
    evidence = rows([3, 4] * 4, sigma=1)
    baseline, _ = call(evidence)
    evidence += [CandidateEvidence("L", 20, f"flank{i}", alignment_score=1000)
                 for i in range(100)]
    result, molecules = call(evidence)
    assert result[0]["status"] != "mixed"
    assert result[0]["candidate_distribution"] == baseline[0]["candidate_distribution"]
    assert result[0]["molecule_support"] == 8 and len(molecules) == 8


def test_em_can_use_many_individually_ambiguous_measurements():
    result, _ = call(rows([3] * 60 + [3.5] * 40, sigma=.2))
    assert result[0]["status"] == "mixed"
    assert result[0]["dominant_repeat"] == 3
    assert result[0]["secondary_repeat"] == 3.5
    assert result[0]["dominant_fraction"] == pytest.approx(.6, abs=.03)


def test_unconverged_mixture_is_not_promoted_to_a_confident_single_allele(monkeypatch):
    from mlvamaps.mapping_classification import fit_reference_mixture

    monkeypatch.setattr(
        "mlvamaps.mapping_classification.fit_reference_mixture",
        lambda likelihoods, weights, **kwargs: fit_reference_mixture(
            likelihoods, weights, max_iterations=1, tolerance=1e-6),
    )
    result, _ = call(rows([3] * 12 + [7] * 8, sigma=.12))
    assert not result[0]["em_converged"]
    assert result[0]["status"] == "ambiguous" and result[0]["repeat_count"] == ""


def test_isolate_model_tolerates_an_indel_outlier():
    result, _ = call(rows([5] * 12 + [9]), sample_mode="isolate")
    assert result[0]["repeat_count"] == 5 and result[0]["status"] == "called"


def test_noisier_product_qualities_increase_length_uncertainty():
    from mlvamaps.locus_measurement import measure_locus_product
    from mlvamaps.sequence import revcomp

    locus = replace(LOCUS, forward_primer="AACCGGTTAACCGG",
                    reverse_primer="GGTTCCAAGGTTCC")
    sequence = locus.forward_primer + "ATGC" * 40 + revcomp(locus.reverse_primer)
    high = measure_locus_product(sequence, locus, "I" * len(sequence))
    low = measure_locus_product(sequence, locus, "0" * len(sequence))
    assert high.raw_repeat_count == low.raw_repeat_count
    assert low.metadata["repeat_sigma"] > high.metadata["repeat_sigma"]


def test_consensus_can_resolve_ambiguity_but_never_adds_molecule_depth():
    evidence = rows([5, 5.5, 5, 5.5], sigma=.5)
    thresholds = InferenceThresholds(minimum_probability=1)
    before, molecules = call(evidence, sample_mode="isolate", thresholds=thresholds)
    assert before[0]["status"] == "ambiguous"
    consensus = {"L": {"raw_repeat_count": 5.25, "supporting_reads": 4}}
    result, after_molecules = call(evidence, sample_mode="isolate", thresholds=thresholds,
                                  consensus_measurements=consensus)
    assert result[0]["estimation_method"] == "consensus_supported"
    assert result[0]["repeat_count"] == 5.5 and result[0]["confidence"] == "low"
    assert result[0]["molecule_support"] == 4 and molecules == after_molecules
    consensus["L"]["raw_repeat_count"] = 12
    rejected, _ = call(evidence, sample_mode="isolate", thresholds=thresholds,
                       consensus_measurements=consensus)
    assert rejected[0]["repeat_count"] == ""
    assert rejected[0]["consensus_status"] == "disagrees"
    empty, _ = call([], consensus_measurements=consensus)
    assert empty[0]["repeat_count"] == ""


@pytest.mark.parametrize("sigma", [0, -1, math.inf, math.nan])
def test_invalid_measurement_uncertainty_is_rejected(sigma):
    with pytest.raises(ValueError, match="uncertainty"):
        call(rows([3], sigma=sigma))


def pair_alignments():
    context = CandidateContext("ref", "L", "A" * 220, 5, 4, 100, 120, "", "")
    first = CandidateAlignment(
        molecule_id="m", read_id="m/1", mate=1, locus_id="L", candidate_id="ref",
        repeat_count=5, reference_id="R", alignment_score=90, mapping_quality=60,
        alignment_identity=1, query_coverage=1, reference_coverage=.2,
        query_start=0, query_end=40, reference_start=30, reference_end=70,
        cigar="40M", cs=":40", primary=True, secondary=False, reverse=False,
        query_sequence="A" * 40, query_quality="I" * 40,
    )
    second = replace(first, read_id="m/2", mate=2, reference_start=150,
                     reference_end=190, reverse=True)
    return context, first, second


def test_pair_length_estimation_requires_calibration_and_remains_an_estimate():
    context, first, second = pair_alignments()
    alignments = [replace(row, molecule_id=f"m{i}") for i in range(5) for row in (first, second)]
    evidence = extract_short_read_evidence(alignments, [context], [LOCUS])
    assert call(evidence)[0][0]["repeat_count"] == ""
    # Reference outer span is 160 bp. A 168 bp library mean adds two units.
    evidence = extract_short_read_evidence(alignments, [context], [LOCUS],
                                          insert_mean=168, insert_sd=8)
    result, molecules = call(evidence, sample_mode="isolate")
    assert result[0]["repeat_count"] == 7
    assert result[0]["status"] == "estimated" and result[0]["confidence"] == "low"
    assert result[0]["estimation_method"] == "paired_fragment_length"
    assert result[0]["estimate_support"] == 5 and result[0]["molecule_support"] == 0
    assert molecules == {}
    duplicated, _ = call(evidence * 10, sample_mode="isolate")
    assert result == duplicated
    direct = rows([5] * 3)
    assert call(evidence + direct, sample_mode="isolate")[0][0]["repeat_count"] == 5
    outward = extract_short_read_evidence([replace(first, reverse=True), second], [context], [LOCUS],
                                         insert_mean=168, insert_sd=8)
    assert call(outward)[0][0]["repeat_count"] == ""
    with pytest.raises(ValueError, match="Insert-size"):
        extract_short_read_evidence([], [], [], insert_mean=168)


def test_long_read_measurement_cost_does_not_scale_with_reference_copies(monkeypatch):
    from types import SimpleNamespace
    from mlvamaps.long_read_evidence import extract_long_read_evidence

    measured = []

    def measure(sequence, locus, quality, **kwargs):
        measured.append((sequence, quality))
        return SimpleNamespace(called_allele=5, status="FULL_PRODUCT", confidence=.99,
                               raw_repeat_count=5, metadata={"repeat_sigma": .1})

    monkeypatch.setattr("mlvamaps.long_read_evidence.measure_locus_product", measure)
    context, prototype, _ = pair_alignments()
    contexts = [replace(context, candidate_id=f"ref{i}", repeat_count=i) for i in range(10)]
    alignments = [replace(prototype, molecule_id=f"m{i}", mate=None, candidate_id=c.candidate_id,
                          repeat_count=c.repeat_count, query_quality=("I" if i % 2 else "H") * 40)
                  for c in contexts for i in range(20)]
    evidence = extract_long_read_evidence(alignments, contexts, [LOCUS], "hifi")
    assert len(measured) == 4  # Two qualities, two orientations; one locus.
    assert len(evidence) == 200
    result, molecules = call(evidence, sample_mode="isolate")
    assert result[0]["molecule_support"] == len(molecules) == 20


def test_cli_requires_both_library_parameters(capsys):
    from mlvamaps.cli import main

    with pytest.raises(SystemExit) as error:
        main(["call", "-i", "sr", "--short-insert-mean", "300"])
    assert error.value.code == 2
    assert "must be supplied together" in capsys.readouterr().err
