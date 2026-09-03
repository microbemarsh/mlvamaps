from __future__ import annotations

from mlvamaps.alignment_evidence import CandidateEvidence
from mlvamaps.allele_inference import InferenceThresholds, infer_alleles
from mlvamaps.candidate_contexts import CandidateContext, generate_candidate_contexts
from mlvamaps.models import Locus


def _locus() -> Locus:
    return Locus(
        "L1", forward_primer="AAAA", reverse_primer="CCCC",
        left_flank_sequence="GGGG", right_flank_sequence="TTTT",
        repeat_motif="AT", repeat_unit_length_bp=2,
        expected_min_repeats=2, expected_max_repeats=5, nominal_repeat_units=3,
    )


def test_candidate_generation_is_bounded_and_structured():
    contexts = generate_candidate_contexts([_locus()], maximum=4)
    assert [context.repeat_count for context in contexts] == [2, 3, 4]
    assert all(context.repeat_end - context.repeat_start == context.repeat_count * 2 for context in contexts)
    assert {context.row()["locus_id"] for context in contexts} == {"L1"}


def test_shared_inference_low_coverage_direct_measurement_is_retained():
    locus = _locus()
    contexts = generate_candidate_contexts([locus])
    evidence = [
        CandidateEvidence(
            "L1", context.repeat_count, "m1", alignment_score=100,
            alignment_identity=1.0, direct_product_measurement=True,
            full_repeat_span=True, technology="hifi", measured_repeat_count=4,
            candidate_id=context.candidate_id,
        )
        for context in contexts
    ]
    calls, molecule_calls = infer_alleles(
        evidence, [locus], contexts, "sample", "hifi",
        InferenceThresholds(minimum_molecules=3),
    )
    assert calls[0]["repeat_count"] == 4
    assert calls[0]["status"] == "low_coverage"
    assert calls[0]["direct_product_support"] == 1
    assert molecule_calls[("L1", "m1")] == 4


def test_generic_flank_mapping_detects_but_does_not_call():
    locus = _locus()
    contexts = generate_candidate_contexts([locus])
    context = contexts[0]
    evidence = [CandidateEvidence(
        "L1", context.repeat_count, "m1", alignment_score=40,
        alignment_identity=1.0, technology="illumina",
        candidate_id=context.candidate_id,
    )]
    calls, _ = infer_alleles(evidence, [locus], contexts, "sample", "illumina")
    assert calls[0]["repeat_count"] == ""
    assert calls[0]["status"] == "detected_unresolved"


def test_duplicate_reference_contexts_do_not_multiply_state_support():
    locus = _locus()
    contexts = generate_candidate_contexts([locus])
    by_repeat = {context.repeat_count: context for context in contexts}
    evidence = [
        CandidateEvidence(
            "L1", 3, "m1", alignment_score=100, alignment_identity=1.0,
            full_repeat_span=True, technology="hifi", candidate_id=f"duplicate{i}",
        ) for i in range(20)
    ] + [CandidateEvidence(
        "L1", 4, "m1", alignment_score=104, alignment_identity=1.0,
        full_repeat_span=True, technology="hifi", candidate_id=by_repeat[4].candidate_id,
    )]
    calls, molecule_calls = infer_alleles(
        evidence, [locus], contexts, "sample", "hifi",
        InferenceThresholds(minimum_molecules=1, minimum_probability=0.5, minimum_margin=0),
    )
    assert molecule_calls[("L1", "m1")] == 4
    assert calls[0]["best_candidate_repeat"] == 4


def test_independent_molecule_mixture_requires_support():
    locus = _locus()
    contexts = generate_candidate_contexts([locus])
    evidence = []
    for repeat, count in ((3, 6), (5, 4)):
        for index in range(count):
            evidence.append(CandidateEvidence(
                "L1", repeat, f"{repeat}-{index}", alignment_score=100,
                alignment_identity=1.0, direct_product_measurement=True,
                full_repeat_span=True, technology="hifi", measured_repeat_count=repeat,
            ))
    calls, _ = infer_alleles(evidence, [locus], contexts, "sample", "hifi")
    assert calls[0]["status"] == "mixed"
    assert calls[0]["dominant_repeat"] == 3
    assert calls[0]["secondary_repeat"] == 5
    assert calls[0]["dominant_fraction"] == 0.6