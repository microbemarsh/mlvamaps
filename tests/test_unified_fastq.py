from __future__ import annotations

from pathlib import Path

from mlvamaps.alignment_evidence import CandidateEvidence
from mlvamaps.allele_inference import InferenceThresholds, infer_alleles
from mlvamaps.candidate_contexts import (
    CandidateContext,
    candidate_repeat_counts,
    generate_candidate_contexts,
    write_candidate_contexts,
)
from mlvamaps.models import Locus
from mlvamaps.minimap_mapping import build_competitive_indexes
from mlvamaps.unified_fastq import taxonomic_query_sequences


def _locus() -> Locus:
    return Locus(
        "L1", forward_primer="AAAA", reverse_primer="CCCC",
        left_flank_sequence="GGGG", right_flank_sequence="TTTT",
        repeat_motif="AT", repeat_unit_length_bp=2,
        expected_min_repeats=2, expected_max_repeats=5, nominal_repeat_units=3,
    )


def test_taxonomic_query_uses_the_background_supported_by_most_molecules():
    contexts = [
        CandidateContext("c1", "L1", "AAAATATACCCC", 3, 2, 4, 10, "AAAA", "CCCC"),
        CandidateContext("c2", "L1", "GGGGTATATTTT", 3, 2, 4, 10, "GGGG", "TTTT"),
    ]
    evidence = [
        CandidateEvidence("L1", 3, "m1", 40, candidate_id="c1"),
        CandidateEvidence("L1", 3, "m2", 42, candidate_id="c1"),
        CandidateEvidence("L1", 3, "m3", 80, candidate_id="c2"),
    ]

    assert taxonomic_query_sequences(
        [{"locus": "L1", "repeat_count": 3, "status": "called"}],
        evidence,
        {("L1", molecule): 3 for molecule in ("m1", "m2", "m3")},
        contexts,
    ) == {"L1": "AAAATATACCCC"}
    assert taxonomic_query_sequences(
        [{"locus": "L1", "repeat_count": 3, "status": "called"}],
        evidence[::2],
        {("L1", "m1"): 3, ("L1", "m3"): 3},
        contexts,
    ) == {}
    tied = CandidateEvidence(
        "L1", 3, "m1", 40, technology="illumina", candidate_id="c1",
        metadata={"background_alignment_margin": 0},
    )
    assert taxonomic_query_sequences(
        [{"locus": "L1", "repeat_count": 3, "status": "called"}],
        [tied], {("L1", "m1"): 3}, contexts,
    ) == {}


def test_candidate_generation_is_bounded_and_structured():
    contexts = generate_candidate_contexts([_locus()], maximum=4)
    assert [context.repeat_count for context in contexts] == [2, 3, 4]
    assert all(context.repeat_end - context.repeat_start == context.repeat_count * 2 for context in contexts)
    assert {context.row()["locus_id"] for context in contexts} == {"L1"}


def test_invalid_candidate_range_is_rejected():
    locus = Locus(
        "bad", repeat_motif="AT", repeat_unit_length_bp=2,
        expected_min_repeats=8, expected_max_repeats=4,
    )
    import pytest
    with pytest.raises(ValueError, match="Invalid candidate repeat range"):
        candidate_repeat_counts(locus, [6])


def test_candidates_preserve_backgrounds_and_deduplicate_provenance(tmp_path):
    locus = _locus()
    database = tmp_path / "database"
    database.mkdir()
    sequence = "AAAA" + "GGGG" + "AT" * 3 + "TTTT" + "GGGG"
    (database / "L1.fasta").write_text(
        f">R1\n{sequence}\n>R2\n{sequence}\n>R3\n{sequence[:-1]}A\n"
    )
    (database / "reference_metadata.tsv").write_text(
        "reference_id\ttaxon_id\ttaxon_name\n"
        "R1\t1\tone\nR2\t1\tone\nR3\t2\ttwo\n"
    )
    contexts = generate_candidate_contexts([locus], database, expansion=0)
    assert len(contexts) == 8
    assert {context.provenance_count for context in contexts} == {1, 2}
    assert len({context.background_id for context in contexts}) == 2
    paths = write_candidate_contexts(contexts, database / "competitive_mapping")
    assert paths["provenance_table"].is_file()
    assert {context.observed_or_synthetic for context in contexts} == {"observed", "synthetic"}


def test_calibrated_legacy_panel_uses_observed_repeat_sequence_as_template(tmp_path):
    """An N placeholder denotes an unknown motif, not an unusable VNTR."""
    locus = Locus(
        "legacy_4bp_32bp_4U",
        forward_primer="AAAA",
        reverse_primer="CCCC",
        repeat_motif="NNNN",
        repeat_unit_length_bp=4,
        expected_product_size_bp=32,
        nominal_repeat_units=4,
        expected_min_repeats=3,
        expected_max_repeats=5,
    )
    database = tmp_path / "database"
    database.mkdir()
    observed = "AAAA" + "ATGC" * 4 + "GGGG" + "GGGG"
    (database / f"{locus.locus_id}.fasta").write_text(f">R1\n{observed}\n")
    (database / "reference_metadata.tsv").write_text(
        "reference_id\ttaxon_id\ttaxon_name\nR1\t1\tone\n"
    )

    contexts = generate_candidate_contexts([locus], database, expansion=0)

    assert [context.repeat_count for context in contexts] == [3, 4, 5]
    assert next(context for context in contexts if context.repeat_count == 4).sequence == observed
    assert all("N" not in context.sequence for context in contexts)


def test_unusable_candidate_error_names_locus_and_missing_state(tmp_path):
    locus = Locus("unconfigured", forward_primer="AAAA", reverse_primer="CCCC")
    database = tmp_path / "database"
    database.mkdir()
    (database / "unconfigured.fasta").write_text(">R1\nAAAATTTTGGGG\n")

    import pytest
    with pytest.raises(ValueError, match="unconfigured: repeat-unit length is missing"):
        generate_candidate_contexts([locus], database)


def test_short_and_long_indexes_use_distinct_parameters(tmp_path, monkeypatch):
    commands = []
    executable = tmp_path / "minimap2"
    executable.write_text("test")
    monkeypatch.setattr("mlvamaps.minimap_mapping.shutil.which", lambda _name: str(executable))

    def fake_run(command, **kwargs):
        commands.append(command)
        Path(command[command.index("-d") + 1]).write_text("index")
        import subprocess
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("mlvamaps.minimap_mapping.subprocess.run", fake_run)
    fasta = tmp_path / "candidate_contexts.fasta"
    fasta.write_text(">candidate0000001\nAAAA\n")
    indexes = build_competitive_indexes(fasta, tmp_path)
    assert indexes["short"].name == "short.mmi"
    assert indexes["long"].name == "long.mmi"
    assert commands[0][1:5] == ["-k", "21", "-w", "11"]
    assert commands[1][1:5] == ["-k", "11", "-w", "5"]


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
