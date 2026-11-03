"""Shared sequence measurements never become separate locus assignments."""

from dataclasses import replace
import csv
from pathlib import Path
import random
import runpy
import shutil

import pytest

from mlvamaps.alignment_evidence import CandidateAlignment
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.locus_groups import infer_locus_groups
from mlvamaps.models import Locus
from mlvamaps.sequence import revcomp


def inputs(alleles=(5, 5, 5)):
    rng = random.Random(81)
    forward, left, right, reverse = ["".join(rng.choices("ACGT", k=n)) for n in (20, 25, 25, 20)]
    first = Locus("A", forward_primer=forward, reverse_primer=reverse,
        left_flank_sequence=left, right_flank_sequence=right, repeat_motif="ATGC",
        repeat_unit_length_bp=4, expected_product_size_bp=106, nominal_repeat_units=4)
    loci = [first, replace(first, locus_id="B")]
    reference = forward + left + "ATGC" * 4 + right + revcomp(reverse)
    contexts = [CandidateContext(locus.locus_id, locus.locus_id, reference, 4, 4, 45, 61, left, right,
        repeat_motif="ATGC", observed_or_synthetic="observed") for locus in loci]
    alignments = []
    for index, allele in enumerate(alleles):
        product = forward + left + "ATGC" * allele + right + revcomp(reverse)
        for context in contexts:
            alignments.append(CandidateAlignment(str(index), str(index), None, context.locus_id,
                context.candidate_id, 4, "ref", 200, 0, 1, 1, 1,
                0, len(product), 0, len(reference), f"45M{4 * (allele - 4)}I61M", "",
                context.locus_id == "A", context.locus_id != "A", False,
                query_sequence=product, query_quality="I" * len(product)))
    ambiguous = {str(index): ("A", "B") for index in range(len(alleles))}
    return alignments, contexts, loci, ambiguous


def test_shared_allele_uses_original_molecules_without_mutating_locus_evidence():
    alignments, contexts, loci, ambiguous = inputs()
    original = tuple(alignments)
    row, = infer_locus_groups(alignments, contexts, loci, ambiguous, "s")
    assert row["member_loci"] == "A;B" and row["locus_assignment"] == "AMBIGUOUS_LOCUS"
    assert row["group_status"] == "SHARED_ALLELE" and row["repeat_count"] == 5
    assert row["independent_support"] == row["presence_support"] == row["primary_molecule_support"] == 3
    assert tuple(alignments) == original and {context.locus_id for context in contexts} == {"A", "B"}
    duplicate_contexts = contexts + [replace(context, candidate_id=context.candidate_id + "copy") for context in contexts]
    duplicates = alignments + [replace(row, candidate_id=row.candidate_id + "copy") for row in alignments]
    assert infer_locus_groups(duplicates * 3, duplicate_contexts, list(reversed(loci)), ambiguous, "s") == [row]
    assert infer_locus_groups(list(reversed(alignments)), list(reversed(contexts)), loci, ambiguous, "s") == [row]


def test_group_alleles_can_differ_between_members_without_becoming_a_single_average():
    data = inputs((5,) * 6 + (9,) * 4)
    row, = infer_locus_groups(*data, "s")
    assert row["group_status"] == "MULTIPLE_ALLELES"
    assert row["repeat_count"] == 5 and row["secondary_repeat"] == 9
    assert row["independent_support"] == 10
    assert row["primary_molecule_support"] == 6 and row["secondary_molecule_support"] == 4
    assert "do not establish strain fractions" in row["reason"]


@pytest.mark.parametrize("field,value", [
    ("nominal_repeat_units", 5), ("repeat_unit_length_bp", 8),
    ("expected_product_size_bp", 110), ("left_flank_sequence", "ACGT" * 7),
    ("forward_primer", "ACGT" * 5),
])
def test_identical_catalog_sequences_with_different_measurement_rules_do_not_group(field, value):
    alignments, contexts, loci, ambiguous = inputs()
    loci[1] = replace(loci[1], **{field: value})
    assert infer_locus_groups(alignments, contexts, loci, ambiguous, "s") == []


@pytest.mark.parametrize("change", ["boundary", "catalog", "synthetic"])
def test_groups_require_matching_observed_catalogs_and_boundaries(change):
    alignments, contexts, loci, ambiguous = inputs()
    if change == "boundary":
        contexts[1] = replace(contexts[1], repeat_start=46)
    elif change == "catalog":
        contexts.append(replace(contexts[1], candidate_id="extra", sequence="A" + contexts[1].sequence[1:]))
    else:
        contexts = [replace(context, observed_or_synthetic="synthetic") for context in contexts]
    assert infer_locus_groups(alignments, contexts, loci, ambiguous, "s") == []


def test_outside_competitors_and_absent_evidence_never_supply_a_group_guess():
    alignments, contexts, loci, ambiguous = inputs()
    ambiguous = {name: (*members, "outside") for name, members in ambiguous.items()}
    for evidence, assignments in ((alignments, ambiguous), (alignments, {}), ([], {})):
        row, = infer_locus_groups(evidence, contexts, loci, assignments, "s")
        assert row["repeat_count"] == "" and row["group_status"] == "UNRESOLVED"
        assert row["independent_support"] == row["presence_support"] == 0


def test_partial_group_sequence_is_a_bound_not_a_full_allele():
    alignments, contexts, loci, ambiguous = inputs()
    partial = loci[0].left_flank_sequence + "ATGC" * 12
    alignments = [replace(row, query_sequence=partial, query_quality="I" * len(partial),
        query_end=len(partial), reference_start=20, reference_end=61,
        query_coverage=41 / len(partial), cigar="41M32S") for row in alignments]
    row, = infer_locus_groups(alignments, contexts, loci, ambiguous, "s")
    assert row["repeat_count"] == "" and row["repeat_lower_bound"] == 12
    assert row["group_status"] == "LOWER_BOUND" and row["independent_support"] == 0
    assert row["presence_support"] == 3


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"),
                    reason="Native tools are needed for the full shared-locus workflow")
def test_native_challenge_exports_shared_allele_without_locus_or_classification_assignment(tmp_path, monkeypatch):
    from mlvamaps.cli import main
    from mlvamaps import mapping_classification

    fixture = runpy.run_path(str(Path(__file__).resolve().parents[1] / "scripts/fastq_challenge_fixture.py"))
    argv, inputs, _ = fixture["prepare_challenge"](tmp_path / "inputs", "paired", copies=2)
    output = Path(argv[argv.index("-o") + 1])
    measurements = []
    classify = mapping_classification.run_mapping_classification

    def record(**kwargs):
        measurements.append(dict(kwargs.get("measured_repeats") or {}))
        return classify(**kwargs)

    monkeypatch.setattr(mapping_classification, "run_mapping_classification", record)
    assert main(argv) == 0

    def table(name):
        with (output / name).open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    group, = table("locus_group_calls.tsv")
    assert group["member_loci"] == "ambiguous_a;ambiguous_b"
    assert group["repeat_count"] == "5" and group["independent_support"] == "40"
    assert group["group_status"] == "SHARED_ALLELE"
    per_locus = {row["locus_id"]: row for row in table("calls.tsv")}
    for name in ("ambiguous_a", "ambiguous_b"):
        assert per_locus[name]["repeat_count"] == "4" and per_locus[name]["genotype_status"] == "AMBIGUOUS_LOCUS"
    assert all(locus not in {"ambiguous_a", "ambiguous_b", group["group_id"]}
               for locus, molecule in measurements[0])
    report = (output / "report.html").read_text()
    assert "Shared Alleles Across Indistinguishable Loci" in report and "SHARED_ALLELE" in report

    for path in inputs:
        Path(path).write_text("")
    assert main(argv) == 0
    empty, = table("locus_group_calls.tsv")
    assert empty["repeat_count"] == "" and empty["independent_support"] == "0"
    assert empty["group_status"] == "UNRESOLVED" and empty["presence_support"] == "0"
    assert not measurements[-1]
    assert "SHARED_ALLELE" not in (output / "report.html").read_text()
