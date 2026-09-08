from dataclasses import replace
import csv
import json
from pathlib import Path

import numpy as np
import pytest

from mlvamaps.alignment_evidence import CandidateAlignment
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.mapping_classification import (
    UNKNOWN, alignment_statistics, molecule_log_likelihoods,
    fit_reference_mixture, classify_molecules, run_mapping_classification,
)
from mlvamaps.models import Locus
from mlvamaps.phylogeny import _parse_newick, _tip_patristic_distances


def context(name, count=4, locus="L1"):
    return CandidateContext(name, locus, "ACGTAC" + "GA" * count + "TGCATG",
        count, 2, 6, 6 + 2 * count, "ACGTAC", "TGCATG")


def alignment(c, *, molecule="m1", cs=None, cigar=None, mate=None, end=None, score=100):
    end = len(c.sequence) if end is None else end
    return CandidateAlignment(
        molecule, molecule, mate, c.locus_id, c.candidate_id, c.repeat_count, "",
        score, 0, 1, 1, 1, 0, end, 0, end,
        cigar or f"{end}=", cs if cs is not None else f":{end}", True, False, False,
    )


def read_tsv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_repeat_indels_count_once_and_partial_flanks_do_not_imply_length():
    a, b = context("a"), context("b", 6)
    rows = [alignment(a), alignment(b, cs=":6-gaga:14", cigar="6=4D14=")]
    stats = alignment_statistics(rows[1], b)
    assert stats["repeat_delta"] == -2
    assert stats["deletions"] == 0
    scores, _ = molecule_log_likelihoods(rows, [a, b], {"a": ["A"], "b": ["B"]})
    assert scores[("L1", "m1")]["A"] == 0
    assert scores[("L1", "m1")]["B"] == -2
    paired = [replace(row, mate=mate) for row in rows for mate in (1, 2)]
    paired_scores, _ = molecule_log_likelihoods(paired, [a, b], {"a": ["A"], "b": ["B"]})
    assert paired_scores[("L1", "m1")]["B"] == -2
    partial, _ = molecule_log_likelihoods(
        [alignment(a, cs=":6", end=6), alignment(b, cs=":6", end=6)],
        [a, b], {"a": ["A"], "b": ["B"]},
    )
    assert partial[("L1", "m1")]["A"] == partial[("L1", "m1")]["B"] == 0


def test_sequence_mismatches_and_nonrepeat_gaps_still_contribute():
    a = context("a")
    perfect = alignment(a)
    mismatch = alignment(a, molecule="snp", cs=":2*ac:17")
    flank_gap = alignment(a, molecule="gap", cs=":2-a:17")
    scores, _ = molecule_log_likelihoods([perfect, mismatch, flank_gap], [a], {"a": ["A"]})
    assert scores[("L1", "m1")]["A"] == 0
    assert scores[("L1", "snp")]["A"] < 0
    assert scores[("L1", "gap")]["A"] < 0
    assert alignment_statistics(flank_gap, a)["repeat_delta"] == 0


def test_duplicate_alignments_and_reference_copies_do_not_add_support():
    a, b = context("a"), context("b", 6)
    reads = [alignment(a), alignment(b, cs=":6-gaga:14")]
    scores, _ = molecule_log_likelihoods(reads, [a, b], {"a": ["A"], "b": ["B"]})
    duplicated, _ = molecule_log_likelihoods(reads * 4, [a, b], {"a": ["A", "A_copy"], "b": ["B"]})
    original, _ = classify_molecules(scores, ["A", "B"])
    copies, _ = classify_molecules(duplicated, ["A", "A_copy", "B"])
    assert copies[0]["references"] == ["A", "A_copy"]
    assert copies[0]["model_weight"] == pytest.approx(original[0]["model_weight"])


def test_em_resolves_mixture_and_has_monotone_objective():
    values = np.array([[0, -12]] * 7 + [[-12, 0]] * 3 + [[0, 0]] * 5)
    fractions, history = fit_reference_mixture(values, np.ones(len(values)))
    assert fractions == pytest.approx([0.7, 0.3], abs=1e-4)
    assert np.diff(history).min() >= -1e-8
    assert fractions.sum() == pytest.approx(1)
    assert fit_reference_mixture(values, np.ones(len(values)))[0] == pytest.approx(fractions)


def test_joint_reference_identity_cannot_switch_between_loci():
    scores = {
        ("L1", "m1"): {"A": 0, "B": -8, UNKNOWN: -20},
        ("L2", "m2"): {"A": -8, "B": 0, UNKNOWN: -20},
    }
    groups, _ = classify_molecules(scores, ["A", "B"])
    assert groups[0]["model_weight"] == pytest.approx(groups[1]["model_weight"])
    assert groups[0]["log_likelihood"] == -8
    mixture, diagnostics = classify_molecules(scores, ["A", "B"], sample_mode="metagenome")
    assert mixture[0]["locus_balanced_fraction"] == pytest.approx(0.5, abs=1e-4)
    assert diagnostics["converged"]


def test_unclassified_component_prevents_forcing_poor_matches():
    groups, _ = classify_molecules({("L1", "m"): {"A": -50, UNKNOWN: -8}}, ["A"])
    assert groups[0]["references"] == [UNKNOWN]
    assert groups[0]["model_weight"] > 0.99
    assert classify_molecules({}, ["A"])[0] == []


def test_cross_locus_mapping_ties_do_not_duplicate_a_molecule():
    a, b = context("a"), context("b", locus="L2")
    scores, _ = molecule_log_likelihoods([alignment(a), alignment(b)], [a, b], {"a": ["A"], "b": ["B"]})
    assert not scores


def database(tmp_path):
    root = tmp_path / "db"
    root.mkdir()
    loci = [Locus(locus_id=f"L{i}", forward_primer="ACG", reverse_primer="TTA",
        left_flank_sequence="TT", right_flank_sequence="CC", repeat_motif="GA",
        repeat_unit_length_bp=2, expected_min_repeats=1, expected_max_repeats=10,
    ) for i in range(1, 4)]
    for locus in loci:
        (root / f"{locus.locus_id}.fasta").write_text(
            ">A\nACGTTGAGACCTAA\n>B\nACGTTGAGAGAGAGAGAGAGAGACCTAA\n"
        )
    (root / "reference_metadata.tsv").write_text("reference_id\ttaxon_id\ttaxon_name\nA\t1\tTaxon A\nB\t2\tTaxon B\n")
    return root, loci


def test_assembly_mapping_and_profile_tree_need_no_phylogenetic_executables(tmp_path):
    root, loci = database(tmp_path)
    result = run_mapping_classification(database_path=root, loci=loci, outdir=tmp_path / "out",
        sample_id="sample", query_sequences={l.locus_id: "ACGTTGAGACCTAA" for l in loci},
        query_repeat_counts={l.locus_id: 2 for l in loci})
    matches = read_tsv(result["mapping_reference_matches"])
    assert matches[0]["reference_id"] == "A"
    assert matches[0]["distance"] == "-0.0"
    summary = read_tsv(result["taxonomic_identification"])[0]
    assert summary["best_taxon"] == "1"
    assert summary["assignment_status"] == "SUPPORTED"
    tree = result["mlva_profile_tree"].read_text()
    assert tree.endswith(";\n") or tree.endswith(";")
    pairs = _tip_patristic_distances(_parse_newick(tree))
    assert pairs[("A", "sample")] == pytest.approx(0)
    metadata = read_tsv(result["mlva_profile_tree_metadata"])
    assert {r["sample_id"] for r in metadata} == {"A", "B", "sample"}
    # Missing query calls are not imputed to zero merely to obtain a tree.
    no_tree = run_mapping_classification(database_path=root, loci=loci, outdir=tmp_path / "out",
        sample_id="sample", query_sequences={"L1": "ACGTTGAGACCTAA"}, query_repeat_counts={"L1": 2})
    assert "mlva_profile_tree" not in no_tree
    assert not result["mlva_profile_tree"].exists()


def test_raw_read_mapping_preserves_alternatives_and_low_coverage(tmp_path, monkeypatch):
    root, loci = database(tmp_path)
    seen = {}
    def fake_mapping(reference, reads1, reads2, contexts, bam, threads, technology, **kwargs):
        seen.update(kwargs)
        # One molecule per locus still supports classification; no coverage gate.
        rows = []
        for c in contexts:
            if c.repeat_count == 2:
                rows.append(alignment(c, molecule=c.locus_id))
            else:
                delta = len(c.sequence) - len("ACGTTGAGACCTAA")
                rows.append(alignment(c, molecule=c.locus_id, cs=f":5-{'GA' * (delta // 2)}:9"))
        return rows
    monkeypatch.setattr("mlvamaps.mapping_classification.map_reads_to_candidates_bam", fake_mapping)
    result = run_mapping_classification(database_path=root, loci=loci, outdir=tmp_path / "out",
        sample_id="reads", reads1="raw.fastq", technology="illumina",
        locus_quality={l.locus_id: {"depth": 1, "status": "low_coverage"} for l in loci})
    assert seen["retain_all_competitors"]
    assert seen["include_unknown_repeats"]
    matches = read_tsv(result["mapping_reference_matches"])
    assert matches[0]["reference_id"] == "A"
    assert matches[0]["missing_locus_gate_passed"] == "no"
    assert float(matches[0]["missing_locus_penalty"]) == 0
    assert read_tsv(result["taxonomic_identification"])[0]["best_taxon"] == "1"
    assert json.loads(result["classification_details"].read_text())["molecules"] == 3


def test_real_minimap_low_coverage_and_mixture(tmp_path):
    import os
    import random
    import shutil

    executable = os.environ.get("MINIMAP2_BIN") or shutil.which("minimap2")
    if not executable:
        pytest.skip("minimap2 is not installed")
    rng = random.Random(17)
    root = tmp_path / "observed"
    root.mkdir()
    loci, products = [], {}
    for i in range(3):
        left = "".join(rng.choices("ACGT", k=80))
        right = "".join(rng.choices("ACGT", k=80))
        locus = Locus(locus_id=f"L{i}", forward_primer="ACG", reverse_primer="TTA",
            left_flank_sequence=left, right_flank_sequence=right,
            repeat_motif="GA", repeat_unit_length_bp=2,
            expected_min_repeats=1, expected_max_repeats=12)
        loci.append(locus)
        a, b = "ACG" + left + "GA" * 4 + right + "TAA", "ACG" + left + "GA" * 10 + right + "TAA"
        products[locus.locus_id] = (a, b)
        (root / f"{locus.locus_id}.fasta").write_text(f">A\n{a}\n>B\n{b}\n")
    (root / "reference_metadata.tsv").write_text("reference_id\ttaxon_id\ttaxon_name\nA\t1\tTaxon A\nB\t2\tTaxon B\n")
    for depth in (1, 3):
        reads = tmp_path / f"depth{depth}.fastq"
        reads.write_text("".join(f"@{locus}_{i}\n{pair[0]}\n+\n{'I' * len(pair[0])}\n" for locus, pair in products.items() for i in range(depth)))
        result = run_mapping_classification(database_path=root, loci=loci, outdir=tmp_path / f"out{depth}",
            sample_id=f"depth{depth}", reads1=reads, minimap2_bin=executable, technology="hifi")
        summary = read_tsv(result["taxonomic_identification"])[0]
        assert summary["best_taxon"] == "1"
        assert float(summary["model_support"]) > 0.99
    mixed = tmp_path / "mixed.fastq"
    mixed.write_text("".join(f"@{locus}_{i}\n{seq}\n+\n{'I' * len(seq)}\n" for locus, pair in products.items() for i, seq in enumerate([pair[0]] * 7 + [pair[1]] * 3)))
    result = run_mapping_classification(database_path=root, loci=loci, outdir=tmp_path / "mixture",
        sample_id="mixed", reads1=mixed, minimap2_bin=executable, technology="hifi", sample_mode="metagenome")
    matches = read_tsv(result["mapping_reference_matches"])
    fractions = {row["reference_id"]: float(row["locus_balanced_fraction"]) for row in matches}
    assert fractions["A"] == pytest.approx(0.7, abs=0.03)
    assert fractions["B"] == pytest.approx(0.3, abs=0.03)
    assert read_tsv(result["taxonomic_identification"])[0]["assignment_status"] == "MIXED_REFERENCES"


def test_repeat_sequence_changes_with_equal_length_are_not_erased():
    c = context("a")
    row = alignment(c, cs=":6+ga:2-ga:10")
    stats = alignment_statistics(row, c)
    assert stats["repeat_delta"] == 0
    assert stats["insertions"] == stats["deletions"] == 2
    scores, _ = molecule_log_likelihoods([row], [c], {"a": ["A"]})
    assert scores[("L1", "m1")]["A"] < 0


def test_mapping_report_and_profile_match_export_keep_new_score_semantics(tmp_path):
    from mlvamaps.report import write_assembly_report
    from mlvamaps.profile_matching import sequence_reference_match_rows

    root, loci = database(tmp_path)
    result = run_mapping_classification(database_path=root, loci=loci, outdir=tmp_path / "out",
        sample_id="sample", query_sequences={l.locus_id: "ACGTTGAGACCTAA" for l in loci},
        query_repeat_counts={l.locus_id: 2 for l in loci})
    matches = read_tsv(result["mapping_reference_matches"])
    normalized = sequence_reference_match_rows(matches)
    assert normalized[0]["match_type"] == "mapping_reference"
    assert normalized[0]["distance"] == matches[0]["distance"]
    write_assembly_report(tmp_path / "out", "sample", [], [], loci=loci, phylogenetic_rows=matches)
    report = (tmp_path / "out" / "report.html").read_text()
    assert "Original mapping evidence" in report
    assert "classification/mlva_profiles.tree" in report
    assert "weighted nearest-reference averages" not in report
    assert "Technical marker-distance components" not in report


@pytest.mark.parametrize("same_taxon", [False, True])
def test_ambiguous_reference_group_does_not_force_a_reference(tmp_path, same_taxon):
    root, loci = database(tmp_path)
    if same_taxon:
        (root / "reference_metadata.tsv").write_text("reference_id\ttaxon_id\ttaxon_name\nA\t1\tTaxon A\nB\t1\tTaxon A\n")
    for locus in loci:
        (root / f"{locus.locus_id}.fasta").write_text(
            ">A\nACGTTGAGACCTAA\n>B\nACGTTGAGACCTAA\n"
        )
    result = run_mapping_classification(database_path=root, loci=loci, outdir=tmp_path / "out",
        sample_id="sample", query_sequences={l.locus_id: "ACGTTGAGACCTAA" for l in loci})
    summary = read_tsv(result["taxonomic_identification"])[0]
    assert summary["assignment_status"] == "AMBIGUOUS_REFERENCES"
    assert summary["assignment"] == ("Taxon A" if same_taxon else "Ambiguous taxa: Taxon A; Taxon B")
    assert summary["equivalent_references"] == "A;B"
    assert summary["best_taxon"] == ("1" if same_taxon else "")
    assert float(summary["unclassified_fraction"]) < 0.01
    assert read_tsv(result["mapping_reference_matches"])[0]["equivalent_references"] == "A;B"


def test_zero_repeat_reference_is_a_measured_allele_not_missing_data():
    zero = context("zero", 0)
    row = alignment(zero, cs=":6+gagagaga:6")
    stats = alignment_statistics(row, zero)
    assert stats["repeat_delta"] == 4
    assert stats["insertions"] == 0
    scores, _ = molecule_log_likelihoods([row], [zero], {"zero": ["Z"]})
    assert scores[("L1", "m1")]["Z"] == -4


def test_coverage_penalty_is_reference_specific_and_applied_once(tmp_path, monkeypatch):
    root, loci = database(tmp_path)
    (root / "L3.fasta").write_text(">A\nACGTTGAGACCTAA\n")
    def fake_mapping(reference, reads1, reads2, contexts, bam, threads, technology, **kwargs):
        return [alignment(c, molecule=c.locus_id) for c in contexts if c.locus_id != "L3"]
    monkeypatch.setattr("mlvamaps.mapping_classification.map_reads_to_candidates_bam", fake_mapping)
    quality = {l.locus_id: {"depth": 1, "status": "called"} for l in loci[:2]}
    quality["L3"] = {"depth": 0, "status": "not_found"}
    kwargs = dict(database_path=root, loci=loci, sample_id="reads", reads1="raw.fastq",
        technology="illumina", locus_quality=quality, missing_locus_min_fraction=2 / 3)
    low = run_mapping_classification(**kwargs, outdir=tmp_path / "low")
    low_groups = json.loads(low["classification_details"].read_text())["groups"]
    assert next(g for g in low_groups if "A" in g["references"])["references"] == ["A", "B"]
    for locus in ("L1", "L2"):
        quality[locus]["depth"] = 3
    high = run_mapping_classification(**kwargs, outdir=tmp_path / "high")
    matches = {r["reference_id"]: r for r in read_tsv(high["mapping_reference_matches"])}
    assert float(matches["A"]["log_likelihood"]) == -1
    assert float(matches["B"]["log_likelihood"]) == 0
    assert float(matches["A"]["missing_locus_penalty"]) == 1
    assert float(matches["B"]["missing_locus_penalty"]) == 0


@pytest.mark.parametrize("sample_mode", ["isolate", "metagenome"])
def test_identification_follows_reference_not_species_total(tmp_path, monkeypatch, sample_mode):
    from mlvamaps.report import write_assembly_report

    root, loci = database(tmp_path)
    for locus in loci:
        with (root / f"{locus.locus_id}.fasta").open("a") as handle:
            handle.write(">C\nACGTTGAGAGACCTAA\n")
    (root / "reference_metadata.tsv").write_text(
        "reference_id\ttaxon_id\ttaxon_name\nA\t1\tTaxon A\nB\t2\tTaxon B\nC\t2\tTaxon B\n")
    if sample_mode == "isolate":
        likelihoods = {(locus.locus_id, locus.locus_id): {
            "A": np.log(0.4) / 3, "B": np.log(0.35) / 3,
            "C": np.log(0.25) / 3, UNKNOWN: -100,
        } for locus in loci}
    else:
        likelihoods = {(locus.locus_id, str(i)): {
            ref: 0.0 if ref == source else -100.0 for ref in ("A", "B", "C", UNKNOWN)
        } for locus in loci for i, source in enumerate(["A"] * 4 + ["B"] * 3 + ["C"] * 3)}
    monkeypatch.setattr("mlvamaps.mapping_classification.molecule_log_likelihoods",
                        lambda *args: (likelihoods, {}))
    result = run_mapping_classification(database_path=root, loci=loci, outdir=tmp_path / "out",
        sample_id="sample", sample_mode=sample_mode)
    matches = read_tsv(result["mapping_reference_matches"])
    summary = read_tsv(result["taxonomic_identification"])[0]
    evidence = read_tsv(result["taxonomic_identification_evidence"])
    assert summary["assignment"] == "Taxon A"
    assert summary["reference_id"] == summary["closest_reference"] == matches[0]["reference_id"] == "A"
    assert summary["best_taxon"] == "1"
    assert float(summary["model_support"]) == pytest.approx(0.4)
    assert [r["reference_id"] for r in evidence] == [r["reference_id"] for r in matches]
    assert sum(float(r["model_support"]) for r in evidence if r["taxon_id"] == "2") == pytest.approx(0.6)
    assert summary["assignment_status"] == ("MIXED_REFERENCES" if sample_mode == "metagenome" else "CLOSEST_REFERENCE_LOW_CONFIDENCE")
    write_assembly_report(tmp_path / "out", "sample", [], [], loci=loci, phylogenetic_rows=matches)
    report = (tmp_path / "out" / "report.html").read_text()
    assert '<div class="taxon-call">Taxon A</div>' in report
    assert "Closest reference IDs" in report
    assert "<th>Reference IDs</th><th>Taxon annotation</th>" in report
    assert "support is never summed across a species" in report


def test_reference_identification_without_taxon_metadata(tmp_path):
    root, loci = database(tmp_path)
    (root / "reference_metadata.tsv").unlink()
    result = run_mapping_classification(database_path=root, loci=loci, outdir=tmp_path / "out",
        sample_id="sample", taxon_identification=True,
        query_sequences={l.locus_id: "ACGTTGAGACCTAA" for l in loci})
    summary = read_tsv(result["taxonomic_identification"])[0]
    assert summary["assignment"] == "Taxonomy unavailable"
    assert summary["reference_id"] == "A"
    assert summary["best_taxon"] == summary["best_species"] == ""
    assert summary["assignment_status"] == "SUPPORTED"
    assert float(summary["unclassified_fraction"]) < 0.01
