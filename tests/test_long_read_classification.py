"""Distinct regions of one long read can inform different loci."""
from dataclasses import replace

import pytest

from mlvamaps import mapping_classification
from mlvamaps.models import Locus
from test_mapping_classification import alignment, context


def positioned_alignment(candidate, start, *, reverse=False, windowed=False, score=100, mate=None):
    row = alignment(candidate, score=score, mate=mate)
    first = 1000 - start - len(candidate.sequence) if reverse else start
    offset = first if windowed else 0
    return replace(row, query_start=first - offset, query_end=first + len(candidate.sequence) - offset,
                   query_offset=offset, query_length=1000, reverse=reverse)


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("windowed", [False, True])
def test_disjoint_long_read_loci_keep_independent_likelihoods(reverse, windowed):
    a, b = context("a"), context("b", locus="L2")
    rows = [positioned_alignment(a, 10), positioned_alignment(b, 100, reverse=reverse, windowed=windowed)]
    members = {"a": ["A"], "b": ["B"]}
    scores, _ = mapping_classification.molecule_log_likelihoods(rows, [a, b], members, allow_multi_locus=True)
    assert set(scores) == {("L1", "m1"), ("L2", "m1")}
    assert scores[("L1", "m1")]["A"] == scores[("L2", "m1")]["B"] == 0
    duplicate_scores, _ = mapping_classification.molecule_log_likelihoods(
        rows * 3, [a, b], members, allow_multi_locus=True)
    assert duplicate_scores == scores
    assert not mapping_classification.molecule_log_likelihoods(rows, [a, b], members)[0]


@pytest.mark.parametrize("first_score, expected", [(100, {"L3"}), (110, {"L1", "L3"})])
def test_overlapping_cross_locus_competitors_still_share_one_assignment(first_score, expected):
    a, b, c = context("a"), context("b", locus="L2"), context("c", locus="L3")
    rows = [positioned_alignment(a, 10, score=first_score),
            positioned_alignment(b, 10, reverse=True, windowed=True), positioned_alignment(c, 100)]
    scores, _ = mapping_classification.molecule_log_likelihoods(
        rows, [a, b, c], {"a": ["A"], "b": ["B"], "c": ["C"]}, allow_multi_locus=True)
    assert {locus for locus, _ in scores} == expected


def test_a_repeated_reference_match_cannot_hide_a_cross_locus_tie_in_another_window():
    a, b = context("a"), context("b", locus="L2")
    rows = [positioned_alignment(a, 10, score=110), positioned_alignment(a, 100),
            positioned_alignment(b, 100)]
    scores, _ = mapping_classification.molecule_log_likelihoods(
        rows, [a, b], {"a": ["A"], "b": ["B"]}, allow_multi_locus=True)
    assert set(scores) == {("L1", "m1")}


@pytest.mark.parametrize("score, expected", [(100, set()), (110, {("L1", "m1")})])
def test_paired_reads_keep_one_locus_assignment(score, expected):
    a, b = context("a"), context("b", locus="L2")
    rows = [positioned_alignment(a, 10, score=score, mate=1), positioned_alignment(b, 100, mate=2)]
    scores, _ = mapping_classification.molecule_log_likelihoods(rows, [a, b], {"a": ["A"], "b": ["B"]})
    assert set(scores) == expected


@pytest.mark.parametrize("technology, reads, expected", [("hifi", True, True), ("ont", True, True),
                                                       ("illumina", True, False), ("hifi", False, False)])
def test_classification_enables_multiple_loci_only_for_long_reads(tmp_path, monkeypatch, technology, reads, expected):
    a, b = context("a"), context("b", locus="L2")
    contexts, members = [a, b], {"a": ["R"], "b": ["R"]}
    monkeypatch.setattr(mapping_classification, "observed_reference_contexts",
                        lambda *args: (tmp_path, contexts, members))
    original = mapping_classification.molecule_log_likelihoods
    selected = []

    def record_option(alignments, contexts, members, repeat_scale=1, measured_repeats=None, allow_multi_locus=False):
        selected.append(allow_multi_locus)
        return original(alignments, contexts, members, repeat_scale, measured_repeats, allow_multi_locus)

    monkeypatch.setattr(mapping_classification, "molecule_log_likelihoods", record_option)
    mapping_classification.run_mapping_classification(
        database_path=tmp_path, loci=[Locus("L1"), Locus("L2")], outdir=tmp_path / "out", sample_id="s",
        technology=technology, reads1=tmp_path / "reads.fastq" if reads else None,
        read_alignments=[positioned_alignment(a, 10, mate=1), positioned_alignment(b, 100, mate=2)])
    assert selected == [expected]
