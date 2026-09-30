"""Guard redundant work without machine-dependent timing thresholds."""

from dataclasses import replace
import random
from types import SimpleNamespace

from mlvamaps.alignment_evidence import CandidateAlignment
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.locus_measurement import find_anchor
from mlvamaps.models import Locus
from mlvamaps.short_read_evidence import extract_short_read_evidence


def test_impossible_anchors_do_not_start_native_searches(monkeypatch):
    def unexpected_searcher(*args):
        raise AssertionError("an impossible anchor launched a native search")

    monkeypatch.setattr("mlvamaps.locus_measurement._sassy_searcher", unexpected_searcher)
    pattern = "ACGTACGACAGTCGATGCAA"
    assert find_anchor(pattern, "T" * 100, 3) is None
    assert find_anchor(pattern, "ACGT", 3) is None
    assert find_anchor(pattern, "T" * 100, 0) is None
    assert find_anchor(pattern, pattern, 3, search_start=5, search_end=10) is None


def test_anchor_prefilter_retains_substitutions_insertions_deletions_and_iupac(monkeypatch):
    pattern = "ACGTACGACAGTCGATGCAA"
    rng = random.Random(12)
    calls = []

    def search(anchor, text, k):
        calls.append((anchor, text, k))
        return [SimpleNamespace(
            text_start=2, text_end=len(text) - 2, cost=k,
            pattern_start=0, pattern_end=len(anchor), cigar=f"{len(anchor)}M",
        )]

    monkeypatch.setattr(
        "mlvamaps.locus_measurement._sassy_searcher",
        lambda *args: SimpleNamespace(search=search),
    )
    # Edit locations include partition boundaries and the ends of the anchor.
    for edits in (1, 2, 3):
        for _ in range(80):
            observed = list(pattern)
            for _ in range(edits):
                position = rng.randrange(len(observed))
                operation = rng.choice(("insert", "delete", "substitute"))
                if operation == "insert":
                    observed.insert(position, rng.choice("ACGT"))
                elif operation == "delete":
                    observed.pop(position)
                else:
                    observed[position] = rng.choice("ACGT")
            assert find_anchor(pattern, "GG" + "".join(observed) + "CC", edits) is not None
    assert calls
    calls.clear()
    assert find_anchor("ACGNACGACAGTCGATGCAA", "GG" + pattern + "CC", 1) is not None
    assert len(calls) == 1  # Ambiguous anchors still go to the native matcher.


def test_measurements_are_reused_without_collapsing_molecules_or_quality(monkeypatch):
    calls = []

    def measure(sequence, locus, quality, *, source, round_tolerance):
        calls.append((sequence, locus.locus_id, quality, round_tolerance))
        allele = (3 if locus.locus_id == "L1" else 4) if round_tolerance == .25 else 3.5
        return SimpleNamespace(
            called_allele=allele if quality.startswith("H") else None,
            status="FULL_PRODUCT" if quality.startswith("H") else "PRESENCE_ONLY",
            confidence=.9,
        )

    monkeypatch.setattr("mlvamaps.short_read_evidence.measure_locus_product", measure)
    loci = [Locus("L1"), Locus("L2")]
    contexts = [
        CandidateContext(f"{locus.locus_id}_{repeat}", locus.locus_id,
                         "A" * 30, repeat, 4, 10, 22, "", "")
        for locus in loci for repeat in (3, 4, 5)
    ]
    prototype = CandidateAlignment(
        molecule_id="", read_id="", mate=None, locus_id="", candidate_id="",
        repeat_count=3, reference_id="ref", alignment_score=30, mapping_quality=0,
        alignment_identity=1, query_coverage=1, reference_coverage=.4,
        query_start=0, query_end=12, reference_start=0, reference_end=12,
        cigar="12M", cs="", primary=True, secondary=False, reverse=False,
        query_sequence="AACCGGTTAACC", query_quality="I" * 12,
    )
    alignments = [
        replace(prototype, molecule_id=f"m{index}", read_id=f"m{index}",
                locus_id=context.locus_id, candidate_id=context.candidate_id,
                repeat_count=context.repeat_count, query_quality=("I" if index < 2 else "H") * 12)
        for context in contexts for index in range(4) for _duplicate in range(2)
    ]
    evidence = extract_short_read_evidence(alignments, contexts, loci)
    assert len(calls) == 8  # 2 loci x 2 qualities x 2 orientations, including misses.
    assert len(evidence) == 24
    for row in evidence:
        expected = None if row.molecule_id in {"m0", "m1"} else 3 if row.locus_id == "L1" else 4
        assert row.measured_repeat_count == expected
    assert {row.molecule_id for row in evidence if row.direct_product_measurement} == {"m2", "m3"}
    second = extract_short_read_evidence(alignments, contexts, loci, round_tolerance=.4)
    assert len(calls) == 16  # No cross-sample or cross-setting cached results.
    assert {row.measured_repeat_count for row in second} == {None, 3.5}
