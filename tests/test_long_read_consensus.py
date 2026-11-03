"""Isolate consensus repairs small sizing noise before assembly PCR."""
from dataclasses import replace
import shutil

import pytest

from mlvamaps.alignment_evidence import CandidateEvidence
from mlvamaps.allele_inference import InferenceThresholds, infer_alleles
from mlvamaps.long_read_evidence import extract_long_read_evidence
from mlvamaps.mapped_reconstruction import resolve_query_products
from mlvamaps.models import Locus
from test_mapping_reconstruction import mapped_fixture


@pytest.mark.skipif(not shutil.which("sassy"), reason="native PCR required")
@pytest.mark.parametrize("sample_mode,secondary,recovers", [
    ("isolate", False, True), ("metagenome", False, False), ("isolate", True, False)])
def test_adjacent_sizing_noise_does_not_block_isolate_sequence_reconstruction(sample_mode, secondary, recovers):
    locus, context, prototype = mapped_fixture()
    reads = []
    for index in range(10):
        sequence = list(context.sequence)
        position = 25 + index
        sequence[position] = "ACGT"[("ACGT".index(sequence[position]) + 1) % 4]
        sequence = "".join(sequence)
        if index >= 6:
            position = 30 + index
            sequence = sequence[:position] + "A" + sequence[position:]
        cigar = "140M" if index < 6 else f"{position}M1I{140 - position}M"
        reads.append(replace(prototype, molecule_id=str(index), read_id=str(index),
            query_sequence=sequence, query_quality="I" * len(sequence), query_end=len(sequence), cigar=cigar))
    evidence = extract_long_read_evidence(reads, [context], [locus], "hifi")
    thresholds = InferenceThresholds()
    calls, molecules = infer_alleles(evidence, [locus], [context], "s", "hifi", thresholds,
                                    require_measurement=True, sample_mode="isolate")
    assert calls[0]["repeat_count"] == 5 and not calls[0]["possible_secondary"]
    calls[0]["possible_secondary"] = secondary
    products, audit, updated, changed = resolve_query_products(
        calls, evidence, [locus], [], "s", thresholds, .25, sample_mode=sample_mode)
    if recovers:
        assert products == {"L": (context.sequence, 5)}
        assert audit[0]["method"] == "mapped_flank_consensus"
        assert audit[0]["reconstruction_molecules"] == 10
    else:
        assert products == {}
    assert len(molecules) == 10 and updated == evidence and not changed


@pytest.mark.parametrize("conflict", ["flag", "products"])
def test_conflicting_molecule_products_never_seed_consensus(monkeypatch, conflict):
    calls = [{"locus": "L", "repeat_count": 5, "repeat_count_raw": 5,
              "status": "ambiguous", "possible_secondary": False}]
    evidence = [CandidateEvidence("L", 5, str(i), direct_product_measurement=True,
        measured_repeat_count=5, raw_repeat_count=5, repeat_sigma=.1,
        metadata={"product_sequence": "ACGT"}) for i in range(3)]
    if conflict == "flag":
        evidence = [replace(row, quality_flags=("conflicting_measurements",)) for row in evidence]
    else:
        evidence += [replace(row, metadata={"product_sequence": "ACGA"}) for row in evidence]
    monkeypatch.setattr("mlvamaps.mapped_reconstruction.call_query_products", lambda queries, *args:
        pytest.fail("Conflicting molecule products reached PCR") if queries else {})
    products, audit, _updated, _changed = resolve_query_products(
        calls, evidence, [Locus("L")], [], "s", InferenceThresholds(), .25, sample_mode="isolate")
    assert not products and audit[0]["reconstruction_molecules"] == 0

