from dataclasses import replace

from mlvamaps.alignment_evidence import CandidateAlignment, CandidateEvidence
from mlvamaps.allele_inference import InferenceThresholds
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.mapped_reconstruction import aligned_flanks, recruit_mapping_fragments, resolve_query_products
from mlvamaps.models import Locus, ReadRecord


def context(locus="L", identifier="c"):
    return CandidateContext(identifier, locus, "A" * 150, 5, 10, 50, 100, "A" * 50, "A" * 50)


def alignment(molecule="m", **kwargs):
    row = CandidateAlignment(molecule, molecule + "/1", 1, "L", "c", 5, "ref", 100, 0,
                             1, 1, 2 / 3, 0, 100, 0, 100, "100M", "", True, False, False,
                             query_sequence="A" * 100, query_quality="I" * 100)
    return replace(row, **kwargs)


def test_mapping_recruitment_rescues_unmapped_mate_and_keeps_equal_alleles_within_locus():
    contexts = [context(), context(identifier="alternative")]
    rows = [alignment(), alignment(candidate_id="alternative", secondary=True, primary=False)]
    mate = ReadRecord("m/2", "C" * 80, "H" * 80)
    pools = recruit_mapping_fragments(rows, contexts, unmapped_mates={("m", 2): mate})
    assert set(pools["L"]) == {"m"}
    assert pools["L"]["m"].pair.read2 == mate
    assert pools["L"]["m"].rescued
    assert pools["L"]["m"].source == "mapping"


def test_mapping_recruitment_rejects_cross_locus_ties_conflicting_mates_and_repeat_only():
    contexts = [context(), context("other", "other")]
    other = alignment(locus_id="other", candidate_id="other")
    assert not recruit_mapping_fragments([alignment(), other], contexts)
    assert not recruit_mapping_fragments([alignment(), replace(other, mate=2)], contexts)
    repetitive = alignment(reference_start=50, reference_end=100, query_end=50, cigar="50M")
    assert not recruit_mapping_fragments([repetitive], contexts)
    deletion = alignment(reference_start=0, reference_end=100, query_end=50, cigar="50D50M")
    assert aligned_flanks(deletion, contexts[0]) == (0, 0)
    assert not recruit_mapping_fragments([deletion], contexts)


def test_bam_parser_recovers_only_mates_of_mapped_targets(tmp_path):
    from mlvamaps.minimap_mapping import parse_candidate_alignments
    sam = tmp_path / "mapped.sam"
    sam.write_text("@HD\tVN:1.6\n@SQ\tSN:c\tLN:150\n"
                   "m\t73\tc\t1\t60\t100M\t=\t1\t0\t" + "A" * 100 + "\t" + "I" * 100 + "\tNM:i:0\n"
                   "m\t133\tc\t1\t0\t*\t=\t1\t0\tCCCC\tHHHH\n"
                   "absent\t77\t*\t0\t0\t*\t*\t0\t0\tGGGG\tIIII\n")
    mates = {}
    rows = parse_candidate_alignments(sam, [context()], unmapped_mates=mates)
    assert len(rows) == 1
    assert mates == {("m", 2): ReadRecord("m", "CCCC", "HHHH")}


def test_query_reconstruction_requires_independent_sequence_support(monkeypatch):
    from mlvamaps import mapped_reconstruction
    calls = [{"locus": "L", "repeat_count": 5, "status": "called"}]
    loci = [Locus("L")]
    evidence = [CandidateEvidence("L", 5, "one", direct_product_measurement=True,
                                 measured_repeat_count=5, metadata={"product_sequence": "ACGT"})] * 10
    monkeypatch.setattr(mapped_reconstruction, "call_query_products", lambda queries, *args: {} if not queries else (_ for _ in ()).throw(AssertionError("unsupported PCR")))
    products, audit, _, _ = resolve_query_products(calls, evidence, loci, [], "s", InferenceThresholds(), .25)
    assert not products
    assert audit[0]["complete_molecules"] == 1
    # Readless reference candidates never become query products.
    assert not resolve_query_products(calls, [], loci, [], "s", InferenceThresholds(), .25)[0]


def test_query_reconstruction_uses_assembly_calibration_on_observed_sequence(monkeypatch):
    from mlvamaps import mapped_reconstruction
    sequence = "ACGTACGT"
    calls = [{"locus": "L", "repeat_count": 5, "status": "called"}]
    evidence = [CandidateEvidence("L", 5, str(i), direct_product_measurement=True,
                                 measured_repeat_count=5, raw_repeat_count=5,
                                 metadata={"product_sequence": sequence}) for i in range(4)]
    def pcr(queries, sample, tolerance):
        if not queries:
            return {}
        assert queries[0][1] == sequence
        return {"L": ({"status": "PASS", "repeat_count": 6, "repeat_count_raw": "6"}, sequence)}
    monkeypatch.setattr(mapped_reconstruction, "call_query_products", pcr)
    products, audit, updated, changed = resolve_query_products(
        calls, evidence, [Locus("L")], [], "s", InferenceThresholds(), .25)
    assert products == {"L": (sequence, 6)}
    assert changed and all(row.measured_repeat_count == 6 for row in updated)
    calls[0]["status"] = "mixed"
    assert not resolve_query_products(calls, evidence, [Locus("L")], [], "s", InferenceThresholds(), .25)[0]
