from dataclasses import replace
import shutil

import pytest

from mlvamaps.alignment_evidence import CandidateAlignment, CandidateEvidence
from mlvamaps.allele_inference import InferenceThresholds
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.mapped_reconstruction import aligned_flanks, recruit_mapping_fragments, resolve_query_products
from mlvamaps.models import Locus, ReadRecord


def mapped_fixture():
    import random
    from mlvamaps.sequence import revcomp
    rng = random.Random(402)
    forward, left, right, reverse = ["".join(rng.choices("ACGT", k=n)) for n in (20, 40, 40, 20)]
    product = forward + left + "ATGC" * 5 + right + revcomp(reverse)
    locus = Locus("L", forward_primer=forward, reverse_primer=reverse,
                  left_flank_sequence=left, right_flank_sequence=right, repeat_motif="ATGC",
                  repeat_unit_length_bp=4, expected_product_size_bp=len(product), nominal_repeat_units=5,
                  expected_min_repeats=1, expected_max_repeats=10)
    ctx = CandidateContext("c", "L", product, 5, 4, 60, 80, product[:60], product[80:])
    row = alignment(mate=None, query_sequence=product, query_quality="I" * len(product),
                    query_end=len(product), reference_end=len(product), cigar=f"{len(product)}M")
    return locus, ctx, row


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


def test_repeat_expansion_is_not_rejected_as_low_identity():
    locus, ctx, row = mapped_fixture()
    expanded = replace(row, query_sequence=ctx.sequence[:60] + "ATGC" * 8 + ctx.sequence[80:],
                       query_end=152, cigar="60M12I80M", alignment_identity=136 / 152)
    assert "m" in recruit_mapping_fragments([expanded], [ctx])["L"]
    assert not recruit_mapping_fragments([replace(row, alignment_identity=.88)], [ctx])


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
    def pcr(queries, sample, tolerance, max_primer_mismatches=2):
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


def test_guided_window_ignores_remote_pcr_decoy_and_preserves_orientation_and_insertions():
    from mlvamaps.locus_measurement import mapped_locus_window, measure_locus_product
    from mlvamaps.sequence import revcomp
    locus, ctx, row = mapped_fixture()
    # A valid, different PCR product occurs elsewhere on this read.
    prefix = ctx.sequence[:60] + "ATGC" * 2 + ctx.sequence[80:] + "C" * 500
    product = ctx.sequence[:60] + "ATGC" * 8 + ctx.sequence[80:]
    original = prefix + product
    for reverse in (False, True):
        mapped = replace(row, reverse=reverse, query_sequence=revcomp(original) if reverse else original,
                         query_quality="I" * len(original), query_start=len(prefix), query_end=len(original),
                         cigar=f"{len(prefix)}S60M12I80M")
        sequence, quality, offset, spanning = mapped_locus_window(mapped, ctx)
        assert spanning and offset > len(prefix) - 30
        assert product in sequence
        assert measure_locus_product(sequence, locus, quality).called_allele == 8
    partial = replace(row, query_sequence=row.query_sequence[:-30], query_quality=None,
                      query_end=110, reference_end=110, cigar="110M")
    assert not mapped_locus_window(partial, ctx)[3]
    assert not mapped_locus_window(replace(row, cigar="10H140M"), ctx)[3]


def test_guided_consensus_corrects_primers_before_pcr_without_creating_molecules():
    from mlvamaps.long_read_evidence import extract_long_read_evidence
    from mlvamaps.allele_inference import infer_alleles
    from mlvamaps.mapped_reconstruction import apply_reconstructed_calls
    locus, ctx, row = mapped_fixture()
    reads = []
    for i in range(10):
        sequence = list(row.query_sequence)
        for position in ((i + 5 * j) % 20 for j in range(4)):
            sequence[position] = "ACGT"[("ACGT".index(sequence[position]) + 1) % 4]
        reads.append(replace(row, molecule_id=f"read{i}", read_id=f"read{i}", query_sequence="".join(sequence)))
    evidence = extract_long_read_evidence(reads, [ctx], [locus], "ont")
    assert not any(item.direct_product_measurement for item in evidence)
    thresholds = InferenceThresholds()
    calls, molecules = infer_alleles(evidence, [locus], [ctx], "s", "ont", thresholds,
                                    require_measurement=True, sample_mode="isolate")
    products, audit, updated, changed = resolve_query_products(calls, evidence, [locus], [], "s", thresholds, .25)
    assert products == {"L": (ctx.sequence, 5)}
    assert audit[0]["method"] == "mapped_flank_consensus"
    assert audit[0]["reconstruction_molecules"] == 10
    assert audit[0]["supporting_complete_molecules"] == 0
    assert not changed and updated == evidence
    before = calls[0]["molecule_support"]
    apply_reconstructed_calls(calls, audit, thresholds)
    assert calls[0]["repeat_count"] == 5
    assert calls[0]["molecule_support"] == before
    assert len(molecules) == 10
    calls[0].update(status="called", repeat_count=6)
    apply_reconstructed_calls(calls, audit, thresholds)
    assert calls[0]["status"] == "ambiguous" and calls[0]["repeat_count"] == ""
    calls[0].update(status="mixed", repeat_count=6)
    apply_reconstructed_calls(calls, audit, thresholds)
    assert calls[0]["status"] == "mixed" and calls[0]["repeat_count"] == 6


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"), reason="native tools unavailable")
@pytest.mark.parametrize("technology", ["long", "short"])
def test_primer_error_reconstruction_reaches_public_calls(tmp_path, monkeypatch, technology):
    import csv
    from mlvamaps.io import read_fasta, write_fastq, write_tsv
    from mlvamaps.pipeline import run_call
    from mlvamaps.short_reads import run_short_read_call
    locus, ctx, row = mapped_fixture()
    records = []
    for i in range(10):
        sequence = list(row.query_sequence)
        for position in ((i + 5 * j) % 20 for j in range(4)):
            sequence[position] = "ACGT"[("ACGT".index(sequence[position]) + 1) % 4]
        records.append(ReadRecord(f"read{i}", "".join(sequence), row.query_quality))
    reads = tmp_path / "reads.fastq"
    write_fastq(records, reads)
    panel = tmp_path / "panel.tsv"
    write_tsv([vars(locus)], panel, list(vars(locus)))
    database = tmp_path / "database"
    database.mkdir()
    (database / "L.fasta").write_text(f">reference\n{ctx.sequence}\n")
    def duplicate_mapping(*args, **kwargs):
        pytest.fail("Recruitment must reuse the first BAM")
    monkeypatch.setattr("mlvamaps.recruitment.run_minimap2_command", duplicate_mapping)
    kwargs = dict(loci_path=str(panel), database_path=str(database), outdir=str(tmp_path / "out"),
                  sample_id="s", threads=2, show_progress=False, sample_mode="isolate")
    if technology == "long":
        paths = run_call(reads_path=str(reads), **kwargs)
    else:
        paths = run_short_read_call(reads1_path=str(reads), reads2_path=None, **kwargs)
    with paths["calls"].open() as handle:
        call = next(csv.DictReader(handle, delimiter="\t"))
    assert call["status"] == "PASS" and float(call["repeat_count"]) == 5
    assert dict(read_fasta(tmp_path / "out/query_loci.fasta")) == {"L": ctx.sequence}
    with (tmp_path / "out/query_locus_reconstruction.tsv").open() as handle:
        audit = next(csv.DictReader(handle, delimiter="\t"))
    assert audit["method"] == "mapped_flank_consensus" and audit["exported"] == "yes"
