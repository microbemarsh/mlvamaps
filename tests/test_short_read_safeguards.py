from dataclasses import replace
import gzip

import pytest

from mlvamaps.alignment_evidence import CandidateAlignment, CandidateEvidence
from mlvamaps.allele_inference import infer_alleles, estimate_unresolved_calls
from mlvamaps.candidate_contexts import CandidateContext
from mlvamaps.io import read_fastq
from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.short_read_evidence import estimate_insert_distribution, extract_short_read_evidence
from mlvamaps.short_read_fragments import RecruitedFragment


def test_insert_calibration_uses_independent_flank_pairs_and_rejects_repeat_geometry():
    context = CandidateContext("c", "L", "A" * 2000, 5, 4, 1000, 1020, "", "")
    first = CandidateAlignment("", "", 1, "L", "c", 5, "ref", 100, 60, 1, 1, 1,
                               0, 100, 200, 300, "100M", "", True, False, False)
    rows = [replace(first, molecule_id=str(i)) for i in range(120)]
    rows += [replace(first, molecule_id=str(i), mate=2, reverse=True,
                     reference_start=400 + i % 20, reference_end=500 + i % 20) for i in range(120)]
    fit = estimate_insert_distribution(rows, [context])
    assert fit[0] == 309.5 and 5 < fit[1] < 7 and fit[2] == 120
    assert estimate_insert_distribution(rows[:10], [context]) is None
    assert estimate_insert_distribution(rows, [replace(context, repeat_start=350)]) is None
    assert estimate_insert_distribution([replace(row, mapping_quality=0) for row in rows], [context]) is None


def test_unspanned_long_array_reports_a_bound_and_never_copies_the_reference_count():
    locus = Locus("L", left_flank_sequence="GACTGACCTGATCGACGTTAGCTA", repeat_motif="ATGC", repeat_unit_length_bp=4)
    read = ReadRecord("m", locus.left_flank_sequence + "ATGC" * 12, "I" * 72)
    fragments = {"L": {"m": RecruitedFragment(ReadPair("m", read), True, False, False)}}
    evidence = extract_short_read_evidence([], [], [locus], recruited_fragments=fragments)
    calls, molecules = infer_alleles(evidence, [locus], [], "s", "illumina", require_measurement=True)
    estimate_unresolved_calls(calls, evidence, [])
    assert not molecules and calls[0]["repeat_count"] == ""
    assert calls[0]["repeat_lower_bound"] == 12
    assert evidence[0].repeat_count is None


def test_reconstructed_product_does_not_multiply_likelihood_or_classification_votes():
    def run(n):
        evidence = [CandidateEvidence("L", 5, str(i), measured_repeat_count=5, raw_repeat_count=5,
                                      repeat_sigma=.3, evidence_tier="microassembly") for i in range(n)]
        return infer_alleles(evidence, [Locus("L", repeat_unit_length_bp=4)], [], "s", "illumina", require_measurement=True)
    three, molecules = run(3)
    thirty, more = run(30)
    assert not molecules and not more
    assert three[0]["best_probability"] == thirty[0]["best_probability"]
    assert thirty[0]["derived_molecule_support"] == 30


def test_ambiguous_boundaries_survive_fragment_extraction():
    from mlvamaps.sequence import revcomp
    locus = Locus("L", forward_primer="ACGTTAGCTAGTCGTACGAT", reverse_primer="TGCATCAGTCGATGCTAGCA",
        left_flank_sequence="GACTGACCTGATCGACGTTAGCTA", right_flank_sequence="CCGTAACCGGTTAGGCTACTGAC",
        repeat_motif="ATGC", repeat_unit_length_bp=4)
    sequence = (locus.forward_primer + locus.left_flank_sequence + "ATGC" * 2
                + locus.left_flank_sequence + "ATGC" * 5 + locus.right_flank_sequence + revcomp(locus.reverse_primer))
    read = ReadRecord("m", sequence, "I" * len(sequence))
    fragments = {"L": {"m": RecruitedFragment(ReadPair("m", read), True, True, False)}}
    evidence = extract_short_read_evidence([], [], [locus], recruited_fragments=fragments)
    assert evidence[0].measured_repeat_count is not None
    assert "boundary_ambiguous" in evidence[0].quality_flags
    assert evidence[0].repeat_sigma >= 1


@pytest.mark.parametrize("content", [b"@\nACGT\n+\nIIII\n", b"@x\n\n+\n\n", b"@x\nACGT\n+\nIII\n"])
def test_malformed_records_fail_instead_of_creating_missing_genotypes(tmp_path, content):
    path = tmp_path / "reads.fastq"
    path.write_bytes(content)
    with pytest.raises(ValueError, match="FASTQ"):
        list(read_fastq(path))


def test_truncated_gzip_is_an_error_and_empty_fastq_has_zero_records(tmp_path):
    path = tmp_path / "reads.fastq.gz"
    path.write_bytes(gzip.compress(b"@x\nACGT\n+\nIIII\n")[:-8])
    with pytest.raises(ValueError, match="FASTQ"):
        list(read_fastq(path))
    path.write_bytes(gzip.compress(b""))
    assert list(read_fastq(path)) == []


def test_duplicate_fastq_ids_are_rejected_at_ingestion(tmp_path):
    from mlvamaps.io import read_fastq_pairs
    path = tmp_path / "reads.fastq"
    path.write_text("@m\nACGT\n+\nIIII\n@m\nTGCA\n+\nIIII\n")
    with pytest.raises(ValueError, match="Duplicate FASTQ molecule"):
        list(read_fastq_pairs(path))
