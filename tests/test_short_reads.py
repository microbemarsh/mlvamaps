from __future__ import annotations

import csv
import shutil
from pathlib import Path

import pytest

from mlvamaps.io import normalize_read_id, read_fastq_pairs
from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.sample_metadata import myoga_sample_row, normalize_metadata_row
from mlvamaps.sequence import revcomp
from mlvamaps.short_reads import merge_read_pair, qc_read_pairs, run_short_read_call
from mlvamaps.validation import compare_call_sets, summarize_validation


def _write_fastq(path: Path, records: list[tuple[str, str]], quality: str = "I") -> None:
    path.write_text("".join(
        f"@{name}\n{sequence}\n+\n{quality * len(sequence)}\n"
        for name, sequence in records
    ))


def _locus() -> Locus:
    return Locus(
        "L1", forward_primer="ACGTTGCAACGTTGCAACGT",
        reverse_primer="AGTCAGTCAGTCAGTCAGTC",
        left_flank_sequence="CCGGAATTCGACCTGA",
        right_flank_sequence="TTAACCGGCTAGGTCA", repeat_motif="ATGC",
        repeat_unit_length_bp=4, expected_min_repeats=2,
        expected_max_repeats=10, nominal_repeat_units=4,
    )


def _product(repeats: int = 4) -> str:
    locus = _locus()
    return (locus.forward_primer + locus.left_flank_sequence
            + locus.repeat_motif * repeats + locus.right_flank_sequence
            + revcomp(locus.reverse_primer))


def _write_panel(path: Path) -> None:
    locus = _locus()
    path.write_text(
        "locus_id\tforward_primer\treverse_primer\tleft_flank_sequence\t"
        "right_flank_sequence\trepeat_motif\trepeat_unit_length_bp\t"
        "expected_min_repeats\texpected_max_repeats\tnominal_repeat_units\n"
        f"L1\t{locus.forward_primer}\t{locus.reverse_primer}\t"
        f"{locus.left_flank_sequence}\t{locus.right_flank_sequence}\tATGC\t4\t2\t10\t4\n"
    )


def test_read_id_normalization_and_paired_fastq_preserve_mates(tmp_path):
    reads1, reads2 = tmp_path / "r1.fastq", tmp_path / "r2.fastq"
    _write_fastq(reads1, [("SRR1.1/1", "ACGT"), ("SRR1.2/1", "TGCA")])
    _write_fastq(reads2, [("SRR1.1/2", "TGCA"), ("SRR1.2/2", "ACGT")])
    pairs = list(read_fastq_pairs(reads1, reads2))
    assert [pair.molecule_id for pair in pairs] == ["SRR1.1", "SRR1.2"]
    assert normalize_read_id("@name/2") == ("name", 2)


def test_paired_fastq_rejects_count_and_name_mismatches(tmp_path):
    reads1, reads2 = tmp_path / "r1.fastq", tmp_path / "r2.fastq"
    _write_fastq(reads1, [("a/1", "ACGT"), ("b/1", "ACGT")])
    _write_fastq(reads2, [("a/2", "ACGT")])
    with pytest.raises(ValueError, match="different record counts"):
        list(read_fastq_pairs(reads1, reads2))
    _write_fastq(reads2, [("z/2", "ACGT"), ("b/2", "ACGT")])
    with pytest.raises(ValueError, match="IDs differ"):
        list(read_fastq_pairs(reads1, reads2))


def test_qc_retains_good_mate_as_orphan():
    pair = ReadPair("m", ReadRecord("m/1", "A" * 50, "I" * 50),
                    ReadRecord("m/2", "A" * 10, "I" * 10))
    retained, metrics = qc_read_pairs([pair], 40, 20, 0, 0.5)
    assert len(retained) == 1 and retained[0].read2 is None
    assert metrics["orphan_reads"] == 1


def test_direct_overlap_merge_reconstructs_product():
    product = _product()
    pair = ReadPair("m", ReadRecord("m/1", product[:60], "I" * 60),
                    ReadRecord("m/2", revcomp(product[-60:]), "I" * 60))
    assert merge_read_pair(pair).sequence == product


def test_metadata_alias_normalization_and_myoga_id_consistency():
    metadata = normalize_metadata_row(
        {"sra_run": "SRR123", "lat": "1.5", "lon": "-2.5", "geo_loc_name": "Test"}
    )
    row = myoga_sample_row("SRR123", metadata,
                           {"read_technology": "illumina", "complete_loci": 1},
                           {"sample_id": "SRR123", "L1": "4"}, 1)
    assert row["genome_id"] == row["sample_id"] == "SRR123"
    assert row["latitude"] == "1.5" and row["longitude"] == "-2.5"


def test_validation_distinguishes_exact_interval_and_false_exact_calls():
    truth = [{"sample_id": "s", "locus_id": name, "repeat_count": count}
             for name, count in (("L1", "4"), ("L2", "8"), ("L3", "5"))]
    observed = [
        {"sample_id": "s", "locus_id": "L1", "repeat_count": "4"},
        {"sample_id": "s", "locus_id": "L2", "repeat_count": "",
         "repeat_count_min": "7", "repeat_count_max": "9"},
        {"sample_id": "s", "locus_id": "L3", "repeat_count": "6"},
    ]
    details = compare_call_sets(truth, observed, "illumina")
    assert [row["classification"] for row in details] == [
        "exact_match", "within_reported_interval", "incorrect_call"
    ]
    assert summarize_validation(details)[0]["false_exact_call_rate"] == 0.5


@pytest.mark.skipif(
    shutil.which("bowtie2") is None or shutil.which("bowtie2-build") is None,
    reason="Bowtie2 unavailable",
)
def test_canonical_short_read_pipeline_calls_recoverable_product(tmp_path):
    panel = tmp_path / "panel.tsv"
    reads1, reads2 = tmp_path / "reads_1.fastq", tmp_path / "reads_2.fastq"
    _write_panel(panel)
    product = _product()
    _write_fastq(reads1, [(f"m{i}/1", product[:60]) for i in range(4)])
    _write_fastq(reads2, [(f"m{i}/2", revcomp(product[-60:])) for i in range(4)])
    result = run_short_read_call(str(reads1), str(reads2), str(panel),
                                 str(tmp_path / "out"), "sample", min_depth=1)
    with result["calls"].open() as handle:
        call = next(csv.DictReader(handle, delimiter="\t"))
    assert call["repeat_count"] == "4"
    assert call["mlva_method"] == "Bowtie2 short-read mapping"
    assert not (tmp_path / "out" / "short_read_assembly_summary.tsv").exists()