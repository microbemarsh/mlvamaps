from __future__ import annotations

import pytest
import shutil
import csv
import hashlib
from pathlib import Path

from mlvamaps.models import Locus
from mlvamaps.short_read_mapping import (
    LocusContext,
    candidate_repeat_counts,
    estimate_insert_size_distribution,
    expand_candidate_contexts,
    bowtie2_mapping_command,
    build_context_index,
    load_locus_contexts,
)
from mlvamaps.validation import fastq_assembly_concordance


def _locus(unit: int = 4) -> Locus:
    return Locus(
        "L1", forward_primer="ACGT", reverse_primer="TGCA",
        left_flank_sequence="AACCGGTT", right_flank_sequence="GGTTAACC",
        repeat_motif="ATGC" if unit == 4 else "A" * unit,
        repeat_unit_length_bp=unit, expected_min_repeats=2,
        expected_max_repeats=8, nominal_repeat_units=4,
    )


def test_candidate_counts_are_discrete_and_bounded():
    assert candidate_repeat_counts(_locus(), [4, 12], maximum=6) == [2, 3, 4, 5, 6]
    with pytest.raises(ValueError, match="at least 1"):
        candidate_repeat_counts(_locus(), [4], maximum=0)


def test_candidate_context_expansion_changes_only_whole_repeat_units():
    locus = _locus()
    sequence = "ACGTAACCGGTT" + "ATGC" * 4 + "GGTTAACCTGCA"
    context = LocusContext(
        "ctx", "L1", sequence, "REF", expected_repeat_count=4,
        repeat_motif="ATGC", repeat_unit_length_bp=4,
        repeat_start=12, repeat_end=28,
    )
    expanded = expand_candidate_contexts([context], [locus], maximum=6)
    assert [item.expected_repeat_count for item in expanded] == [2, 3, 4, 5, 6]
    assert [len(item.sequence) for item in expanded] == [32, 36, 40, 44, 48]


def test_insert_estimation_is_robust_to_outlier_and_reports_low_support():
    estimate = estimate_insert_size_distribution([295, 300, 305, 310, 5000])
    assert estimate.pairs_used == 4
    assert estimate.median == 302.5
    assert estimate.mean == 302.5
    assert estimate.mad == 5.0
    sparse = estimate_insert_size_distribution([300, 301])
    assert sparse.median is None and sparse.pairs_used == 2


def test_fastq_assembly_concordance_counts_exact_unresolved_and_discordant():
    fastq = [
        {"sample_id": "S", "locus_id": "A", "repeat_count": "4", "product_size_bp": "100", "allele_confidence": "0.99"},
        {"sample_id": "S", "locus_id": "B", "repeat_count": "", "failure_reason": "low coverage"},
        {"sample_id": "S", "locus_id": "C", "repeat_count": "6"},
    ]
    assembly = [
        {"sample_id": "S", "locus_id": "A", "repeat_count": "4", "product_size_bp": "100"},
        {"sample_id": "S", "locus_id": "B", "repeat_count": "5"},
        {"sample_id": "S", "locus_id": "C", "repeat_count": "7"},
    ]
    details, summary = fastq_assembly_concordance(fastq, assembly)
    assert [row["agreement"] for row in details] == ["exact", "fastq_unresolved", "discordant"]
    assert summary == {
        "comparable_loci": 2, "exact_repeat_count_agreement": 1,
        "exact_length_agreement": 1, "repeat_count_concordance": 0.5,
        "fastq_only_calls": 0, "assembly_only_calls": 1,
        "fastq_unresolved": 1, "discordant_calls": 1,
    }


def test_bowtie2_command_is_single_competitive_paired_mapping(tmp_path):
    command = bowtie2_mapping_command(
        "/opt/bowtie2", tmp_path / "index", tmp_path / "r1.fq",
        tmp_path / "r2.fq", tmp_path / "out.sam", 8,
    )
    assert command.count("-x") == 1
    assert command[command.index("-p") + 1] == "8"
    assert command[command.index("-1") + 1].endswith("r1.fq")
    assert command[command.index("-2") + 1].endswith("r2.fq")
    assert "-a" in command


@pytest.mark.skipif(shutil.which("bowtie2-build") is None, reason="bowtie2-build unavailable")
def test_context_index_build_is_versioned_and_reusable(tmp_path):
    locus = _locus()
    sequence = "ACGTAACCGGTT" + "ATGC" * 4 + "GGTTAACCTGCA"
    contexts = expand_candidate_contexts([
        LocusContext("ctx", "L1", sequence, "REF", expected_repeat_count=4,
                     repeat_motif="ATGC", repeat_unit_length_bp=4,
                     repeat_start=12, repeat_end=28)
    ], [locus], maximum=6)
    first = build_context_index(contexts, tmp_path)
    timestamp = first["metadata"].stat().st_mtime_ns
    second = build_context_index(contexts, tmp_path)
    assert second["prefix"] == first["prefix"]
    assert second["metadata"].stat().st_mtime_ns == timestamp
    assert len(list(Path(tmp_path).glob("mlva_contexts.*.bt2*"))) == 6


def test_old_database_without_context_schema_requires_rebuild(tmp_path):
    (tmp_path / "L1.fasta").write_text(">ref\nACGTACGT\n")
    with pytest.raises(ValueError, match="predates the Illumina context schema"):
        load_locus_contexts([_locus()], tmp_path)


def test_context_database_validates_schema_and_sequence_hash(tmp_path):
    sequence = "ACGTAACCGGTT" + "ATGC" * 4 + "GGTTAACCTGCA"
    (tmp_path / "mlva_contexts.fasta.gz").write_bytes(b"")
    # Use the project writer so gzip handling matches production behavior.
    from mlvamaps.io import write_fasta
    write_fasta([("ctx", sequence)], tmp_path / "mlva_contexts.fasta.gz")
    row = LocusContext(
        "ctx", "L1", sequence, "REF", expected_repeat_count=4,
        repeat_motif="ATGC", repeat_unit_length_bp=4,
        repeat_start=12, repeat_end=28,
    ).row()
    row["schema_version"] = "1.0"
    with (tmp_path / "mlva_contexts.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["schema_version", *[key for key in row if key != "schema_version"]], delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    assert load_locus_contexts([_locus()], tmp_path)[0].sequence == sequence
    text = (tmp_path / "mlva_contexts.tsv").read_text()
    (tmp_path / "mlva_contexts.tsv").write_text(text.replace(hashlib.sha256(sequence.encode()).hexdigest(), "0" * 64))
    with pytest.raises(ValueError, match="hash mismatch"):
        load_locus_contexts([_locus()], tmp_path)