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


def test_short_read_wrapper_forwards_automatic_taxon_options(monkeypatch):
    observed = {}

    def fake_call(**kwargs):
        observed.update(kwargs)
        return {}

    monkeypatch.setattr("mlvamaps.short_read_mapping.run_mapping_short_read_call", fake_call)
    run_short_read_call(
        "r1.fastq", "r2.fastq", "panel.tsv", "out", "sample",
        database_path="multi-taxon-db", taxon_identification=True, taxon_min_loci=5,
        classification_repeat_scale=0.2,
    )

    assert observed["database_path"] == "multi-taxon-db"
    assert observed["taxon_identification"] is True
    assert observed["taxon_min_loci"] == 5
    assert observed["classification_repeat_scale"] == 0.2


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
    assert call["mlva_method"] == "competitive_sample_likelihood"
    assert not (tmp_path / "out" / "short_read_assembly_summary.tsv").exists()


def test_legacy_primer_cli_requires_database(tmp_path, monkeypatch, capsys):
    from mlvamaps.cli import main
    def forbidden(*args, **kwargs):
        pytest.fail('Database sequences requested by a primer-only call')
    monkeypatch.setattr('mlvamaps.candidate_contexts._base_contexts', forbidden)
    panel = tmp_path/'primers.tsv'
    locus = _locus()
    product = _product()
    panel.write_text(f'L1_4bp_{len(product)}bp_4U\t{locus.forward_primer}\t{locus.reverse_primer}\n')
    first, second = tmp_path/'r1.fq', tmp_path/'r2.fq'
    _write_fastq(first, [(f'm{i}/1', product[:60]) for i in range(4)])
    _write_fastq(second, [(f'm{i}/2', revcomp(product[-60:])) for i in range(4)])
    out = tmp_path/'out'
    with pytest.raises(SystemExit) as exc:
        main(['call', '-p', str(panel), '-i', 'sr', '--fq1', str(first), '--fq2', str(second),
              '-o', str(out), '--sample-id', 'legacy', '-t', '1'])
    assert exc.value.code == 2
    assert 'requires --database' in capsys.readouterr().err
    assert not out.exists()
