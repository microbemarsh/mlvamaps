from __future__ import annotations

import csv
from pathlib import Path

import pytest

from mlvamaps.io import normalize_read_id, read_fastq_pairs
from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.sample_metadata import myoga_sample_row, normalize_metadata_row
from mlvamaps.short_reads import (
    LocusAssembly,
    RecruitedPair,
    _assemble_one_locus,
    _call_locus,
    _skesa_command,
    check_skesa,
    estimate_insert_size_distribution,
    merge_read_pair,
    qc_read_pairs,
    recruit_read_pairs,
)
from mlvamaps.validation import compare_call_sets, summarize_validation


def _write_fastq(path: Path, records: list[tuple[str, str]], quality: str = "I") -> None:
    path.write_text(
        "".join(
            f"@{name}\n{sequence}\n+\n{quality * len(sequence)}\n"
            for name, sequence in records
        )
    )


def _write_fake_skesa(
    path: Path, contig: str | None = None, assembly_exit_code: int = 0
) -> Path:
    output_text = "" if contig is None else f">Contig_1\n{contig}\n"
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import pathlib, sys\n"
        "if '--version' in sys.argv:\n"
        "    print('SKESA 2.test')\n"
        "    raise SystemExit(0)\n"
        "output = pathlib.Path(sys.argv[sys.argv.index('--contigs_out') + 1])\n"
        f"output.write_text({output_text!r})\n"
        f"raise SystemExit({assembly_exit_code})\n"
    )
    path.chmod(0o755)
    return path


def _locus() -> Locus:
    return Locus(
        "L1",
        forward_primer="ACGTTGCAACGTTGCAACGT",
        reverse_primer="AGTCAGTCAGTCAGTCAGTC",
        left_flank_sequence="CCGGAATTCGACCTGA",
        right_flank_sequence="TTAACCGGCTAGGTCA",
        repeat_motif="ATGC",
        repeat_unit_length_bp=4,
        expected_min_repeats=2,
        expected_max_repeats=10,
        nominal_repeat_units=4,
    )


def _product(repeats: int = 4) -> str:
    locus = _locus()
    return (
        locus.forward_primer
        + locus.left_flank_sequence
        + locus.repeat_motif * repeats
        + locus.right_flank_sequence
        + "GACTGACTGACTGACTGACT"
    )


def test_read_id_normalization_and_paired_fastq_preserve_mates(tmp_path):
    reads1 = tmp_path / "r1.fastq"
    reads2 = tmp_path / "r2.fastq"
    _write_fastq(reads1, [("SRR1.1/1", "ACGT"), ("SRR1.2/1", "TGCA")])
    _write_fastq(reads2, [("SRR1.1/2", "TGCA"), ("SRR1.2/2", "ACGT")])
    pairs = list(read_fastq_pairs(reads1, reads2))
    assert [pair.molecule_id for pair in pairs] == ["SRR1.1", "SRR1.2"]
    assert pairs[0].read1.read_id.endswith("/1")
    assert pairs[0].read2 is not None and pairs[0].read2.read_id.endswith("/2")
    assert normalize_read_id("@name/2") == ("name", 2)


def test_paired_fastq_rejects_count_and_name_mismatches(tmp_path):
    reads1 = tmp_path / "r1.fastq"
    reads2 = tmp_path / "r2.fastq"
    _write_fastq(reads1, [("a/1", "ACGT"), ("b/1", "ACGT")])
    _write_fastq(reads2, [("a/2", "ACGT")])
    with pytest.raises(ValueError, match="different record counts"):
        list(read_fastq_pairs(reads1, reads2))
    _write_fastq(reads2, [("z/2", "ACGT"), ("b/2", "ACGT")])
    with pytest.raises(ValueError, match="IDs differ"):
        list(read_fastq_pairs(reads1, reads2))


def test_qc_retains_pairing_and_reports_orphan():
    pair = ReadPair(
        "m1",
        ReadRecord("m1/1", "A" * 50, "I" * 50),
        ReadRecord("m1/2", "A" * 10, "I" * 10),
    )
    retained, metrics = qc_read_pairs([pair], 40, 20, 0, 0.5)
    assert len(retained) == 1 and retained[0].read2 is None
    assert metrics["orphan_reads"] == 1


def test_competitive_pair_recruitment_has_unique_ambiguous_and_discordant_outcomes():
    index = {"L1": {"AAAAAAAAAAAAAAA"}, "L2": {"CCCCCCCCCCCCCCC"}}
    pairs = [
        ReadPair("unique", ReadRecord("unique/1", "A" * 30, "I" * 30)),
        ReadPair("ambiguous", ReadRecord("ambiguous/1", "A" * 15 + "C" * 15, "I" * 30)),
        ReadPair(
            "discordant",
            ReadRecord("discordant/1", "A" * 30, "I" * 30),
            ReadRecord("discordant/2", "C" * 30, "I" * 30),
        ),
    ]
    outcomes = recruit_read_pairs(pairs, index, k=15, min_seeds=1, min_margin=1)
    assert outcomes[0].outcome == "unique" and outcomes[0].locus_id == "L1"
    assert outcomes[1].outcome == "ambiguous"
    assert outcomes[2].outcome == "discordant"


def test_strict_overlap_merge_and_skesa_command(tmp_path):
    product = _product(4)
    read1 = product[:60]
    read2 = product[-60:]
    pair = ReadPair(
        "m1",
        ReadRecord("m1/1", read1, "I" * len(read1)),
        ReadRecord("m1/2", __import__("mlvamaps.sequence", fromlist=["revcomp"]).revcomp(read2), "I" * len(read2)),
    )
    merged = merge_read_pair(pair)
    assert merged is not None and merged.sequence == product
    executable = _write_fake_skesa(tmp_path / "skesa")
    assert check_skesa(str(executable)) == str(executable)
    command = _skesa_command(
        str(executable),
        (tmp_path / "r1.fastq", tmp_path / "r2.fastq"),
        None,
        tmp_path / "contigs.fasta",
        3,
    )
    assert command[-2:] == ["--reads", f"{tmp_path / 'r1.fastq'},{tmp_path / 'r2.fastq'}"]
    assert command[command.index("--cores") + 1] == "3"


def test_missing_skesa_is_an_explicit_error():
    with pytest.raises(RuntimeError, match="SKESA executable"):
        check_skesa("definitely-not-an-installed-skesa")


def test_skesa_contigs_are_used_and_failure_has_no_fallback(tmp_path):
    locus = _locus()
    product = _product(4)
    recruited = [
        RecruitedPair(
            ReadPair("m", ReadRecord("m/1", product, "I" * len(product))),
            "unique",
            "L1",
            10,
            1.0,
        )
    ]
    executable = _write_fake_skesa(tmp_path / "skesa", contig=product)
    assembly = _assemble_one_locus(
        locus,
        recruited,
        product,
        str(executable),
        tmp_path / "work",
        2,
    )
    assert assembly.status == "ASSEMBLED"
    assert assembly.contigs == (product,)
    assert "--cores 2" in (tmp_path / "work" / "skesa.log").read_text()

    failing = _write_fake_skesa(
        tmp_path / "failing-skesa", assembly_exit_code=7
    )
    with pytest.raises(RuntimeError, match="SKESA failed for locus 'L1'"):
        _assemble_one_locus(
            locus,
            recruited,
            product,
            str(failing),
            tmp_path / "failed-work",
            1,
        )


def test_exact_call_requires_direct_boundary_evidence_and_partial_stays_unresolved():
    locus = _locus()
    product = _product(4)
    exact_pair = ReadPair("exact", ReadRecord("exact/1", product, "I" * len(product)))
    recruited = [RecruitedPair(exact_pair, "unique", "L1", 10, 1.0)]
    exact = _call_locus("sample", locus, recruited, LocusAssembly("L1", (), (), 2.0, "NO_CONTIGS"), 1)
    assert exact["repeat_count"] == 4
    assert exact["evidence_class"] == "BOUNDARY_SPANNING_SINGLE_READ"

    partial_sequence = locus.forward_primer + locus.left_flank_sequence + locus.repeat_motif * 2
    partial_pair = ReadPair("partial", ReadRecord("partial/1", partial_sequence, "I" * len(partial_sequence)))
    partial = _call_locus(
        "sample",
        locus,
        [RecruitedPair(partial_pair, "unique", "L1", 10, 1.0)],
        LocusAssembly("L1", (), (), 1.0, "NO_CONTIGS"),
        1,
    )
    assert partial["repeat_count"] == ""
    assert partial["evidence_class"] in {"PARTIAL_REPEAT_EVIDENCE", "PRESENCE_ONLY"}


def test_mixed_alleles_are_retained_with_informative_fractions():
    locus = _locus()
    recruited = [
        RecruitedPair(ReadPair("a", ReadRecord("a/1", _product(4), "I" * len(_product(4)))), "unique", "L1", 10, 1.0),
        RecruitedPair(ReadPair("b", ReadRecord("b/1", _product(5), "I" * len(_product(5)))), "unique", "L1", 10, 1.0),
    ]
    row = _call_locus("sample", locus, recruited, LocusAssembly("L1", (), (), 2.0, "NO_CONTIGS"), 1)
    assert row["evidence_class"] == "MULTIPLE_ALLELES"
    assert row["num_confirmed_secondary_variants"] == 1
    assert row["estimated_primary_fraction"] == 0.5


def test_insert_size_estimation_requires_opposite_orientations():
    from mlvamaps.sequence import revcomp

    reference = _product(4)
    pair = ReadPair(
        "m",
        ReadRecord("m/1", reference[:35], "I" * 35),
        ReadRecord("m/2", revcomp(reference[-35:]), "I" * 35),
    )
    estimate = estimate_insert_size_distribution(
        [RecruitedPair(pair, "unique", "L1", 10, 1.0)], {"L1": reference}
    )
    assert estimate == (float(len(reference)), 0.0, 1)


def test_metadata_alias_normalization_and_myoga_id_consistency():
    metadata = normalize_metadata_row(
        {"sra_run": "SRR123", "lat": "1.5", "lon": "-2.5", "geo_loc_name": "Test"}
    )
    assert metadata["sample_id"] == "SRR123"
    row = myoga_sample_row(
        "SRR123",
        metadata,
        {"read_technology": "illumina", "complete_loci": 1},
        {"sample_id": "SRR123", "L1": "4"},
        1,
    )
    assert row["genome_id"] == row["sample_id"] == "SRR123"
    assert row["latitude"] == "1.5" and row["longitude"] == "-2.5"


def test_validation_distinguishes_exact_interval_and_false_exact_calls():
    truth = [
        {"sample_id": "s", "locus_id": "L1", "repeat_count": "4"},
        {"sample_id": "s", "locus_id": "L2", "repeat_count": "8"},
        {"sample_id": "s", "locus_id": "L3", "repeat_count": "5"},
    ]
    observed = [
        {"sample_id": "s", "locus_id": "L1", "repeat_count": "4"},
        {"sample_id": "s", "locus_id": "L2", "repeat_count": "", "repeat_count_min": "7", "repeat_count_max": "9"},
        {"sample_id": "s", "locus_id": "L3", "repeat_count": "6"},
    ]
    details = compare_call_sets(truth, observed, "illumina")
    assert [row["classification"] for row in details] == [
        "exact_match",
        "within_reported_interval",
        "incorrect_call",
    ]
    summary = summarize_validation(details)[0]
    assert summary["false_exact_call_rate"] == 0.5


def _write_panel(path: Path) -> None:
    locus = _locus()
    path.write_text(
        "locus_id\tforward_primer\treverse_primer\tleft_flank_sequence\tright_flank_sequence\trepeat_motif\trepeat_unit_length_bp\texpected_min_repeats\texpected_max_repeats\tnominal_repeat_units\n"
        f"{locus.locus_id}\t{locus.forward_primer}\t{locus.reverse_primer}\t{locus.left_flank_sequence}\t{locus.right_flank_sequence}\t{locus.repeat_motif}\t4\t2\t10\t4\n"
    )


def test_short_read_pipeline_agrees_with_assembly_truth_when_product_is_recoverable(tmp_path):
    from mlvamaps.assembly_call import run_assembly_call
    from mlvamaps.sequence import revcomp
    from mlvamaps.short_reads import run_short_read_call

    panel = tmp_path / "panel.tsv"
    assembly = tmp_path / "assembly.fasta"
    reads1 = tmp_path / "reads_1.fastq"
    reads2 = tmp_path / "reads_2.fastq"
    metadata = {"sample_id": "SRR_TEST", "latitude": "10", "longitude": "20"}
    product = _product(4)
    _write_panel(panel)
    assembly.write_text(f">contig\nTTTT{product}CCCC\n")
    _write_fastq(reads1, [(f"m{index}/1", product[:60]) for index in range(4)])
    _write_fastq(reads2, [(f"m{index}/2", revcomp(product[-60:])) for index in range(4)])
    assembly_result = run_assembly_call(
        str(assembly),
        str(panel),
        str(tmp_path / "assembly_result"),
        "SRR_TEST",
        algorithm="legacy",
    )
    short_result = run_short_read_call(
        str(reads1),
        str(reads2),
        str(panel),
        str(tmp_path / "short_result"),
        "SRR_TEST",
        sample_metadata=metadata,
        min_depth=1,
        skesa_bin=str(_write_fake_skesa(tmp_path / "skesa")),
    )
    with assembly_result["calls"].open() as handle:
        assembly_call = next(csv.DictReader(handle, delimiter="\t"))
    with short_result["calls"].open() as handle:
        short_call = next(csv.DictReader(handle, delimiter="\t"))
    assert short_call["repeat_count"] == assembly_call["repeat_count"] == "4"
    assert short_call["evidence_class"] in {
        "COMPLETE_ASSEMBLED_PRODUCT",
        "BOUNDARY_SPANNING_READ_PAIR",
    }
    with short_result["myoga_samples"].open() as handle:
        myoga = next(csv.DictReader(handle))
    assert myoga["genome_id"] == "SRR_TEST"
    assert myoga["latitude"] == "10"


def test_manifest_isolates_failure_and_writes_combined_tables(tmp_path):
    from mlvamaps.cli import main
    from mlvamaps.sequence import revcomp

    panel = tmp_path / "panel.tsv"
    reads1 = tmp_path / "good_1.fastq"
    reads2 = tmp_path / "good_2.fastq"
    manifest = tmp_path / "manifest.tsv"
    product = _product(4)
    _write_panel(panel)
    _write_fastq(reads1, [("m/1", product[:60])])
    _write_fastq(reads2, [("m/2", revcomp(product[-60:]))])
    manifest.write_text(
        "sample_id\treads1\treads2\n"
        f"SRR_GOOD\t{reads1}\t{reads2}\n"
        f"SRR_BAD\t{tmp_path / 'missing.fastq'}\t.\n"
    )
    outdir = tmp_path / "batch"
    assert main(
        [
            "call",
            "-p",
            str(panel),
            "-i",
            "sr",
            "--manifest",
            str(manifest),
            "--read-technology",
            "illumina",
            "--min-depth",
            "1",
            "--skesa-bin",
            str(_write_fake_skesa(tmp_path / "skesa")),
            "-o",
            str(outdir),
        ]
    ) == 0
    with (outdir / "batch_status.tsv").open() as handle:
        status = {row["sample_id"]: row["status"] for row in csv.DictReader(handle, delimiter="\t")}
    assert status == {"SRR_GOOD": "success", "SRR_BAD": "failed"}
    assert (outdir / "calls.tsv").exists()
    assert (outdir / "myoga_samples.csv").exists()
    assert not [
        path
        for path in outdir.rglob("*")
        if path.is_file() and path.suffix.lower() in {".fasta", ".fa", ".fastq", ".fq"}
    ]


def _read_first_call(path: Path) -> dict[str, str]:
    with path.open() as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def test_nonoverlapping_opposite_boundaries_report_interval_not_midpoint(tmp_path):
    from mlvamaps.sequence import revcomp
    from mlvamaps.short_reads import run_short_read_call

    locus = _locus()
    locus = Locus(**{**locus.__dict__, "expected_max_repeats": 40, "nominal_repeat_units": 25})
    product = (
        locus.forward_primer
        + locus.left_flank_sequence
        + locus.repeat_motif * 25
        + locus.right_flank_sequence
        + revcomp(locus.reverse_primer)
    )
    panel = tmp_path / "panel.tsv"
    panel.write_text(
        "locus_id\tforward_primer\treverse_primer\tleft_flank_sequence\tright_flank_sequence\trepeat_motif\trepeat_unit_length_bp\texpected_min_repeats\texpected_max_repeats\tnominal_repeat_units\n"
        f"L1\t{locus.forward_primer}\t{locus.reverse_primer}\t{locus.left_flank_sequence}\t{locus.right_flank_sequence}\tATGC\t4\t2\t40\t25\n"
    )
    reads1, reads2 = tmp_path / "r1.fastq", tmp_path / "r2.fastq"
    _write_fastq(reads1, [("m/1", product[:60])])
    _write_fastq(reads2, [("m/2", revcomp(product[-60:]))])
    result = run_short_read_call(
        str(reads1),
        str(reads2),
        str(panel),
        str(tmp_path / "out"),
        "s",
        min_depth=1,
        skesa_bin=str(_write_fake_skesa(tmp_path / "skesa")),
    )
    row = _read_first_call(result["calls"])
    assert row["repeat_count"] == ""
    assert row["repeat_count_min"] == "2" and row["repeat_count_max"] == "40"
    assert row["evidence_class"] == "BOUNDARY_SPANNING_READ_PAIR"


def test_one_boundary_and_low_depth_are_reported_honestly(tmp_path):
    from mlvamaps.sequence import revcomp
    from mlvamaps.short_reads import run_short_read_call

    panel = tmp_path / "panel.tsv"
    _write_panel(panel)
    long_product = _product(20)
    partial = tmp_path / "partial.fastq"
    _write_fastq(partial, [("p/1", long_product[:60])])
    skesa = _write_fake_skesa(tmp_path / "skesa")
    partial_result = run_short_read_call(
        str(partial),
        None,
        str(panel),
        str(tmp_path / "partial_out"),
        "partial",
        min_depth=1,
        skesa_bin=str(skesa),
    )
    partial_row = _read_first_call(partial_result["calls"])
    assert partial_row["repeat_count"] == ""
    assert partial_row["evidence_class"] in {"PARTIAL_REPEAT_EVIDENCE", "PRESENCE_ONLY"}

    product = _product(4)
    reads1, reads2 = tmp_path / "low_1.fastq", tmp_path / "low_2.fastq"
    _write_fastq(reads1, [("low/1", product[:60])])
    _write_fastq(reads2, [("low/2", revcomp(product[-60:]))])
    low_result = run_short_read_call(
        str(reads1),
        str(reads2),
        str(panel),
        str(tmp_path / "low_out"),
        "low",
        min_depth=3,
        skesa_bin=str(skesa),
    )
    low_row = _read_first_call(low_result["calls"])
    assert low_row["repeat_count"] == "4"
    assert low_row["status"] == "LOW_DEPTH"
    assert low_row["evidence_class"] == "LOW_DEPTH"


def test_pipeline_preserves_two_directly_observed_alleles(tmp_path):
    from mlvamaps.short_reads import run_short_read_call

    panel = tmp_path / "panel.tsv"
    reads = tmp_path / "mixed.fastq"
    _write_panel(panel)
    records = [
        ("a1/1", _product(4)),
        ("a2/1", _product(4)),
        ("b1/1", _product(5)),
    ]
    _write_fastq(reads, records)
    result = run_short_read_call(
        str(reads),
        None,
        str(panel),
        str(tmp_path / "out"),
        "mixed",
        min_depth=1,
        skesa_bin=str(_write_fake_skesa(tmp_path / "skesa")),
    )
    row = _read_first_call(result["calls"])
    assert row["evidence_class"] == "MULTIPLE_ALLELES"
    assert row["repeat_count"] == "4"
    assert row["secondary_alleles"].startswith("5:")


def test_shared_repeat_only_reads_do_not_cross_recruit_similar_loci():
    sequence = "ATGC" * 20
    index = {
        "L1": {"ACGTTGCAACGTTGC"},
        "L2": {"TTGCCATGTTGCCAT"},
    }
    pair = ReadPair("repeat", ReadRecord("repeat/1", sequence, "I" * len(sequence)))
    result = recruit_read_pairs([pair], index, k=15, min_seeds=1)
    assert result[0].outcome == "unassigned"
