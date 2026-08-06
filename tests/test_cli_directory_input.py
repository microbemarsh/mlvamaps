import csv
from pathlib import Path

import pytest

import mlvamaps.cli as cli


def _write_panel(path: Path) -> Path:
    panel = path / "primers.tsv"
    panel.write_text(
        "locus_id\tforward_primer\treverse_primer\n"
        "VNTR_01\tACGTTGCAAC\tTGCATGCAAA\n"
    )
    return panel


def test_input_directory_discovers_supported_files_in_stable_order(tmp_path):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "b.FASTA.GZ").write_bytes(b"")
    (inputs / "A.fastq").write_text("")
    (inputs / "notes.txt").write_text("")
    (inputs / "nested").mkdir()
    (inputs / "nested" / "ignored.fasta").write_text("")

    assert [path.name for path in cli._input_files(str(inputs))] == [
        "A.fastq",
        "b.FASTA.GZ",
    ]


def test_call_directory_dispatches_each_file_to_its_sample_outdir(
    tmp_path, monkeypatch
):
    panel = _write_panel(tmp_path)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "reads.fastq.gz").write_bytes(b"")
    (inputs / "assembly.fna").write_text(">contig\nACGT\n")
    (inputs / "README").write_text("ignored")
    observed = []

    def fake_run(args, input_path, outdir, sample_id):
        observed.append((input_path.name, outdir, sample_id))

    monkeypatch.setattr(cli, "_run_single_input", fake_run)

    assert cli.main(
        ["call", "-p", str(panel), "-i", str(inputs), "--outdir", str(tmp_path / "results")]
    ) == 0
    assert observed == [
        ("assembly.fna", tmp_path / "results" / "assembly", "assembly"),
        ("reads.fastq.gz", tmp_path / "results" / "reads", "reads"),
    ]


def test_call_directory_rejects_duplicate_derived_sample_ids(
    tmp_path, monkeypatch
):
    panel = _write_panel(tmp_path)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "sample.fastq").write_text("")
    (inputs / "sample.fasta").write_text("")
    monkeypatch.setattr(cli, "_run_single_input", lambda *args: None)

    with pytest.raises(SystemExit, match="2"):
        cli.main(["call", "-p", str(panel), "-i", str(inputs)])


def test_empty_input_directory_has_a_clear_error(tmp_path, capsys):
    panel = _write_panel(tmp_path)
    inputs = tmp_path / "inputs"
    inputs.mkdir()

    with pytest.raises(SystemExit, match="2"):
        cli.main(["call", "-p", str(panel), "-i", str(inputs)])

    assert "contains no supported FASTA or FASTQ files" in capsys.readouterr().err


def test_combine_legacy_fingerprints_writes_one_mlva_finder_table(tmp_path):
    first = tmp_path / "first.csv"
    second = tmp_path / "second.csv"
    first.write_text(
        "key,Access_number,vrrA,Bams01\n001,GCF_000001,10,16\n"
    )
    second.write_text(
        "key,Access_number,vrrA,Bams01\n001,GCF_000002,10.5,\n"
    )

    output = cli._combine_legacy_fingerprints(
        [first, second], tmp_path / "MLVA_analysis_genomes.csv"
    )

    with output.open(newline="") as handle:
        assert list(csv.reader(handle)) == [
            ["key", "Access_number", "vrrA", "Bams01"],
            ["001", "GCF_000001", "10", "16"],
            ["002", "GCF_000002", "10.5", ""],
        ]


def test_call_directory_writes_combined_mlva_finder_analysis(
    tmp_path, monkeypatch, capsys
):
    panel = _write_panel(tmp_path)
    inputs = tmp_path / "Ba_ref_genomes"
    inputs.mkdir()
    (inputs / "GCF_000002.fna").write_text(">contig\nACGT\n")
    (inputs / "GCF_000001.fna").write_text(">contig\nACGT\n")

    def fake_run(args, input_path, outdir, sample_id):
        outdir.mkdir(parents=True)
        fingerprint = outdir / "legacy_mlva_analysis.csv"
        fingerprint.write_text(
            "key,Access_number,VNTR_01\n001," + sample_id + ",5\n"
        )
        return {"legacy_fingerprint": fingerprint}

    monkeypatch.setattr(cli, "_run_single_input", fake_run)
    results = tmp_path / "results"

    assert cli.main(
        ["call", "-p", str(panel), "-i", str(inputs), "--outdir", str(results)]
    ) == 0
    combined = results / "MLVA_analysis_Ba_ref_genomes.csv"
    with combined.open(newline="") as handle:
        assert list(csv.reader(handle)) == [
            ["key", "Access_number", "VNTR_01"],
            ["001", "GCF_000001", "5"],
            ["002", "GCF_000002", "5"],
        ]
    assert str(combined) in capsys.readouterr().out
