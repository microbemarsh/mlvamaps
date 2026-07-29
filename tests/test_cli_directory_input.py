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
        ["call", str(panel), str(inputs), "--outdir", str(tmp_path / "results")]
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
        cli.main(["call", str(panel), str(inputs)])


def test_empty_input_directory_has_a_clear_error(tmp_path, capsys):
    panel = _write_panel(tmp_path)
    inputs = tmp_path / "inputs"
    inputs.mkdir()

    with pytest.raises(SystemExit, match="2"):
        cli.main(["call", str(panel), str(inputs)])

    assert "contains no supported FASTA or FASTQ files" in capsys.readouterr().err
