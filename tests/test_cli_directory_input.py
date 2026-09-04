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


def _fake_short_result(outdir: Path, sample_id: str) -> dict[str, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    files = {
        "calls": ("calls.tsv", f"sample_id\tlocus_id\trepeat_count\n{sample_id}\tL1\t4\n"),
        "repeat_counts": ("locus_repeat_counts.tsv", f"sample_id\tlocus_id\trepeat_count\n{sample_id}\tL1\t4\n"),
        "fingerprint": ("mlva_fingerprint.tsv", f"sample_id\tL1\n{sample_id}\t4\n"),
        "sample_summary": ("sample_summary.tsv", f"sample_id\trun_status\n{sample_id}\tsuccess\n"),
        "myoga_samples": ("myoga_samples.csv", f"genome_id,sample_id\n{sample_id},{sample_id}\n"),
        "myoga_loci": ("myoga_loci.csv", f"genome_id,sample_id,locus_id,repeat_count\n{sample_id},{sample_id},L1,4\n"),
    }
    result = {}
    for key, (filename, contents) in files.items():
        result[key] = outdir / filename
        result[key].write_text(contents)
    return result


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


def test_short_read_directory_discovers_exact_pairs_in_prefix_order(tmp_path):
    inputs = tmp_path / "short_reads"
    inputs.mkdir()
    for name in (
        "b_sample_2.fastq.gz",
        "A_sample_1.fastq.gz",
        "b_sample_1.fastq.gz",
        "A_sample_2.fastq.gz",
    ):
        (inputs / name).write_bytes(b"")
    (inputs / "ignored_R1.fastq.gz").write_bytes(b"")
    (inputs / "notes.txt").write_text("ignored")

    rows = cli._short_read_directory_rows(inputs)

    assert [row["sample_id"] for row in rows] == ["A_sample", "b_sample"]
    assert Path(rows[0]["reads1"]).name == "A_sample_1.fastq.gz"
    assert Path(rows[0]["reads2"]).name == "A_sample_2.fastq.gz"


def test_short_read_directory_rejects_orphaned_exact_mates(tmp_path):
    inputs = tmp_path / "short_reads"
    inputs.mkdir()
    (inputs / "sample_1.fastq.gz").write_bytes(b"")

    with pytest.raises(ValueError, match=r"sample \(missing mate 2\)"):
        cli._short_read_directory_rows(inputs)


def test_short_read_directory_cli_routes_pairs_to_batch(tmp_path, monkeypatch):
    panel = _write_panel(tmp_path)
    inputs = tmp_path / "short_reads"
    inputs.mkdir()
    (inputs / "sample_1.fastq.gz").write_bytes(b"")
    (inputs / "sample_2.fastq.gz").write_bytes(b"")
    observed = {}

    def fake_batch(args, parser, rows):
        observed["rows"] = rows
        observed["short_read_mode"] = args.short_read_mode

    monkeypatch.setattr(cli, "_run_short_batch", fake_batch)

    assert cli.main(
        ["call", "-p", str(panel), "-i", str(inputs), "--short-reads"]
    ) == 0
    assert observed["short_read_mode"] is True
    assert [row["sample_id"] for row in observed["rows"]] == ["sample"]


def test_batch_thread_allocation_is_bounded(monkeypatch):
    monkeypatch.setenv("MLVAMAPS_MAX_CONCURRENT_SAMPLES", "3")
    workers, per_sample = cli._batch_allocation(10, 20)
    assert workers == 3
    assert workers * per_sample <= 10


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
    combined = results / "batch_summary" / "MLVA_analysis_Ba_ref_genomes.csv"
    with combined.open(newline="") as handle:
        assert list(csv.reader(handle)) == [
            ["key", "Access_number", "VNTR_01"],
            ["001", "GCF_000001", "5"],
            ["002", "GCF_000002", "5"],
        ]
    assert str(combined) in capsys.readouterr().out
    assert not (results / "MLVA_analysis_Ba_ref_genomes.csv").exists()


@pytest.mark.parametrize(
    ("filename", "contents"),
    [("A.fasta", ">a\nACGT\n"), ("A.fastq", "@a\nACGT\n+\nIIII\n")],
)
def test_normal_directory_batches_keep_samples_separate_from_batch_outputs(
    tmp_path, monkeypatch, filename, contents
):
    panel = _write_panel(tmp_path)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / filename).write_text(contents)
    (inputs / filename.replace("A.", "B.")).write_text(contents.replace("a", "b"))

    def fake_run(args, input_path, outdir, sample_id):
        outdir.mkdir(parents=True)
        (outdir / "calls.tsv").write_text(
            f"sample_id\tlocus_id\trepeat_count\n{sample_id}\tL1\t4\n"
        )
        return {}

    monkeypatch.setattr(cli, "_run_single_input", fake_run)
    results = tmp_path / "results"

    assert cli.main(
        ["call", "-p", str(panel), "-i", str(inputs), "-o", str(results)]
    ) == 0
    assert {path.name for path in results.iterdir()} == {"A", "B", "batch_summary"}
    assert (results / "A" / "calls.tsv").exists()
    assert (results / "B" / "calls.tsv").exists()
    assert not (results / "calls.tsv").exists()


def test_single_input_keeps_user_selected_output_directory(tmp_path, monkeypatch):
    panel = _write_panel(tmp_path)
    assembly = tmp_path / "sample.fasta"
    assembly.write_text(">sample\nACGT\n")
    observed = {}

    def fake_run(args, input_path, outdir, sample_id):
        observed.update(outdir=outdir, sample_id=sample_id)
        outdir.mkdir(parents=True)
        (outdir / "calls.tsv").write_text("sample_id\n")
        return {}

    monkeypatch.setattr(cli, "_run_single_input", fake_run)
    results = tmp_path / "sample_results"

    assert cli.main(
        ["call", "-p", str(panel), "-i", str(assembly), "-o", str(results)]
    ) == 0
    assert observed == {"outdir": results, "sample_id": "sample"}
    assert (results / "calls.tsv").exists()
    assert not (results / "sample").exists()
    assert not (results / "batch_summary").exists()


def test_directory_batch_rejects_reserved_sample_id(tmp_path, capsys):
    panel = _write_panel(tmp_path)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "batch_summary.fasta").write_text(">sample\nACGT\n")

    with pytest.raises(SystemExit, match="2"):
        cli.main(["call", "-p", str(panel), "-i", str(inputs)])

    assert "reserved for batch aggregate outputs" in capsys.readouterr().err


def test_paired_fastq_directory_writes_clean_batch_layout(tmp_path, monkeypatch):
    panel = _write_panel(tmp_path)
    inputs = tmp_path / "reads"
    inputs.mkdir()
    for sample_id in ("sample1", "sample2"):
        (inputs / f"{sample_id}_1.fastq.gz").write_bytes(b"")
        (inputs / f"{sample_id}_2.fastq.gz").write_bytes(b"")

    def fake_run(args, reads1, reads2, outdir, sample_id, metadata):
        return _fake_short_result(outdir, sample_id)

    monkeypatch.setattr(cli, "_run_short_input", fake_run)
    results = tmp_path / "results"

    assert cli.main(
        ["call", "-p", str(panel), "-i", str(inputs), "--short-reads", "-o", str(results)]
    ) == 0
    assert {path.name for path in results.iterdir()} == {
        "sample1", "sample2", "batch_summary"
    }
    for sample_id in ("sample1", "sample2"):
        assert (results / sample_id / "calls.tsv").exists()
        assert (results / sample_id / "myoga_samples.csv").exists()
    for filename in (
        "batch_status.tsv", "calls.tsv", "locus_repeat_counts.tsv",
        "mlva_fingerprint.tsv", "sample_summary.tsv", "myoga_samples.csv",
        "myoga_loci.csv",
    ):
        assert (results / "batch_summary" / filename).exists()
        assert not (results / filename).exists()


def test_manifest_batch_layout_resumes_from_sample_directory(tmp_path, monkeypatch):
    panel = _write_panel(tmp_path)
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "sample_id\treads1\n"
        f"sample1\t{reads}\n"
        f"sample2\t{reads}\n"
    )
    calls = []

    def fake_run(args, reads1, reads2, outdir, sample_id, metadata):
        calls.append(sample_id)
        return _fake_short_result(outdir, sample_id)

    monkeypatch.setattr(cli, "_run_short_input", fake_run)
    results = tmp_path / "results"
    command = [
        "call", "-p", str(panel), "-i", "sr", "--manifest", str(manifest),
        "-o", str(results),
    ]

    assert cli.main(command) == 0
    assert cli.main(command) == 0
    assert calls == ["sample1", "sample2"]
    with (results / "batch_summary" / "batch_status.tsv").open() as handle:
        statuses = [row["status"] for row in csv.DictReader(handle, delimiter="\t")]
    assert statuses == ["skipped_success", "skipped_success"]
    assert {path.name for path in results.iterdir()} == {
        "sample1", "sample2", "batch_summary"
    }


def test_manifest_rejects_reserved_batch_summary_sample_id(tmp_path, capsys):
    panel = _write_panel(tmp_path)
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(f"sample_id\treads1\nbatch_summary\t{reads}\n")

    with pytest.raises(SystemExit, match="2"):
        cli.main(
            ["call", "-p", str(panel), "-i", "sr", "--manifest", str(manifest)]
        )

    assert "reserved for batch aggregate outputs" in capsys.readouterr().err
