from __future__ import annotations

import csv
import json
import subprocess
import zipfile
from pathlib import Path

import pytest

import mlvamaps.cli as cli
import mlvamaps.reference_pipeline as pipeline
from mlvamaps.io import read_loci
from mlvamaps.reference_pipeline import (
    TaxonReference,
    build_combined_taxon_database,
    build_taxon_references,
    ensure_combined_taxon_database,
    prepare_taxon_reference,
    read_taxon_references,
)


def _write_taxon_database(root, reference_id, sequence, taxid, taxon_name):
    database = root / "database"
    database.mkdir(parents=True)
    (database / "L1.fasta").write_text(f">{reference_id}\n{sequence}\n")
    (database / "reference_metadata.tsv").write_text(
        "reference_id\ttaxid\torganism_name\n"
        f"{reference_id}\t{taxid}\t{taxon_name}\n"
    )
    (database / "reference_assemblies.tsv").write_text(
        "reference_id\tassembly_file\tassembly_sha256\n"
        f"{reference_id}\t/{reference_id}.fna\tdigest-{reference_id}\n"
    )
    return {"database": str(database)}


def test_combined_taxon_database_is_ready_for_default_call(tmp_path, monkeypatch):
    panel = tmp_path / "panel.csv"
    panel.write_text("name,forward,reverse\nL1,AAA,CCC\n")
    results = [
        _write_taxon_database(tmp_path / "a", "R1", "AAATTTCCC", "1", "Species one"),
        _write_taxon_database(tmp_path / "b", "R2", "AAAGGGCCC", "2", "Species two"),
    ]
    observed = {}

    def fake_phylogeny(database, phylogeny, loci, threads, **kwargs):
        observed.update(database=database, phylogeny=phylogeny, loci=loci, threads=threads)
        Path(phylogeny).mkdir(parents=True)
        return {"phylogeny": Path(phylogeny)}

    monkeypatch.setattr(pipeline, "build_reference_phylogenies", fake_phylogeny)
    output = tmp_path / "combined"
    combined = build_combined_taxon_database(
        [TaxonReference("1", "taxon_one"), TaxonReference("2", "taxon_two")],
        results,
        panel,
        output,
        threads=2,
    )

    with (combined["database"] / "reference_metadata.tsv").open() as handle:
        metadata = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["reference_id"], row["taxon_id"], row["taxon_name"]) for row in metadata] == [
        ("R1", "1", "taxon_one"),
        ("R2", "2", "taxon_two"),
    ]
    assert (combined["database"] / "reference_panel.tsv").is_file()
    assert (combined["database"] / "mlva_contexts.tsv").is_file()
    assert (combined["database"] / "mlva_contexts.fasta.gz").is_file()
    records = pipeline.read_fasta(combined["database"] / "L1.fasta.gz")
    assert [name for name, _sequence in records] == [
        "R1",
        "R2",
    ]
    assert observed["database"] == combined["database"]

    args = cli.build_parser().parse_args(
        ["call", "-i", "sample.fasta", "--database", str(output)]
    )
    cli._resolve_call_args(cli.build_parser(), args)
    assert args.loci == str(combined["database"] / "reference_panel.tsv")
    assert args.taxon_identification is None


def test_combined_taxon_database_rejects_cross_taxon_reference_collisions(
    tmp_path, monkeypatch
):
    panel = tmp_path / "panel.csv"
    panel.write_text("name,forward,reverse\nL1,AAA,CCC\n")
    results = [
        _write_taxon_database(tmp_path / "a", "R1", "AAATTTCCC", "1", "Species one"),
        _write_taxon_database(tmp_path / "b", "R1", "AAAGGGCCC", "2", "Species two"),
    ]
    monkeypatch.setattr(pipeline, "build_reference_phylogenies", lambda *args, **kwargs: {})

    with pytest.raises(ValueError, match="occurs in both taxon"):
        build_combined_taxon_database(
            [TaxonReference("1", "taxon_one"), TaxonReference("2", "taxon_two")],
            results,
            panel,
            tmp_path / "combined",
        )


def test_combined_taxon_database_excludes_metadata_without_usable_loci(
    tmp_path, monkeypatch
):
    panel = tmp_path / "panel.csv"
    panel.write_text("name,forward,reverse\nL1,AAA,CCC\n")
    usable = _write_taxon_database(
        tmp_path / "a", "R1", "AAATTTCCC", "1", "Species one"
    )
    unusable = _write_taxon_database(
        tmp_path / "b", "R2", "AAAGGGCCC", "2", "Species two"
    )
    Path(unusable["database"], "L1.fasta").unlink()
    monkeypatch.setattr(
        pipeline,
        "build_reference_phylogenies",
        lambda database, phylogeny, *args, **kwargs: {"phylogeny": Path(phylogeny)},
    )

    combined = build_combined_taxon_database(
        [TaxonReference("1", "taxon_one"), TaxonReference("2", "taxon_two")],
        [usable, unusable],
        panel,
        tmp_path / "combined",
    )
    with (combined["database"] / "reference_metadata.tsv").open() as handle:
        metadata = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["reference_id"] for row in metadata] == ["R1"]


def test_existing_multi_taxid_build_is_upgraded_for_call(tmp_path, monkeypatch):
    root = tmp_path / "references"
    panel = root / "taxon_one" / "reference" / "database" / "reference_panel.tsv"
    results = [
        _write_taxon_database(root / "taxon_one" / "reference", "R1", "AAATTTCCC", "1", "one"),
        _write_taxon_database(root / "taxon_two" / "reference", "R2", "AAAGGGCCC", "2", "two"),
    ]
    panel.write_text("name,forward,reverse\nL1,AAA,CCC\n")
    (root / "reference_pipeline_manifest.json").write_text(
        json.dumps(
            {
                "references": [
                    {"taxid": "1", "name": "taxon_one", "database": results[0]["database"]},
                    {"taxid": "2", "name": "taxon_two", "database": results[1]["database"]},
                ]
            }
        )
    )
    monkeypatch.setattr(
        pipeline,
        "build_reference_phylogenies",
        lambda database, phylogeny, *args, **kwargs: {"phylogeny": Path(phylogeny)},
    )

    assert ensure_combined_taxon_database(root) == root.resolve()
    assert (root / "database" / "reference_panel.tsv").is_file()
    assert ensure_combined_taxon_database(root) == root.resolve()


def test_call_help_keeps_automatic_taxon_workflow_simple(capsys):
    with pytest.raises(SystemExit, match="0"):
        cli.build_parser().parse_args(["call", "--help"])
    help_text = capsys.readouterr().out
    assert "--database DATABASE" in help_text
    assert "--taxon-identification" not in help_text
    assert "--target-taxon-id" not in help_text
    assert "--taxon-calibration" not in help_text


def test_read_taxon_references_from_csv_with_optional_names(tmp_path):
    taxids = tmp_path / "taxids.csv"
    taxids.write_text(
        "taxid,name\n86661,b_cereus_group\n1280,staphylococcus_aureus\n"
    )

    assert read_taxon_references(taxids_csv=taxids) == [
        TaxonReference("86661", "b_cereus_group"),
        TaxonReference("1280", "staphylococcus_aureus"),
    ]
    assert read_taxon_references(taxid="001280") == [
        TaxonReference("1280", "taxid_1280")
    ]


def test_taxid_csv_rejects_duplicate_taxids_and_unsafe_names(tmp_path):
    duplicate = tmp_path / "duplicate.csv"
    duplicate.write_text("taxid\n1280\n1280\n")
    with pytest.raises(ValueError, match="duplicate taxids"):
        read_taxon_references(taxids_csv=duplicate)

    unsafe = tmp_path / "unsafe.csv"
    unsafe.write_text("taxid,name\n1280,../shared\n")
    with pytest.raises(ValueError, match="invalid reference name"):
        read_taxon_references(taxids_csv=unsafe)


def test_rich_loci_csv_uses_the_same_schema_as_call(tmp_path):
    loci_csv = tmp_path / "loci.csv"
    loci_csv.write_text(
        "name,forward,reverse,repeat_motif,"
        "expected_min_repeats,expected_max_repeats\n"
        "VNTR_1,AAACCC,GGGTTT,AT,2,20\n"
    )

    loci = read_loci(loci_csv)

    assert len(loci) == 1
    assert loci[0].locus_id == "VNTR_1"
    assert loci[0].forward_primer == "AAACCC"
    assert loci[0].reverse_primer == "GGGTTT"
    assert loci[0].repeat_motif == "AT"
    assert loci[0].expected_min_repeats == 2
    assert loci[0].expected_max_repeats == 20


def test_loci_csv_reports_available_columns_instead_of_key_error(tmp_path):
    loci_csv = tmp_path / "loci.csv"
    loci_csv.write_text("marker,forward,reverse\nVNTR_1,AAA,CCC\n")

    with pytest.raises(ValueError, match=r"needs a locus_id.*found: marker"):
        read_loci(loci_csv)


def test_ncbi_download_retries_after_a_partial_stream_failure(tmp_path, monkeypatch):
    package = tmp_path / "ncbi_dataset.zip"
    command = ["datasets", "download", "--filename", str(package)]
    attempts = []

    def fake_run(_command):
        attempts.append(len(attempts) + 1)
        if len(attempts) == 1:
            package.write_bytes(b"partial download")
            raise RuntimeError("stream error: INTERNAL_ERROR")
        assert not package.exists()
        with zipfile.ZipFile(package, "w") as archive:
            archive.writestr("ncbi_dataset/data/README.md", "ok")
        return subprocess.CompletedProcess(_command, 0, "", "")

    monkeypatch.setattr(pipeline, "_run_command", fake_run)
    monkeypatch.setattr(pipeline.time, "sleep", lambda _seconds: None)

    pipeline._download_package(command, package, attempts=3)

    assert attempts == [1, 2]
    assert zipfile.is_zipfile(package)


def test_prepare_taxon_reference_downloads_and_normalizes_ncbi_package(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(
        pipeline.shutil, "which", lambda executable: f"/tools/{executable}"
    )
    monkeypatch.setattr(pipeline, "_tool_version", lambda executable: "test 1.0")

    raw_metadata = (
        "Assembly Accession\tOrganism Name\tOrganism Taxonomic ID\tAssembly Level\t"
        "Assembly BioSample Accession\tAssembly BioSample Geographic Location\t"
        "Assembly BioSample Latitude Longitude\n"
        "GCF_000001.1\tExample bacterium\t1280\tComplete Genome\tSAMN1\t"
        "USA: New York\t40.7 N 74.0 W\n"
    )

    def fake_run(command):
        if command[1:4] == ["download", "genome", "taxon"]:
            package = Path(command[command.index("--filename") + 1])
            package.parent.mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(package, "w") as archive:
                archive.writestr(
                    "ncbi_dataset/data/GCF_000001.1/GCF_000001.1_genomic.fna",
                    ">contig\nACGT\n",
                )
            return subprocess.CompletedProcess(command, 0, "", "")
        return subprocess.CompletedProcess(command, 0, raw_metadata, "")

    monkeypatch.setattr(pipeline, "_run_command", fake_run)
    result = prepare_taxon_reference(
        TaxonReference("1280", "s_aureus"),
        tmp_path / "prepared",
    )

    with result["metadata"].open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["accession"] == "GCF_000001.1"
    assert rows[0]["assembly_file"].endswith("GCF_000001.1_genomic.fna")
    assert rows[0]["latitude"] == "40.7"
    assert rows[0]["longitude"] == "-74.0"
    manifest = json.loads(result["manifest"].read_text())
    assert manifest["selection"] == {
        "assembly_source": "refseq",
        "taxid": "1280",
    }
    assert manifest["assembly_count"] == 1
    assert manifest["metadata_count"] == 1


def test_build_taxon_references_keeps_each_database_isolated(tmp_path, monkeypatch):
    def fake_prepare(reference, output, **kwargs):
        output = Path(output)
        genomes = output / "package"
        genomes.mkdir(parents=True)
        metadata = output / "metadata.tsv"
        metadata.write_text("accession\tassembly_file\n")
        return {
            "outdir": output,
            "genomes": genomes,
            "metadata": metadata,
            "manifest": output / "download_manifest.json",
            "taxid": reference.taxid,
            "name": reference.name,
        }

    observed = []

    def fake_builder(**kwargs):
        observed.append(kwargs)
        output = Path(kwargs["outdir"])
        return {
            "outdir": output,
            "database": output / "database",
            "manifest": output / "reference_build_manifest.tsv",
        }

    monkeypatch.setattr(pipeline, "prepare_taxon_reference", fake_prepare)
    result = build_taxon_references(
        [
            TaxonReference("1280", "staph"),
            TaxonReference("86661", "bcg"),
        ],
        "primers.csv",
        tmp_path / "references",
        builder=fake_builder,
        threads=4,
    )

    assert [Path(call["outdir"]) for call in observed] == [
        tmp_path / "references" / "staph" / "reference",
        tmp_path / "references" / "bcg" / "reference",
    ]
    assert observed[0]["metadata_path"] == (
        tmp_path / "references" / "staph" / "prepared" / "metadata.tsv"
    )
    assert len(json.loads(result["manifest"].read_text())["references"]) == 2
    assert [entry["status"] for entry in result["references"]] == ["BUILT", "BUILT"]


def test_multi_taxon_build_continues_after_no_usable_loci(tmp_path, monkeypatch):
    references = [
        TaxonReference("1", "taxon_a"),
        TaxonReference("2", "taxon_b"),
        TaxonReference("3", "taxon_c"),
    ]

    def fake_prepare(reference, output, **kwargs):
        output = Path(output)
        genomes = output / "package"
        genomes.mkdir(parents=True)
        metadata = output / "metadata.tsv"
        metadata.write_text("accession\tassembly_file\n")
        return {"outdir": output, "genomes": genomes, "metadata": metadata}

    built_names = []

    def fake_builder(**kwargs):
        output = Path(kwargs["outdir"])
        name = output.parent.name
        built_names.append(name)
        return {
            "status": "NO_USABLE_LOCI" if name == "taxon_b" else "BUILT",
            "outdir": output,
            "database": output / "database",
            "manifest": output / "reference_build_manifest.tsv",
        }

    monkeypatch.setattr(pipeline, "prepare_taxon_reference", fake_prepare)
    result = build_taxon_references(
        references,
        "primers.csv",
        tmp_path / "references",
        builder=fake_builder,
    )

    assert built_names == ["taxon_a", "taxon_b", "taxon_c"]
    assert [entry["status"] for entry in result["references"]] == [
        "BUILT",
        "NO_USABLE_LOCI",
        "BUILT",
    ]
    manifest = json.loads(result["manifest"].read_text())
    assert manifest["references"][1]["name"] == "taxon_b"
    assert manifest["references"][1]["status"] == "NO_USABLE_LOCI"


def test_multi_taxon_amplifiability_summaries_and_console(tmp_path, monkeypatch, capsys):
    references = [TaxonReference("1", "taxon_a"), TaxonReference("2", "taxon_b")]

    def fake_prepare(reference, output, **kwargs):
        output = Path(output)
        genomes = output / "package"
        genomes.mkdir(parents=True)
        metadata = output / "metadata.tsv"
        metadata.write_text("accession\tassembly_file\n")
        return {"outdir": output, "genomes": genomes, "metadata": metadata}

    def fake_builder(**kwargs):
        output = Path(kwargs["outdir"])
        rows = {
            "taxon_a": [
                {"locus_id": "L1", "genomes_examined": 10, "genomes_with_valid_amplicon": 10,
                 "valid_amplicons": 10, "percent_genomes_amplifiable": 100.0,
                 "amplifiable": "TRUE", "tree_status": "BUILT"},
                {"locus_id": "L2", "genomes_examined": 10, "genomes_with_valid_amplicon": 8,
                 "valid_amplicons": 8, "percent_genomes_amplifiable": 80.0,
                 "amplifiable": "TRUE", "tree_status": "BUILT"},
            ],
            "taxon_b": [
                {"locus_id": "L1", "genomes_examined": 4, "genomes_with_valid_amplicon": 2,
                 "valid_amplicons": 2, "percent_genomes_amplifiable": 50.0,
                 "amplifiable": "TRUE", "tree_status": "INSUFFICIENT_REFERENCES"},
                {"locus_id": "L2", "genomes_examined": 4, "genomes_with_valid_amplicon": 0,
                 "valid_amplicons": 0, "percent_genomes_amplifiable": 0.0,
                 "amplifiable": "FALSE", "tree_status": "NO_AMPLICONS"},
            ],
        }[output.parent.name]
        return {
            "outdir": output,
            "database": output / "database",
            "manifest": output / "reference_build_manifest.tsv",
            "locus_summary_rows": rows,
        }

    monkeypatch.setattr(pipeline, "prepare_taxon_reference", fake_prepare)
    result = build_taxon_references(
        references, "primers.csv", tmp_path / "references", builder=fake_builder,
        show_progress=True,
    )

    taxon_rows = list(csv.DictReader(result["taxon_summary"].open(), delimiter="\t"))
    locus_rows = list(csv.DictReader(result["locus_amplifiability"].open(), delimiter="\t"))
    assert list(taxon_rows[0]) == pipeline.TAXON_REFERENCE_SUMMARY_FIELDS
    assert list(locus_rows[0]) == pipeline.TAXON_LOCUS_AMPLIFIABILITY_FIELDS
    assert [(row["taxon"], row["loci_amplifiable"], row["status"]) for row in taxon_rows] == [
        ("taxon_a", "2", "BUILT"),
        ("taxon_b", "1", "PARTIAL"),
    ]
    assert taxon_rows[1]["percent_loci_amplifiable"] == "50.0"
    assert taxon_rows[1]["total_valid_amplicons"] == "2"
    assert len(locus_rows) == 4
    output = capsys.readouterr().out
    assert "taxon_b [taxid 2]" in output
    assert "Amplifiable loci: 1 / 2 (50.0%)" in output
    assert "Non-amplifiable loci: L2" in output


def test_quiet_taxon_build_suppresses_console_summary(tmp_path, monkeypatch, capsys):
    def fake_prepare(reference, output, **kwargs):
        output = Path(output)
        output.mkdir(parents=True)
        return {"outdir": output, "genomes": output, "metadata": output / "metadata.tsv"}

    def fake_builder(**kwargs):
        output = Path(kwargs["outdir"])
        return {
            "outdir": output,
            "database": output / "database",
            "manifest": output / "manifest.tsv",
            "locus_summary_rows": [
                {"locus_id": "L1", "genomes_examined": 1, "genomes_with_valid_amplicon": 1,
                 "valid_amplicons": 1, "percent_genomes_amplifiable": 100.0,
                 "amplifiable": "TRUE", "tree_status": "INSUFFICIENT_REFERENCES"}
            ],
        }

    monkeypatch.setattr(pipeline, "prepare_taxon_reference", fake_prepare)
    build_taxon_references(
        [TaxonReference("1", "taxon_a")], "primers.csv", tmp_path / "references",
        builder=fake_builder, show_progress=False,
    )
    assert capsys.readouterr().out == ""


def test_multi_taxon_build_does_not_swallow_unexpected_builder_failure(
    tmp_path, monkeypatch
):
    def fake_prepare(reference, output, **kwargs):
        output = Path(output)
        output.mkdir(parents=True)
        return {
            "outdir": output,
            "genomes": output / "package",
            "metadata": output / "metadata.tsv",
        }

    def fail_builder(**kwargs):
        raise RuntimeError("programming error")

    monkeypatch.setattr(pipeline, "prepare_taxon_reference", fake_prepare)
    with pytest.raises(RuntimeError, match="programming error"):
        build_taxon_references(
            [TaxonReference("1", "taxon_a")],
            "primers.csv",
            tmp_path / "references",
            builder=fail_builder,
        )


def test_cli_dispatches_taxid_reference_pipeline(tmp_path, monkeypatch):
    observed = {}
    panel = tmp_path / "primers.csv"
    panel.write_text("locus_id,forward_primer,reverse_primer\nL1,AAA,CCC\n")

    def fake_build(references, primers_path, outdir, **kwargs):
        observed.update(
            references=references,
            primers_path=primers_path,
            outdir=outdir,
            kwargs=kwargs,
        )
        return {
            "manifest": tmp_path / "reference_pipeline_manifest.json",
            "references": [
                {
                    "taxid": "1280",
                    "name": "taxid_1280",
                    "database": str(tmp_path / "taxid_1280" / "reference" / "database"),
                }
            ],
        }

    monkeypatch.setattr(cli, "build_taxon_references", fake_build)
    assert (
        cli.main(
            [
                "build-reference",
                "--taxid",
                "1280",
                    "-p",
                    str(panel),
                "--output",
                str(tmp_path / "references"),
                "--quiet",
            ]
        )
        == 0
    )
    assert observed["references"] == [TaxonReference("1280", "taxid_1280")]
    assert observed["outdir"] == str(tmp_path / "references")
    assert observed["kwargs"]["show_progress"] is False
