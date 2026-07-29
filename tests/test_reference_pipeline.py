from __future__ import annotations

import csv
import json
import subprocess
import zipfile
from pathlib import Path

import pytest

import mlvamaps.cli as cli
import mlvamaps.reference_pipeline as pipeline
from mlvamaps.reference_pipeline import (
    TaxonReference,
    build_taxon_references,
    prepare_taxon_reference,
    read_taxon_references,
)


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


def test_cli_dispatches_taxid_reference_pipeline(tmp_path, monkeypatch):
    observed = {}

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
                "--primers",
                "primers.csv",
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
