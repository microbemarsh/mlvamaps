from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

import mlvamaps.reference_builder as reference_builder
from mlvamaps.io import read_fasta
from mlvamaps.reference_builder import build_reference_database

from test_phylogeny import _fake_mafft, _fake_raxml_ng


def _rows(path: Path, delimiter: str = "\t") -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


@pytest.fixture(autouse=True)
def _stub_external_mapping_resources(monkeypatch):
    def fake_resources(*, database_dir, **kwargs):
        competitive = Path(database_dir) / "competitive_mapping"
        deacon = Path(database_dir) / "deacon"
        competitive.mkdir(parents=True, exist_ok=True)
        deacon.mkdir(parents=True, exist_ok=True)
        paths = {
            "fasta": competitive / "candidate_contexts.fasta",
            "metadata": competitive / "candidate_metadata.tsv",
            "provenance": competitive / "candidate_provenance.json",
            "provenance_table": competitive / "candidate_provenance.tsv",
            "short_index": competitive / "short.mmi",
            "long_index": competitive / "long.mmi",
            "deacon_reference": deacon / "reference_genomes.fasta",
            "deacon_index": deacon / "target_recruitment.idx",
        }
        for path in paths.values():
            path.write_text("test\n")
        return {**paths, "minimap2_version": "test", "deacon_version": "test"}

    monkeypatch.setattr(reference_builder, "build_mapping_resources", fake_resources)
    monkeypatch.setattr(reference_builder, "_tool_version", lambda _executable: "test")


def _reference_inputs(tmp_path: Path) -> tuple[Path, Path, Path]:
    assemblies = tmp_path / "assemblies"
    assemblies.mkdir()
    (assemblies / "R1.fasta").write_text(">R1\nACGT\n")
    primers = tmp_path / "primers.csv"
    primers.write_text("name,forward,reverse\nL1,AAA,CCC\n")
    metadata = tmp_path / "metadata.csv"
    metadata.write_text("reference_id\nR1\n")
    return assemblies, primers, metadata


def _mock_empty_extraction(monkeypatch) -> None:
    monkeypatch.setattr(
        reference_builder,
        "run_in_silico_pcr_loci",
        lambda assembly, loci, outdir, **kwargs: {
            "stats": Path(outdir) / "stats.tsv",
            "products": Path(outdir) / "products.fasta",
        },
    )
    monkeypatch.setattr(reference_builder, "read_pcr_results", lambda *args: [])
    monkeypatch.setattr(reference_builder, "pcr_rows_to_products", lambda *args: [])


def test_build_reference_database_from_assemblies_and_metadata(tmp_path, monkeypatch):
    assemblies = tmp_path / "assemblies"
    assemblies.mkdir()
    for sample in ("R1", "R2", "R3"):
        (assemblies / f"{sample}.fasta").write_text(f">{sample}\nACGT\n")
    primers = tmp_path / "primers.csv"
    primers.write_text("name,forward,reverse\nL1,AAA,CCC\nL2,GGG,TTT\n")
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(
        "sample_id,collection_date,latitude,longitude\n"
        "R1,2020-01-01,1,2\nR2,2020-01-02,3,4\nR3,2020-01-03,5,6\n"
    )

    monkeypatch.setattr(
        "mlvamaps.reference_builder.run_in_silico_pcr_loci",
        lambda assembly, loci, outdir, **kwargs: {
            "stats": Path(outdir) / "stats.tsv",
            "products": Path(outdir) / "products.fasta",
        },
    )
    monkeypatch.setattr("mlvamaps.reference_builder.read_pcr_results", lambda *args: [])

    def fake_products(rows, loci, sample_id):
        products = []
        for index, locus in enumerate(loci):
            # Make R3/L2 ambiguous so the default QC policy excludes it.
            copies = 2 if sample_id == "R3" and locus.locus_id == "L2" else 1
            for copy in range(copies):
                products.append(
                    {
                        "locus_id": locus.locus_id,
                        "product_id": f"{locus.locus_id}|{sample_id}|{copy}",
                        "sequence": "AAACCC" if index == 0 else "GGGTTT",
                        "product_size_bp": 6 + copy,
                        "forward_mismatches": 0,
                        "reverse_mismatches": 0,
                        "primer_error_round": 0,
                    }
                )
        return products

    monkeypatch.setattr("mlvamaps.reference_builder.pcr_rows_to_products", fake_products)
    result = build_reference_database(
        assemblies,
        primers,
        metadata,
        tmp_path / "reference",
        threads=1,
        min_references_per_tree=2,
        mafft_bin=str(_fake_mafft(tmp_path)),
        raxml_ng_bin=str(_fake_raxml_ng(tmp_path)),
    )

    assert len(list(read_fasta(result["database"] / "L1.fasta.gz"))) == 3
    assert len(list(read_fasta(result["database"] / "L2.fasta.gz"))) == 2
    assert (result["phylogeny"] / "L1" / "reference_tree.nwk").exists()
    assert (result["phylogeny"] / "L2" / "reference_tree.nwk").exists()
    assert (result["phylogeny"] / "L1.tree").exists()
    assert len(_rows(result["reference_assemblies"])) == 3
    assert _rows(result["reference_assemblies"])[0]["assembly_sha256"]
    manifest = _rows(result["build_qc"])
    ambiguous = [row for row in manifest if row["reference_id"] == "R3" and row["locus_id"] == "L2"]
    assert ambiguous[0]["status"] == "AMBIGUOUS_EXCLUDED"
    assert ambiguous[0]["best_product_count"] == "2"
    uncompressed_sequences = [
        path
        for path in result["outdir"].rglob("*")
        if path.is_file() and path.suffix.lower() in {".fasta", ".fa", ".fastq", ".fq"}
        and path.name not in {"candidate_contexts.fasta", "reference_genomes.fasta"}
    ]
    assert uncompressed_sequences == []
    assert _rows(result["metadata"])[0]["reference_id"] == "R1"
    assert _rows(result["myoga_metadata"], ",")[0]["genome_id"] == "R1"
    assert (result["database"] / "competitive_mapping" / "candidate_metadata.tsv").is_file()
    assert (result["database"] / "competitive_mapping" / "candidate_contexts.fasta").is_file()
    assert (result["database"] / "competitive_mapping" / "short.mmi").is_file()
    assert (result["database"] / "competitive_mapping" / "long.mmi").is_file()
    assert (result["database"] / "deacon" / "target_recruitment.idx").is_file()
    assert result["status"] == "BUILT"
    build_manifest = json.loads(result["manifest"].read_text())
    assert build_manifest["schema_version"] == "2.0"
    assert build_manifest["status"] == "complete"
    assert build_manifest["competitive_mapping"]["indexes"] == {
        "long": {"k": 11, "w": 5},
        "short": {"k": 21, "w": 11},
    }
    assert build_manifest["phylogeny"]["sequence_type"] == (
        "real_observed_locus_sequences_only"
    )
    locus_summary = _rows(result["locus_amplifiability"])
    assert list(locus_summary[0]) == reference_builder.REFERENCE_LOCUS_AMPLIFIABILITY_FIELDS
    assert [(row["locus_id"], row["valid_amplicons"], row["amplifiable"]) for row in locus_summary] == [
        ("L1", "3", "TRUE"),
        ("L2", "2", "TRUE"),
    ]
    assert [row["percent_genomes_amplifiable"] for row in locus_summary] == ["100.0", "66.7"]
    assert {row["tree_status"] for row in locus_summary} == {"BUILT"}


def test_empty_reference_database_skips_phylogeny_and_preserves_qc_outputs(
    tmp_path, monkeypatch, capsys
):
    assemblies, primers, metadata = _reference_inputs(tmp_path)
    _mock_empty_extraction(monkeypatch)

    def unexpected_phylogeny(*args, **kwargs):
        pytest.fail("phylogeny construction must not run without locus FASTAs")

    monkeypatch.setattr(reference_builder, "build_reference_phylogenies", unexpected_phylogeny)

    result = build_reference_database(
        assemblies,
        primers,
        metadata,
        tmp_path / "reference",
        threads=1,
        show_progress=True,
    )

    assert result["status"] == "NO_USABLE_LOCI"
    assert result["phylogeny"] is None
    assert not list(result["database"].glob("*.fasta*"))
    assert result["manifest"].is_file()
    assert result["metadata"].is_file()
    assert _rows(result["build_qc"])[0]["status"] == "NOT_FOUND"
    assert _rows(result["locus_amplifiability"])[0] == {
        "locus_id": "L1",
        "genomes_examined": "1",
        "genomes_with_valid_amplicon": "0",
        "valid_amplicons": "0",
            "genomes_with_multiple_products": "0",
            "genomes_failing_product_constraints": "1",
        "percent_genomes_amplifiable": "0.0",
        "amplifiable": "FALSE",
        "tree_status": "NO_AMPLICONS",
    }
    assert (
        "No usable reference amplicons were recovered; skipping phylogeny construction."
        in capsys.readouterr().err
    )


def test_unexpected_phylogeny_failure_is_not_swallowed(tmp_path, monkeypatch):
    assemblies, primers, metadata = _reference_inputs(tmp_path)
    _mock_empty_extraction(monkeypatch)
    monkeypatch.setattr(
        reference_builder,
        "pcr_rows_to_products",
        lambda rows, loci, sample_id: [
            {
                "locus_id": "L1",
                "product_id": "L1|R1|0",
                "sequence": "AAACCC",
                "product_size_bp": 6,
                "forward_mismatches": 0,
                "reverse_mismatches": 0,
                "primer_error_round": 0,
            }
        ],
    )

    def fail_phylogeny(*args, **kwargs):
        raise RuntimeError("unrelated tree failure")

    monkeypatch.setattr(reference_builder, "build_reference_phylogenies", fail_phylogeny)

    with pytest.raises(RuntimeError, match="unrelated tree failure"):
        build_reference_database(
            assemblies,
            primers,
            metadata,
            tmp_path / "reference",
            threads=1,
        )


def test_valid_amplicon_below_tree_minimum_remains_amplifiable(tmp_path, monkeypatch):
    assemblies, primers, metadata = _reference_inputs(tmp_path)
    _mock_empty_extraction(monkeypatch)
    monkeypatch.setattr(
        reference_builder,
        "pcr_rows_to_products",
        lambda rows, loci, sample_id: [
            {
                "locus_id": "L1",
                "product_id": "L1|R1|0",
                "sequence": "AAACCC",
                "product_size_bp": 6,
                "forward_mismatches": 0,
                "reverse_mismatches": 0,
                "primer_error_round": 0,
            }
        ],
    )
    monkeypatch.setattr(
        reference_builder,
        "build_reference_phylogenies",
        lambda *args, **kwargs: {"phylogeny": Path(args[1])},
    )

    result = build_reference_database(
        assemblies,
        primers,
        metadata,
        tmp_path / "reference",
        threads=1,
        min_references_per_tree=3,
    )

    row = _rows(result["locus_amplifiability"])[0]
    assert row["valid_amplicons"] == "1"
    assert row["amplifiable"] == "TRUE"
    assert row["tree_status"] == "INSUFFICIENT_REFERENCES"
    assert result["status"] == "PARTIAL"


def test_cli_exposes_reference_builder():
    from mlvamaps.cli import build_parser

    args = build_parser().parse_args(
        [
            "build-reference",
            "-i",
            "assemblies",
            "-p",
            "primers.csv",
            "--metadata",
            "metadata.csv",
        ]
    )
    assert args.multiple_products == "exclude"
    assert args.min_references_per_tree == 3
    assert args.raxml_model == "DNA"
    assert args.quiet is False


def test_reference_build_keeps_global_threads_but_pins_raxml_in_manifest(
    tmp_path, monkeypatch
):
    assemblies, primers, metadata = _reference_inputs(tmp_path)
    _mock_empty_extraction(monkeypatch)
    monkeypatch.setattr(
        reference_builder,
        "pcr_rows_to_products",
        lambda rows, loci, sample_id: [{
            "locus_id": "L1", "product_id": "L1|R1|0", "sequence": "AAACCC",
            "product_size_bp": 6, "forward_mismatches": 0,
            "reverse_mismatches": 0, "primer_error_round": 0,
        }],
    )
    observed = {}

    def fake_resources(**kwargs):
        observed["mapping_threads"] = kwargs["threads"]
        return {}

    def fake_phylogeny(database, outdir, loci, threads, **kwargs):
        observed["phylogeny_threads"] = threads
        return {"phylogeny": Path(outdir)}

    monkeypatch.setattr(reference_builder, "build_mapping_resources", fake_resources)
    monkeypatch.setattr(reference_builder, "build_reference_phylogenies", fake_phylogeny)
    result = build_reference_database(
        assemblies, primers, metadata, tmp_path / "reference", threads=32
    )

    assert observed == {"mapping_threads": 32, "phylogeny_threads": 32}
    manifest = json.loads(result["manifest"].read_text())
    assert manifest["build_threads"] == 32
    assert manifest["phylogeny"]["raxml_ng_threads_per_process"] == 1


def test_cli_help_lists_only_supported_commands(capsys):
    from mlvamaps.cli import build_parser

    with pytest.raises(SystemExit, match="0"):
        build_parser().parse_args(["--help"])
    help_text = capsys.readouterr().out
    assert "{call,export-myoga,build-reference,calibrate-taxa}" in help_text


def test_build_reference_accepts_primers_and_taxids_spelling():
    from mlvamaps.cli import build_parser

    args = build_parser().parse_args(
        [
            "build-reference", "--primers", "primer.csv", "--taxids", "taxids.csv",
            "--threads", "32",
        ]
    )
    assert args.panel_path == "primer.csv"
    assert args.taxids_csv == "taxids.csv"
    assert args.threads == 32


def test_reference_builder_normalizes_and_checks_taxonomy(tmp_path, monkeypatch):
    assemblies, primers, metadata = _reference_inputs(tmp_path)
    metadata.write_text("reference_id,taxid,organism_name\nR1,123,Species one\n")
    _mock_empty_extraction(monkeypatch)
    result = build_reference_database(
        assemblies, primers, metadata, tmp_path / "reference", threads=1
    )
    rows = _rows(result["metadata"])
    assert rows[0]["taxon_id"] == "123"
    assert rows[0]["taxon_name"] == "Species one"
    assert (result["database"] / "reference_panel.tsv").is_file()

    metadata.write_text("reference_id,taxid\nR1,\n")
    with pytest.raises(ValueError, match="blank taxon identifiers"):
        build_reference_database(
            assemblies, primers, metadata, tmp_path / "invalid", threads=1
        )


def test_cli_exposes_package_version(capsys):
    from mlvamaps import __version__
    from mlvamaps.cli import build_parser

    with pytest.raises(SystemExit, match="0"):
        build_parser().parse_args(["--version"])
    assert capsys.readouterr().out.strip() == f"mlvamaps {__version__}"


def test_reference_extraction_uses_multiple_processes_and_reports_progress(
    tmp_path, capsys
):
    assemblies = tmp_path / "assemblies"
    assemblies.mkdir()
    product = "ACGTACGTACGT" + "T" * 50 + "TAGCTAGCTAGC"
    for sample in ("R1", "R2"):
        (assemblies / f"{sample}.fasta").write_text(f">{sample}\n{product}\n")
    primers = tmp_path / "primers.csv"
    primers.write_text(
        "name,forward,reverse\nL1,ACGTACGTACGT,GCTAGCTAGCTA\n"
    )
    metadata = tmp_path / "metadata.csv"
    metadata.write_text("reference_id\nR1\nR2\n")

    result = build_reference_database(
        assemblies,
        primers,
        metadata,
        tmp_path / "reference",
        threads=2,
        min_references_per_tree=2,
        mafft_bin=str(_fake_mafft(tmp_path)),
        raxml_ng_bin=str(_fake_raxml_ng(tmp_path)),
        show_progress=True,
    )

    assert len(list(read_fasta(result["database"] / "L1.fasta.gz"))) == 2
    progress = capsys.readouterr().err
    assert "with 2 worker(s)" in progress
    assert "Extracted assemblies: 2/2 (100.0%)" in progress
    assert "Processed tree loci: 1/1 (100.0%)" in progress
