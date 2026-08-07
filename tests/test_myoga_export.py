from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

import numpy as np
import pytest

import mlvamaps.cli as cli
from mlvamaps.myoga_export import (
    SampleCalls,
    calculate_pairwise_distances,
    export_myoga,
    read_export_metadata,
)
from mlvamaps.combined_marker_export import recover_masked_sequences
from mlvamaps.models import Locus
from mlvamaps.phylogeny import (
    _parse_newick,
    _tip_patristic_distance_matrix,
    neighbor_joining_tree_from_matrix,
)


CALL_FIELDS = [
    "sample_id",
    "locus_id",
    "present",
    "repeat_count",
    "repeat_count_raw",
    "product_size_bp",
    "read_depth",
    "primary_read_depth",
    "mean_coverage",
    "allele_confidence",
    "status",
    "evidence",
]


def _write_calls(root: Path, sample_id: str, alleles: dict[str, object]) -> Path:
    sample = root / sample_id
    sample.mkdir(parents=True, exist_ok=True)
    path = sample / "calls.tsv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CALL_FIELDS, delimiter="\t")
        writer.writeheader()
        for locus, allele in alleles.items():
            exact = allele not in (None, "")
            writer.writerow(
                {
                    "sample_id": sample_id,
                    "locus_id": locus,
                    "present": "yes",
                    "repeat_count": "" if allele is None else allele,
                    "repeat_count_raw": "" if allele is None else allele,
                    "product_size_bp": 100 if exact else "",
                    "read_depth": 10,
                    "primary_read_depth": 8 if exact else 0,
                    "mean_coverage": 20,
                    "allele_confidence": 0.99 if exact else "",
                    "status": "PASS" if exact else "PRESENCE_ONLY",
                    "evidence": "COMPLETE_ASSEMBLED_PRODUCT" if exact else "PRESENCE_ONLY",
                }
            )
    return path


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _metadata(path: Path, delimiter: str = "\t") -> Path:
    path.write_text(
        delimiter.join(["shared_identifier", "latitude", "longitude", "host", "country"])
        + "\n"
        + delimiter.join(["S1", "10", "20", "human", "US"])
        + "\n"
        + delimiter.join(["S2", "", "", "cow", "CA"])
        + "\n"
        + delimiter.join(["S3", "95", "200", "bird", "MX"])
        + "\n",
        encoding="utf-8",
    )
    return path


def _basic_results(tmp_path: Path) -> tuple[Path, Path]:
    results = tmp_path / "results"
    _write_calls(results, "S1", {"L1": 1, "L2": 2, "L3": 3})
    _write_calls(results, "S2", {"L1": 1, "L2": 4, "L3": None})
    _write_calls(results, "S3", {"L1": 2, "L2": 2, "L3": 5})
    return results, _metadata(tmp_path / "metadata.tsv")


def test_successful_export_builds_long_wide_metadata_distances_and_tree(tmp_path):
    results, metadata = _basic_results(tmp_path)
    output = tmp_path / "export"

    result = export_myoga(
        results,
        metadata,
        output,
        min_callable_fraction=2 / 3,
        min_pairwise_fraction=2 / 3,
    )

    assert result["tree_samples"] == 3
    profiles = _read_tsv(output / "mlva_profiles.tsv")
    assert [row["sample_id"] for row in profiles] == ["S1", "S2", "S3"]
    assert profiles[1] == {"sample_id": "S2", "L1": "1", "L2": "4", "L3": ""}
    long_rows = _read_tsv(output / "mlva_calls_long.tsv")
    assert len(long_rows) == 9
    unresolved = next(row for row in long_rows if row["sample_id"] == "S2" and row["locus_id"] == "L3")
    assert unresolved["present"] == "yes"
    assert unresolved["repeat_count"] == ""
    assert unresolved["status"] == "PRESENCE_ONLY"

    pairs = {(row["sample_1"], row["sample_2"]): row for row in _read_tsv(output / "mlva_pairwise_distances.tsv")}
    assert pairs[("S1", "S2")]["loci_compared"] == "2"
    assert pairs[("S1", "S2")]["categorical_differences"] == "1"
    assert pairs[("S1", "S2")]["categorical_distance"] == "0.50000000"
    assert pairs[("S1", "S2")]["repeat_distance_raw"] == "2.00000000"
    assert pairs[("S1", "S2")]["repeat_distance"] == "1.00000000"

    matrix = _read_tsv(output / "mlva_distance_matrix.tsv")
    assert [row["sample_id"] for row in matrix] == ["S1", "S2", "S3"]
    for row in matrix:
        assert row[row["sample_id"]] == "0.00000000"
    assert matrix[0]["S3"] == matrix[2]["S1"] == "1.00000000"
    tips = set(re.findall(r"'([^']+)'", (output / "mlva_nj.tree").read_text()))
    metadata_ids = {row["sample_id"] for row in _read_tsv(output / "myoga_metadata.tsv")}
    assert tips == metadata_ids == {"S1", "S2", "S3"}


def test_metadata_tsv_preserves_fields_and_records_missing_invalid_coordinates(tmp_path):
    results, metadata = _basic_results(tmp_path)
    output = tmp_path / "export"
    export_myoga(results, metadata, output, min_callable_fraction=2 / 3)

    rows = {row["sample_id"]: row for row in _read_tsv(output / "myoga_metadata.tsv")}
    assert rows["S1"]["host"] == "human"
    assert rows["S1"]["latitude"] == "10"
    assert rows["S2"]["latitude"] == ""
    assert rows["S3"]["latitude"] == ""
    assert rows["S3"]["metadata_latitude"] == "95"
    issues = {(row["sample_id"], row["reason"], row["scope"]) for row in _read_tsv(output / "samples_excluded.tsv")}
    assert ("S2", "MISSING_COORDINATES", "geography") in issues
    assert ("S3", "INVALID_COORDINATES", "geography") in issues
    assert {row["sample_id"] for row in _read_tsv(output / "samples_used.tsv")} == {"S1", "S2", "S3"}


def test_metadata_csv_custom_id_and_coordinate_columns(tmp_path):
    results = tmp_path / "results"
    _write_calls(results, "SRR1", {"L1": 2})
    metadata = tmp_path / "metadata.csv"
    metadata.write_text("run,lat,lon,note\nSRR1,1.5,-2.5,contains comma-like text\n")
    output = tmp_path / "export"

    export_myoga(
        results,
        metadata,
        output,
        metadata_id="run",
        latitude="lat",
        longitude="lon",
        min_callable_fraction=1,
    )

    row = _read_tsv(output / "myoga_metadata.tsv")[0]
    assert row["sample_id"] == "SRR1"
    assert row["run"] == "SRR1"
    assert row["latitude"] == "1.5"
    assert row["longitude"] == "-2.5"


def test_tabs_take_precedence_even_for_csv_extension_and_commas_in_values(tmp_path):
    metadata = tmp_path / "metadata.csv"
    metadata.write_text("shared_identifier\tlatitude\tlongitude\tnote\nS1\t1\t2\ta,b,c\n")
    table = read_export_metadata(metadata)
    assert table.rows_by_id["S1"]["note"] == "a,b,c"


def test_recorded_sample_id_not_directory_basename_drives_metadata_join(tmp_path):
    results = tmp_path / "results"
    path = _write_calls(results, "directory_name", {"L1": 3})
    text = path.read_text().replace("directory_name", "RECORDED_ID")
    path.write_text(text)
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("shared_identifier\tlatitude\tlongitude\nRECORDED_ID\t1\t2\n")
    output = tmp_path / "export"

    export_myoga(results, metadata, output)

    assert _read_tsv(output / "samples_used.tsv")[0]["sample_id"] == "RECORDED_ID"
    assert "'RECORDED_ID'" in (output / "mlva_nj.tree").read_text()


def test_missing_metadata_keeps_tree_sample_and_writes_blank_metadata_row(tmp_path):
    results = tmp_path / "results"
    _write_calls(results, "S1", {"L1": 1})
    _write_calls(results, "S_MISSING", {"L1": 2})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("shared_identifier\tlatitude\tlongitude\nS1\t1\t2\n")
    output = tmp_path / "export"

    export_myoga(results, metadata, output)

    rows = {row["sample_id"]: row for row in _read_tsv(output / "myoga_metadata.tsv")}
    assert set(rows) == {"S1", "S_MISSING"}
    assert rows["S_MISSING"]["latitude"] == ""
    assert any(
        row["sample_id"] == "S_MISSING" and row["reason"] == "METADATA_NOT_FOUND" and row["scope"] == "geography"
        for row in _read_tsv(output / "samples_excluded.tsv")
    )


def test_conflicting_metadata_sample_id_is_preserved_without_overwriting_tree_id(tmp_path):
    results = tmp_path / "results"
    _write_calls(results, "TREE_ID", {"L1": 1})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "shared_identifier\tsample_id\tlatitude\tlongitude\n"
        "TREE_ID\tORIGINAL_METADATA_ID\t1\t2\n"
    )
    output = tmp_path / "export"

    export_myoga(results, metadata, output)

    row = _read_tsv(output / "myoga_metadata.tsv")[0]
    assert row["sample_id"] == "TREE_ID"
    assert row["metadata_sample_id"] == "ORIGINAL_METADATA_ID"
    assert "'TREE_ID'" in (output / "mlva_nj.tree").read_text()


def test_callable_threshold_uses_exact_numeric_repeat_not_presence(tmp_path):
    results = tmp_path / "results"
    _write_calls(results, "GOOD", {"L1": 1, "L2": 2, "L3": 3, "L4": 4})
    _write_calls(results, "PRESENCE", {"L1": 1, "L2": 2, "L3": None, "L4": None})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("shared_identifier\tlatitude\tlongitude\nGOOD\t1\t2\nPRESENCE\t3\t4\n")
    output = tmp_path / "export"

    export_myoga(results, metadata, output, min_callable_fraction=0.75)

    assert [row["sample_id"] for row in _read_tsv(output / "samples_used.tsv")] == ["GOOD"]
    excluded = _read_tsv(output / "samples_excluded.tsv")
    row = next(row for row in excluded if row["sample_id"] == "PRESENCE")
    assert row["reason"] == "TOO_FEW_CALLABLE_LOCI"
    assert row["callable_loci"] == "2"


def test_min_callable_loci_and_fraction_both_apply(tmp_path):
    results = tmp_path / "results"
    _write_calls(results, "S1", {"L1": 1, "L2": 2, "L3": None, "L4": None})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("shared_identifier\tlatitude\tlongitude\nS1\t1\t2\n")
    output = tmp_path / "export"

    export_myoga(
        results,
        metadata,
        output,
        min_callable_fraction=0.25,
        min_callable_loci=3,
    )

    assert _read_tsv(output / "samples_used.tsv") == []
    assert not (output / "mlva_nj.tree").exists()
    assert _read_tsv(output / "samples_excluded.tsv")[0]["reason"] == "TOO_FEW_CALLABLE_LOCI"


def test_pairwise_missing_loci_are_shared_only_and_never_imputed_zero():
    calls = np.asarray([[1.0, np.nan, 4.0], [3.0, 99.0, np.nan]])
    rows, categorical, repeat, overlap = calculate_pairwise_distances(
        ["A", "B"], calls, min_pairwise_loci=1, min_pairwise_fraction=0
    )
    assert rows[0]["loci_compared"] == 1
    assert rows[0]["categorical_differences"] == 1
    assert rows[0]["repeat_distance_raw"] == "2.00000000"
    assert repeat[0, 1] == 2
    assert categorical[0, 1] == 1
    assert overlap[0, 1] == 1


def test_pairwise_overlap_filter_writes_na_and_prunes_poorly_connected_sample(tmp_path):
    results = tmp_path / "results"
    _write_calls(results, "A", {"L1": 1, "L2": 2, "L3": None, "L4": None})
    _write_calls(results, "B", {"L1": 2, "L2": 2, "L3": None, "L4": None})
    _write_calls(results, "C", {"L1": None, "L2": None, "L3": 3, "L4": 4})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("shared_identifier\tlatitude\tlongitude\nA\t1\t2\nB\t3\t4\nC\t5\t6\n")
    output = tmp_path / "export"

    export_myoga(
        results,
        metadata,
        output,
        min_callable_fraction=0.5,
        min_pairwise_loci=2,
        min_pairwise_fraction=0.5,
    )

    pairs = _read_tsv(output / "mlva_pairwise_distances.tsv")
    unsupported = next(row for row in pairs if {row["sample_1"], row["sample_2"]} == {"A", "C"})
    assert unsupported["comparison_status"] == "insufficient_overlap"
    assert unsupported["repeat_distance"] == ""
    assert {row["sample_id"] for row in _read_tsv(output / "samples_used.tsv")} == {"A", "B"}
    assert any(row["sample_id"] == "C" and row["reason"] == "INSUFFICIENT_PAIRWISE_OVERLAP" for row in _read_tsv(output / "samples_excluded.tsv"))


def test_categorical_distance_can_drive_matrix_and_two_sample_tree(tmp_path):
    results = tmp_path / "results"
    _write_calls(results, "A", {"L1": 1, "L2": 2})
    _write_calls(results, "B", {"L1": 9, "L2": 2})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("shared_identifier\tlatitude\tlongitude\nA\t1\t2\nB\t3\t4\n")
    output = tmp_path / "export"

    export_myoga(results, metadata, output, distance="categorical")

    matrix = _read_tsv(output / "mlva_distance_matrix.tsv")
    assert matrix[0]["B"] == matrix[1]["A"] == "0.50000000"
    tree = (output / "mlva_nj.tree").read_text()
    assert set(re.findall(r"'([^']+)'", tree)) == {"A", "B"}
    assert all(float(value) >= 0 for value in re.findall(r":([0-9.]+)", tree))


def test_one_sample_tree_is_valid_and_zero_sample_writes_no_tree(tmp_path):
    results = tmp_path / "results"
    _write_calls(results, "A", {"L1": 1})
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("shared_identifier\tlatitude\tlongitude\nA\t1\t2\n")
    one = tmp_path / "one"
    export_myoga(results, metadata, one)
    assert (one / "mlva_nj.tree").read_text() == "('A':0.00000000);\n"

    zero = tmp_path / "zero"
    export_myoga(results, metadata, zero, min_callable_loci=2)
    assert not (zero / "mlva_nj.tree").exists()
    assert _read_tsv(zero / "myoga_metadata.tsv") == []


def test_failed_and_incomplete_results_are_reported_without_crashing(tmp_path):
    results = tmp_path / "results"
    _write_calls(results, "GOOD", {"L1": 1})
    _write_calls(results, "FAILED_WITH_CALLS", {"L1": 2})
    partial = results / "PARTIAL"
    partial.mkdir(parents=True)
    (partial / "locus_repeat_counts.tsv").write_text("sample_id\tlocus_id\trepeat_count\nPARTIAL\tL1\t1\n")
    (results / "batch_status.tsv").write_text(
        "sample_id\tstatus\tmessage\n"
        "GOOD\tsuccess\t\n"
        "FAILED_WITH_CALLS\tfailed\tassembler failed\n"
        "FAILED_NO_DIR\tfailed\tmissing reads\n"
    )
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("shared_identifier\tlatitude\tlongitude\nGOOD\t1\t2\n")
    output = tmp_path / "export"

    export_myoga(results, metadata, output)

    reasons = {(row["sample_id"], row["reason"]) for row in _read_tsv(output / "samples_excluded.tsv")}
    assert ("FAILED_WITH_CALLS", "FAILED_BATCH_SAMPLE") in reasons
    assert ("FAILED_NO_DIR", "FAILED_BATCH_SAMPLE") in reasons
    assert ("PARTIAL", "MISSING_CALLS_FILE") in reasons
    assert [row["sample_id"] for row in _read_tsv(output / "samples_used.tsv")] == ["GOOD"]


def test_batch_root_combined_calls_do_not_duplicate_leaf_samples(tmp_path):
    results = tmp_path / "results"
    leaf = _write_calls(results, "S1", {"L1": 1})
    (results / "calls.tsv").write_text(leaf.read_text())
    (results / "batch_status.tsv").write_text("sample_id\tstatus\tmessage\nS1\tsuccess\t\n")
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("shared_identifier\tlatitude\tlongitude\nS1\t1\t2\n")
    output = tmp_path / "export"

    export_myoga(results, metadata, output)

    assert [row["sample_id"] for row in _read_tsv(output / "samples_used.tsv")] == ["S1"]
    assert not any(row["reason"] == "DUPLICATE_SAMPLE_ID" for row in _read_tsv(output / "samples_excluded.tsv"))


def test_batch_root_combined_calls_are_a_fallback_when_leaf_files_are_absent(tmp_path):
    staging = tmp_path / "staging"
    first = _write_calls(staging, "S1", {"L1": 1, "L2": 2})
    second = _write_calls(staging, "S2", {"L1": 2, "L2": 3})
    first_lines = first.read_text().splitlines()
    second_lines = second.read_text().splitlines()
    results = tmp_path / "results"
    results.mkdir()
    (results / "calls.tsv").write_text(
        "\n".join([*first_lines, *second_lines[1:]]) + "\n"
    )
    (results / "batch_status.tsv").write_text(
        "sample_id\tstatus\tmessage\nS1\tsuccess\t\nS2\tskipped_success\t\n"
    )
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "shared_identifier\tlatitude\tlongitude\nS1\t1\t2\nS2\t3\t4\n"
    )
    output = tmp_path / "export"

    export_myoga(results, metadata, output)

    assert [row["sample_id"] for row in _read_tsv(output / "samples_used.tsv")] == ["S1", "S2"]
    assert _read_tsv(output / "export_summary.tsv")[0] == {
        "metric": "result_directories_discovered",
        "value": "2",
    }


def test_outputs_are_deterministic_and_force_is_required_for_overwrite(tmp_path):
    results, metadata = _basic_results(tmp_path)
    first = tmp_path / "first"
    export_myoga(results, metadata, first, min_callable_fraction=2 / 3)
    names = (
        "myoga_metadata.tsv", "mlva_profiles.tsv", "mlva_calls_long.tsv",
        "mlva_pairwise_distances.tsv", "mlva_distance_matrix.tsv", "mlva_nj.tree",
        "samples_used.tsv", "samples_excluded.tsv", "export_summary.tsv", "export_summary.txt",
    )
    snapshot = {name: (first / name).read_bytes() for name in names}
    with pytest.raises(ValueError, match="--force"):
        export_myoga(results, metadata, first, min_callable_fraction=2 / 3)
    export_myoga(results, metadata, first, min_callable_fraction=2 / 3, force=True)
    assert {name: (first / name).read_bytes() for name in names} == snapshot


@pytest.mark.parametrize(
    "contents, message",
    [
        ("sample_id\tlocus_id\nS1\tL1\n", "missing required calls columns"),
        ("sample_id\tlocus_id\trepeat_count\nS1\tL1\tunknown\n", "non-numeric repeat_count"),
        ("sample_id\tlocus_id\trepeat_count\nS1\tL1\t1\nS1\tL1\t2\n", "duplicate loci"),
    ],
)
def test_malformed_calls_are_classified(tmp_path, contents, message):
    results = tmp_path / "results"
    sample = results / "S1"
    sample.mkdir(parents=True)
    (sample / "calls.tsv").write_text(contents)
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("shared_identifier\tlatitude\tlongitude\nS1\t1\t2\n")
    output = tmp_path / "export"

    export_myoga(results, metadata, output)

    row = _read_tsv(output / "samples_excluded.tsv")[0]
    assert row["reason"] == "MALFORMED_RESULTS"
    assert message in row["details"]


@pytest.mark.parametrize(
    "contents, message",
    [
        ("wrong_id\tlatitude\tlongitude\nS1\t1\t2\n", "was not found"),
        ("shared_identifier\tlatitude\tlongitude\nS1\t1\t2\nS1\t3\t4\n", "Duplicate metadata identifier"),
        ("shared_identifier latitude longitude\nS1 1 2\n", "comma- or tab-delimited"),
    ],
)
def test_malformed_metadata_has_clear_errors(tmp_path, contents, message):
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(contents)
    with pytest.raises(ValueError, match=message):
        read_export_metadata(metadata)


def test_export_myoga_cli_routes_all_options(monkeypatch, capsys):
    observed = {}

    def fake_export(results, metadata, outdir, **options):
        observed.update(
            {"results": results, "metadata": metadata, "outdir": outdir, **options}
        )
        return {
            "metadata": Path(outdir) / "myoga_metadata.tsv",
            "distance_matrix": Path(outdir) / "mlva_distance_matrix.tsv",
            "tree": Path(outdir) / "mlva_nj.tree",
            "summary": Path(outdir) / "export_summary.tsv",
            "combined_marker_tree": Path(outdir) / "combined_marker_nj.tree",
        }

    monkeypatch.setattr(cli, "export_myoga", fake_export)
    assert cli.main(
        [
            "export-myoga", "--results", "results", "--metadata", "meta.csv",
            "--metadata-id", "run", "--latitude", "lat", "--longitude", "lon",
            "--min-callable-fraction", "0.75", "--min-callable-loci", "3",
            "--min-pairwise-loci", "2", "--min-pairwise-fraction", "0.6",
            "--distance", "categorical", "--combined-markers", "--loci", "panel.tsv",
            "--phylogeny-snp-weight", "2", "--phylogeny-repeat-weight", "3",
            "-t", "4", "--mafft-bin", "mafft-x", "--raxml-ng-bin", "raxml-x",
            "--raxml-model", "GTR+G", "--force", "-o", "export",
        ]
    ) == 0
    assert observed == {
        "results": "results",
        "metadata": "meta.csv",
        "outdir": "export",
        "metadata_id": "run",
        "latitude": "lat",
        "longitude": "lon",
        "min_callable_fraction": 0.75,
        "min_callable_loci": 3,
        "min_pairwise_loci": 2,
        "min_pairwise_fraction": 0.6,
        "distance": "categorical",
        "combined_markers": True,
        "loci_path": "panel.tsv",
        "snp_weight": 2.0,
        "repeat_weight": 3.0,
        "threads": 4,
        "mafft_bin": "mafft-x",
        "raxml_ng_bin": "raxml-x",
        "raxml_model": "GTR+G",
        "force": True,
    }
    stream = capsys.readouterr().out
    assert "Wrote MLVA relatedness tree" in stream
    assert "Wrote combined SNP/repeat relatedness tree" in stream


def _write_gzip_fasta(path: Path, name: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(f">{name}\n{sequence}\n")


def _fake_export_mafft(tmp_path: Path) -> Path:
    executable = tmp_path / "mafft"
    executable.write_text(
        """#!/usr/bin/env python3
import pathlib, sys
if '--version' in sys.argv:
    print('v7.0', file=sys.stderr)
    raise SystemExit(0)
print(pathlib.Path(sys.argv[-1]).read_text().rstrip())
"""
    )
    executable.chmod(0o755)
    return executable


def _fake_export_raxml(tmp_path: Path) -> Path:
    executable = tmp_path / "raxml-ng"
    executable.write_text(
        """#!/usr/bin/env python3
import pathlib, sys
if '--version' in sys.argv:
    print('RAxML-NG v1.2.2')
    raise SystemExit(0)
args = sys.argv[1:]
prefix = args[args.index('--prefix') + 1]
msa = pathlib.Path(args[args.index('--msa') + 1])
names = [line[1:].split()[0] for line in msa.read_text().splitlines() if line.startswith('>')]
pathlib.Path(prefix + '.raxml.bestTree').write_text(
    '(' + ','.join(f'{name}:0.1' for name in names) + ');\\n'
)
pathlib.Path(prefix + '.raxml.bestModel').write_text('GTR+G\\n')
"""
    )
    executable.chmod(0o755)
    return executable


def test_combined_marker_export_reuses_precomputed_masked_queries(tmp_path):
    results = tmp_path / "results"
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "shared_identifier\tlatitude\tlongitude\n"
        "S1\t1\t2\nS2\t3\t4\nS3\t5\t6\n"
    )
    for sample_id, repeat, sequence in (
        ("S1", 1, "AAAA"),
        ("S2", 2, "AAAT"),
        ("S3", 1, "AAAA"),
    ):
        path = _write_calls(results, sample_id, {"L1": repeat})
        _write_gzip_fasta(
            path.parent / "phylogeny" / "L1" / "query.fasta.gz",
            f"QUERY__{sample_id}",
            sequence,
        )
    output = tmp_path / "export"

    result = export_myoga(
        results,
        metadata,
        output,
        combined_markers=True,
        mafft_bin=str(_fake_export_mafft(tmp_path)),
    )

    assert result["combined_marker_tree"] == output / "combined_marker_nj.tree"
    status = _read_tsv(output / "locus_tree_status.tsv")
    assert status[0]["status"] == "TWO_HAPLOTYPE_DISTANCE"
    pairs = {
        (row["sample_1"], row["sample_2"]): row
        for row in _read_tsv(output / "combined_marker_pairwise_distances.tsv")
    }
    assert pairs[("S1", "S2")]["mean_normalized_snp_distance"] == "1.00000000"
    assert pairs[("S1", "S2")]["mean_normalized_repeat_distance"] == "2.00000000"
    assert pairs[("S1", "S2")]["combined_marker_distance"] == "3.00000000"
    assert pairs[("S1", "S3")]["combined_marker_distance"] == "0.00000000"
    assert {
        row["sample_id"] for row in _read_tsv(output / "combined_marker_metadata.tsv")
    } == {"S1", "S2", "S3"}


def test_combined_marker_export_runs_raxml_for_three_snp_haplotypes(tmp_path):
    results = tmp_path / "results"
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "shared_identifier\tlatitude\tlongitude\n"
        "S1\t1\t2\nS2\t3\t4\nS3\t5\t6\n"
    )
    for sample_id, sequence in (("S1", "AAAA"), ("S2", "AAAT"), ("S3", "AATT")):
        path = _write_calls(results, sample_id, {"L1": 1})
        _write_gzip_fasta(
            path.parent / "phylogeny" / "L1" / "query.fasta.gz",
            f"QUERY__{sample_id}",
            sequence,
        )
    output = tmp_path / "export"

    export_myoga(
        results,
        metadata,
        output,
        combined_markers=True,
        mafft_bin=str(_fake_export_mafft(tmp_path)),
        raxml_ng_bin=str(_fake_export_raxml(tmp_path)),
    )

    status = _read_tsv(output / "locus_tree_status.tsv")[0]
    assert status["status"] == "RAXML_NG"
    assert status["snp_haplotypes"] == "3"
    assert (output / "locus_trees" / "L1" / "haplotypes.raxml.tree").is_file()
    tips = set(re.findall(r"'([^']+)'", (output / "locus_trees" / "L1" / "samples.tree").read_text()))
    assert tips == {"S1", "S2", "S3"}


def test_retained_assembly_evidence_is_masked_with_rich_panel(tmp_path):
    sample_dir = tmp_path / "results" / "S1"
    calls_path = _write_calls(tmp_path / "results", "S1", {"L1": 2})
    rows = _read_tsv(calls_path)
    rows[0]["evidence"] = "L1|contig|+|1-14"
    with calls_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CALL_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    _write_gzip_fasta(
        sample_dir / "assembly_amplicons.fasta.gz",
        "L1|contig|+|1-14",
        "AAAGGATATCCGGG",
    )
    sample = SampleCalls("S1", calls_path, rows, ["L1"])
    locus = Locus(
        locus_id="L1",
        forward_primer="AAA",
        reverse_primer="CCC",
        left_flank_sequence="GG",
        right_flank_sequence="CC",
        repeat_motif="AT",
        repeat_unit_length_bp=2,
    )

    recovered, statuses = recover_masked_sequences([sample], ["L1"], {"L1": locus})

    assert recovered[("S1", "L1")].sequence == "AAAGGCCGGG"
    assert recovered[("S1", "L1")].masking_method == "flank_bounded"
    assert statuses[0]["status"] == "RECOVERED"


def test_neighbor_joining_preserves_quoted_sample_ids_exactly():
    labels = ["O'Brien sample", "plain"]
    tree = neighbor_joining_tree_from_matrix(
        labels, np.asarray([[0.0, 1.0], [1.0, 0.0]])
    )
    parsed_labels, _matrix = _tip_patristic_distance_matrix(_parse_newick(tree))
    assert set(parsed_labels) == set(labels)
