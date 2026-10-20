import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location("compare_short_reads", Path(__file__).resolve().parents[1] / "scripts/compare_short_read_results.py")
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def test_metrics_include_ties_zero_counts_and_missing_calls():
    assembly = {("a", "L"): 0, ("b", "L"): 1, ("c", "L"): 2, ("d", "L"): 4, ("e", "L"): 5}
    short = {("a", "L"): 0, ("b", "L"): 1, ("c", "L"): 1, ("d", "L"): 7, ("e", "L"): None}
    result = comparison.compare(assembly, short)
    overall = result["overall"]
    assert overall["n"] == 4
    assert overall["missing_short_read"] == 1
    assert overall["exact_agreement"] == .5
    assert overall["within_1"] == overall["within_2"] == .75
    assert overall["mae"] == 1
    assert overall["median_absolute_error"] == .5
    assert overall["spearman"] == pytest.approx(.948683298)
    assert result["per_locus"]["L"] == overall
    assert comparison.metrics([])["spearman"] is None
    assert comparison.metrics([(1, 1), (2, 1)])["spearman"] is None


def test_duplicate_calls_are_rejected_and_missing_is_not_zero(tmp_path):
    calls = tmp_path / "calls.tsv"
    calls.write_text("sample_id\tlocus_id\trepeat_count\ns\tL\t\n")
    assert comparison.read_calls(calls) == {("s", "L"): None}
    calls.write_text("sample_id\tlocus_id\trepeat_count\ns\tL\t0\ns\tL\t1\n")
    with pytest.raises(ValueError, match="Duplicate"):
        comparison.read_calls(calls)


def test_result_tree_does_not_count_aggregate_rows_twice(tmp_path):
    text = "sample_id\tlocus_id\trepeat_count\ns\tL\t5\n"
    for folder in ("s", "batch_summary"):
        path = tmp_path / folder
        path.mkdir()
        (path / "calls.tsv").write_text(text)
    assert comparison.read_calls(tmp_path) == {("s", "L"): 5}
    assert comparison.read_calls(tmp_path / "batch_summary") == {("s", "L"): 5}


def test_audit_separates_estimates_and_reports_sample_mismatches(tmp_path):
    header = "sample_id\tlocus_id\trepeat_count\tstatus\testimation_method\n"
    assembly = tmp_path / "assembly.tsv"
    reads = tmp_path / "reads.tsv"
    assembly.write_text(header + "s\tL1\t0\tPASS\t\ns\tL2\t5\tPASS\t\ns\tL3\t8\tPASS\t\nonly_assembly\tL1\t3\tPASS\t\n")
    reads.write_text(header + "s\tL1\t0\tPASS\t\ns\tL2\t3\tESTIMATED\treference_conditioned\ns\tL3\t\tNOT_FOUND\t\nonly_reads\tL1\t3\tPASS\t\n")
    result = comparison.audit_calls(comparison.read_calls(assembly, with_metadata=True),
                                    comparison.read_calls(reads, with_metadata=True))
    assert result["overall"]["n"] == 2
    assert result["overall"]["mean_signed_error"] == -1
    assert result["by_read_status"]["PASS"]["overall"]["exact_agreement"] == 1
    assert result["by_estimation_method"]["reference_conditioned"]["overall"]["mae"] == 2
    assert result["sample_matching"] == {"shared": ["s"], "assembly_only": ["only_assembly"], "reads_only": ["only_reads"]}
    assert result["per_sample"]["s"]["missing_short_read"] == 1
    assert result["disagreements"][0]["locus_id"] == "L2"
    assert result["disagreements"][0]["signed_error"] == -2
    assert result["disagreements"][0]["call_file"] == str(reads)
    assert all(row["locus_id"] != "L1" or row["sample_id"] != "s" for row in result["disagreements"])


def test_read_alias_produces_a_long_read_audit(tmp_path):
    import json
    import subprocess
    import sys

    calls = tmp_path / "calls.tsv"
    calls.write_text("sample_id\tlocus_id\trepeat_count\tstatus\ns\tL\t0\tPASS\n")
    output = subprocess.check_output([sys.executable, comparison.__file__, "compare", "--assembly", str(calls),
                                      "--reads", str(calls)], text=True)
    assert json.loads(output)["updated"]["by_read_status"]["PASS"]["overall"]["exact_agreement"] == 1
