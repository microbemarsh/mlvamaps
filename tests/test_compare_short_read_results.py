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
    assert overall["exact_agreement_all_assembly_numeric"] == .4
    assert overall["greater_than_2"] == .25
    assert overall["max_absolute_error"] == 3
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


def test_audit_joins_optional_diagnostics_by_sample_and_locus_without_changing_calls(tmp_path):
    header = "sample_id\tlocus_id\trepeat_count\tstatus\tpresent\testimation_method\n"
    assembly = tmp_path / "assembly.tsv"
    calls = tmp_path / "calls.tsv"
    assembly.write_text(header + "s\tL\t40\tPASS\tyes\t\ns\tM\t\tNOT_FOUND\tno\t\ns\tN\t4\tPASS\tyes\t\n")
    calls.write_text(header + "s\tL\t17.5\tESTIMATED\tyes\tpartial_repeat_likelihood\ns\tM\t16\tESTIMATED\tyes\treference_conditioned\ns\tN\t4\tNOT_FOUND\tno\tobserved_reference_prior\n")
    (tmp_path / "genotype_evidence.tsv").write_text(
        "sample\tlocus\trepeat_count\tgenotype_status\tfull_span_support\tprior_only\n"
        "s\tL\t17.5\tLOWER_BOUND\t0\tFalse\ns\tN\t4\tPRIOR_ONLY\t0\tTrue\n"
        "other\tL\t100\tEXACT\t100\tFalse\n")
    (tmp_path / "short_read_diagnostics.tsv").write_text(
        "sample_id\tlocus_id\tspanning_fragment_count\tmicroassembly_reason\n"
        "s\tL\t0\trepeat_cycle\n")
    (tmp_path / "query_locus_reconstruction.tsv").write_text(
        "sample_id\tlocus_id\tstatus\tmethod\n"
        "s\tL\tNO_SUPPORTED_PRODUCT\t\n")
    (tmp_path / "library_insert_distribution.tsv").write_text(
        "sample\tsource\tmean_bp\tsd_bp\nother\tuser_supplied\t400\t30\ns\tunavailable\t\t\n")
    (tmp_path / "sample_summary.tsv").write_text("sample_id\tsample_mode\ns\tisolate\n")
    reads = comparison.read_calls(calls, with_metadata=True)
    assert reads[("s", "L")]["repeat_count"] == 17.5
    assert reads[("s", "L")]["status"] == "ESTIMATED"
    assert reads[("s", "L")]["diagnostics"]["genotype_evidence"]["genotype_status"] == "LOWER_BOUND"
    assert "genotype_evidence" not in reads[("s", "M")]["diagnostics"]
    assert reads[("s", "M")]["diagnostics"]["insert_calibration"]["source"] == "unavailable"
    result = comparison.audit_calls(comparison.read_calls(assembly, with_metadata=True), reads)
    assert result["overall"]["exact_agreement_all_assembly_numeric"] == .5
    assert result["presence"] == {"compared": 3, "unknown": 0, "agreements": 1, "disagreements": 2}
    mismatches = {row["locus_id"]: row for row in result["disagreements"]}
    assert mismatches["M"]["missing_side"] == "assembly"
    assert mismatches["M"]["presence_comparison"] == "disagreement"
    assert mismatches["N"]["comparison"] == "presence_disagreement"
    assert mismatches["N"]["signed_error"] == 0
    assert "microassembly:repeat_cycle" in mismatches["L"]["evidence_observations"]
    assert mismatches["L"]["diagnostics"]["run_summary"]["sample_mode"] == "isolate"
    assert result["by_evidence_observation"]["genotype:LOWER_BOUND"]["overall"]["mae"] == 22.5
    assert result["by_evidence_observation"]["insert_calibration:unavailable"]["overall"]["total_sample_loci"] == 3


def test_optional_diagnostic_duplicates_are_rejected_only_when_loading_metadata(tmp_path):
    calls = tmp_path / "calls.tsv"
    calls.write_text("sample_id\tlocus_id\trepeat_count\ns\tL\t1\n")
    (tmp_path / "genotype_evidence.tsv").write_text("sample\tlocus\ns\tL\ns\tL\n")
    assert comparison.read_calls(calls) == {("s", "L"): 1}
    with pytest.raises(ValueError, match="Duplicate diagnostic"):
        comparison.read_calls(calls, with_metadata=True)


def test_benchmark_passes_identical_explicit_library_settings_to_both_sources(tmp_path, monkeypatch):
    import argparse
    import json

    commands = []
    def run(command, **kwargs):
        commands.append(command)
        return argparse.Namespace(stdout=json.dumps({"seconds": 1., "peak_rss_mib": 10}))
    monkeypatch.setattr(comparison.subprocess, "run", run)
    args = argparse.Namespace(outdir=str(tmp_path / "benchmark"), trials=2, baseline_source="old", source="new",
        reads1="R1.fastq", reads2="R2.fastq", threads=2, sample_id="s", panel="panel.tsv", database="db",
        assembly=None, sample_mode="metagenome", assembly_round_tolerance=0., short_insert_mean=400., short_insert_sd=30.)
    result = comparison.benchmark(args)
    assert [command[command.index("--source") + 1] for command in commands] == ["old", "new", "new", "old"]
    for command in commands:
        for option, expected in (("--sample-mode", "metagenome"), ("--assembly-round-tolerance", "0.0"),
                                 ("--short-insert-mean", "400.0"), ("--short-insert-sd", "30.0")):
            assert command[command.index(option) + 1] == expected
    assert result["parameters"]["sample_mode"] == "metagenome"


@pytest.mark.parametrize("options", [["--short-insert-mean", "400"],
    ["--short-insert-mean", "inf", "--short-insert-sd", "30"], ["--assembly-round-tolerance", "nan"]])
def test_benchmark_rejects_invalid_library_settings_before_running(tmp_path, options):
    import subprocess
    import sys

    result = subprocess.run([sys.executable, comparison.__file__, "benchmark", "--reads1", "r1", "--reads2", "r2",
        "--panel", "panel", "--baseline-source", "old", "--outdir", str(tmp_path / "out"), *options], capture_output=True, text=True)
    assert result.returncode == 2
    assert not (tmp_path / "out").exists()
