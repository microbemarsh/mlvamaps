"""Comparison must preserve status changes, missing calls, and both directions."""
from pathlib import Path
import json
import shutil
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from benchmark_fastq_inference import challenge_metrics, file_sha256, version_disagreements


def test_version_audit_retains_uncertainty_and_does_not_pick_a_winner(tmp_path):
    before = {("s", "a"): {"repeat_count": 5, "status": "PASS"},
              ("s", "b"): {"repeat_count": None, "status": "UNRESOLVED"}}
    after = {("s", "a"): {"repeat_count": 5, "status": "ESTIMATED"},
             ("s", "c"): {"repeat_count": 7, "status": "PASS"}}
    result = version_disagreements(before, after)
    assert [row["locus"] for row in result] == ["a", "b", "c"]
    assert result[0]["baseline"] == before[("s", "a")]
    assert result[0]["updated"] == after[("s", "a")]
    assert version_disagreements(before, before) == []
    path = tmp_path / "reads.fastq"
    path.write_bytes(b"abc")
    assert file_sha256(path) == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"


def test_challenge_metrics_keep_missing_calls_in_denominator(tmp_path):
    output = tmp_path / "updated_0" / "calls"
    output.mkdir(parents=True)
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    (inputs / "read_truth.json").write_text(json.dumps({"a": "correct", "b": "wrong", "c": "missing"}))
    (output / "calls.tsv").write_text(
        "sample_id\tlocus_id\trepeat_count\tallele_confidence\n"
        "s\tcorrect\t5\t.99\ns\twrong\t4\t.6\ns\tmissing\t\t0\ns\tabsent\t\t0\n")
    (output / "molecule_candidate_evidence.tsv").write_text(
        "molecule\tlocus\tquality_flags\na\tcorrect\t\na\tcorrect\t\n"
        "b\tcorrect\t\nc\twrong\tambiguous_locus\n")
    truth = {name: dict(repeat_count=value, alleles=[value], coverage=3,
                       motif_length=4, repeat_bases=20, error_rate=0)
             for name, value in (("correct", 5), ("wrong", 9), ("missing", 7), ("absent", None))}
    result = challenge_metrics(output, truth)
    assert result["call_rate"] == 2 / 3
    assert result["accuracy_among_calls"] == .5
    assert result["mae"] == 2.5
    assert result["exact_accuracy_all_truth_known"] == 1 / 3
    assert result["numerical_completion_rate"] == .5
    assert result["missing_truth_known"] == 1 and result["mae_all_truth_known"] is None
    assert result["evidence_classes"]["UNREPORTED"]["total"] == 3
    assert result["calibration"]["model_support"]["bins"][-1]["n"] == 1
    assert result["false_positive_genotype_rate"] == 0
    assert result["locus_assignment_accuracy"] == .5
    assert result["locus_assignment_rate"] == 2 / 3
    assert result["accuracy_call_rate_curve"][-1]["accuracy"] == 1
    assert result["accuracy_call_rate_curve"][-1]["call_rate"] == 1 / 3


@pytest.mark.skipif(not shutil.which("minimap2") or not shutil.which("sassy"),
                    reason="minimap2 and sassy are needed for the benchmark workflows")
@pytest.mark.parametrize("baseline", [False, True])
def test_challenge_cli_runs_hard_loci_with_or_without_baseline(tmp_path, monkeypatch, baseline):
    import benchmark_fastq_inference as benchmark

    workdir, output = tmp_path / "retained", tmp_path / "report.json"
    arguments = ["benchmark_fastq_inference.py", "--challenge", "--layouts", "paired",
                 "--copies", "1", "--repetitions", "1", "--sample-mode", "isolate",
                 "--max-candidates", "32", "--max-likelihood-cells", "100000",
                 "--output", str(output)]
    if baseline:
        arguments += ["--baseline-source", str(Path(__file__).resolve().parents[1]), "--workdir", str(workdir)]
    else:
        # The ordinary --challenge invocation retains an automatic workdir.
        original_mkdtemp = benchmark.tempfile.mkdtemp
        monkeypatch.setattr(benchmark.tempfile, "mkdtemp", lambda *args, **kwargs:
            str(workdir) if kwargs.get("prefix") == "mlvamaps-benchmark-" else original_mkdtemp(*args, **kwargs))
    monkeypatch.setattr(sys, "argv", arguments)
    benchmark.main()

    report = json.loads(output.read_text())
    retained = workdir / ("comparison.json" if baseline else "benchmark.json")
    assert output.read_bytes() == retained.read_bytes()
    assert set(report["sources"]) == ({"baseline", "updated"} if baseline else {"updated"})
    assert set(report["workflows"]) == {"paired", "challenge_paired"}
    for layout, workflow in report["workflows"].items():
        trial = workflow["trials"][0]
        assert ("comparison" in trial) == baseline
        for version in report["sources"]:
            measurement = trial[version]
            assert measurement["exit_code"] == 0
            assert measurement["seconds"] > 0 and measurement["peak_rss_mib"] > 0
            assert "competitive_mapping_seconds" in measurement["stage_metrics"]
            assert measurement["artifact_sha256"]["calls.tsv"]
            argv = json.loads((workdir / layout / f"{version}_0" / "arguments.json").read_text())
            assert argv[argv.index("--fastq-max-candidates") + 1] == "32"
            assert argv[argv.index("--fastq-max-likelihood-cells") + 1] == "100000"
            if layout == "challenge_paired":
                assert measurement["challenge_metrics"]["truth_known_loci"] == 14
                assert measurement["challenge_metrics"]["configured_loci"] == 15
                assert (Path(measurement["output"]) / "genotype_evidence.tsv").exists()


def test_benchmark_rejects_nonpositive_limits(monkeypatch):
    import benchmark_fastq_inference as benchmark

    monkeypatch.setattr(sys, "argv", ["benchmark_fastq_inference.py", "--max-candidates", "0"])
    with pytest.raises(SystemExit) as error:
        benchmark.main()
    assert error.value.code == 2


@pytest.mark.parametrize("relative", ["paired/inputs/r1.fastq", "paired/updated_0/calls/calls.tsv",
                                     "paired/updated_0/measurement.json", "challenge_paired/inputs/truth.json", "paired"])
def test_report_cannot_overwrite_retained_inputs_or_artifacts(tmp_path, relative):
    import argparse
    import benchmark_fastq_inference as benchmark

    workdir = tmp_path / "retained"
    args = argparse.Namespace(workdir=workdir, output=workdir / relative)
    with pytest.raises(ValueError, match="top-level JSON file"):
        benchmark.compare_versions(args, {})
    assert not workdir.exists()


def test_report_symlink_cannot_overwrite_retained_artifacts(tmp_path):
    import argparse
    import benchmark_fastq_inference as benchmark

    workdir, output = tmp_path / "retained", tmp_path / "report.json"
    output.symlink_to(workdir / "paired/updated_0/measurement.json")
    with pytest.raises(ValueError, match="top-level JSON file"):
        benchmark.compare_versions(argparse.Namespace(workdir=workdir, output=output), {})
    assert not workdir.exists()


@pytest.mark.parametrize("exit_code,hashes,repeatable", [
    (0, [{"calls.tsv": "a"}, {"calls.tsv": "a"}], True),
    (0, [{"calls.tsv": "a"}, {"calls.tsv": "b"}], False),
    (0, [{}, {}], False), (1, [None, None], False),
    (1, [{"calls.tsv": "a"}, {"calls.tsv": "a"}], False),
])
def test_repeatability_requires_successful_workers_and_nonempty_matching_hashes(
        tmp_path, monkeypatch, exit_code, hashes, repeatable):
    import argparse
    import benchmark_fastq_inference as benchmark

    calls = "sample_id\tlocus_id\trepeat_count\ns\tL\t5\n"
    def fixture(directory, layout, tolerance, patch):
        reads = directory / "r1.fastq"
        reads.write_text("@m\nACGT\n+\nIIII\n")
        (directory / "assembly").mkdir()
        (directory / "assembly/calls.tsv").write_text(calls)
        main(["call", "--sample-mode", "isolate", "--fq1", str(reads), "-o", str(directory / "prepared")])

    monkeypatch.setitem(fixture.__globals__, "main", lambda argv: 0)
    remaining = iter(hashes)
    def worker(command, **kwargs):
        arguments = json.loads(Path(command[4]).read_text())
        output = Path(arguments[arguments.index("-o") + 1])
        output.mkdir()
        (output / "calls.tsv").write_text(calls)
        result = {"exit_code": exit_code, "seconds": 1.}
        artifact_hashes = next(remaining)
        if artifact_hashes is not None:
            result["artifact_sha256"] = artifact_hashes
        Path(command[5]).write_text(json.dumps(result))
        return argparse.Namespace(returncode=exit_code)

    monkeypatch.setattr(benchmark.subprocess, "run", worker)
    workdir = tmp_path / "retained"
    args = argparse.Namespace(workdir=workdir, output=workdir / "custom.json", baseline_source=None,
        copies=1, repetitions=2, sample_mode="isolate", layouts=["paired"], challenge=False,
        max_candidates=None, max_likelihood_cells=None, ablations=[])
    if exit_code:
        with pytest.raises(SystemExit) as error:
            benchmark.compare_versions(args, {"test_short_read_fingerprint_matches_assembly": fixture})
        assert error.value.code == 1
    else:
        benchmark.compare_versions(args, {"test_short_read_fingerprint_matches_assembly": fixture})
    report = json.loads(args.output.read_text())
    assert args.output.read_bytes() == (workdir / "benchmark.json").read_bytes()
    assert report["workflows"]["paired"]["deterministic_artifacts"]["updated"] is repeatable
    assert all(trial["updated"]["exit_code"] == exit_code for trial in report["workflows"]["paired"]["trials"])
