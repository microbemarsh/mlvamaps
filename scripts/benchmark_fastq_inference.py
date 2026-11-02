"""Time real FASTQ workflows against the assembly-concordance fixture.

Run from a checkout with dev dependencies, minimap2 and sassy on PATH:
python scripts/benchmark_fastq_inference.py --repetitions 3 --copies 25

Only FASTQ calling is timed; database construction and assembly truth are not.
The fixture asserts the final fingerprint, including novel/fractional alleles.
"""

import argparse
import contextlib
import csv
import io
import json
import hashlib
import os
from pathlib import Path
import runpy
import statistics
import sys
import subprocess
import tempfile
import time

import pytest

from compare_short_read_results import audit_calls, read_calls


def file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def version_disagreements(baseline, updated):
    fields = ("repeat_count", "status", "estimation_method")
    return [{"sample": key[0], "locus": key[1], "baseline": baseline.get(key),
             "updated": updated.get(key)} for key in sorted(baseline.keys() | updated.keys())
            if key not in baseline or key not in updated
            or any(baseline[key].get(field) != updated[key].get(field) for field in fields)]


def version_worker(source, arguments, output, ablation=""):
    """Fresh process: import exactly one checkout, retaining failures and timings."""
    import resource
    sys.path.insert(0, str(Path(source).resolve()))
    from mlvamaps.cli import main as call
    if ablation == "no_error_model":
        import mlvamaps.long_read_evidence as module
        module.flank_error_model = lambda *args: None
    elif ablation == "no_motif_refinement":
        import mlvamaps.long_read_evidence as module
        module.refine_repeat_sequence = lambda *args: ({}, "")
    elif ablation == "no_partial_bounds":
        import mlvamaps.short_read_evidence as module
        module.partial_repeat_lower_bound = lambda *args: None

    started = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_SELF)
    children_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = {}
    try:
        result["exit_code"] = call(arguments)
    except Exception as exc:
        result.update(exit_code=1, error=f"{type(exc).__name__}: {exc}")
    after = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    result.update(seconds=time.perf_counter() - started,
                  cpu_seconds=after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime
                  + children.ru_utime + children.ru_stime - children_before.ru_utime - children_before.ru_stime,
                  peak_rss_mib=max(after.ru_maxrss, children.ru_maxrss)
                  / (1024 ** 2 if sys.platform == "darwin" else 1024))
    root = Path(arguments[arguments.index("-o") + 1])
    performance = root / "fastq_performance.tsv"
    if performance.exists():
        with performance.open() as handle:
            result["stage_metrics"] = {row["metric"]: float(row["value"]) for row in csv.DictReader(handle, delimiter="\t")}
    result["artifact_sha256"] = {name: file_sha256(root / name) for name in (
        "calls.tsv", "common_locus_calls.tsv", "genotype_evidence.tsv", "sequence_alleles.tsv",
        "mlva_fingerprint.tsv", "locus_snps.tsv", "molecule_candidate_evidence.tsv") if (root / name).exists()}
    Path(output).write_text(json.dumps(result, indent=2) + "\n")


def challenge_metrics(output, truth):
    calls = read_calls(output, with_metadata=True)
    values = {locus: row["repeat_count"] for (_, locus), row in calls.items()}
    present = {locus: row for locus, row in truth.items() if row["repeat_count"] is not None}
    called = {locus: value for locus, value in values.items() if value is not None and locus in present}
    errors = [min(abs(value - allele) for allele in present[locus]["alleles"]) for locus, value in called.items()]
    absent = [locus for locus, row in truth.items() if row["repeat_count"] is None]
    evidence_path = Path(output) / "genotype_evidence.tsv"
    evidence = {}
    if evidence_path.exists():
        with evidence_path.open() as handle:
            evidence = {row["locus"]: row for row in csv.DictReader(handle, delimiter="\t")}
    evidence_classes, calibration = {}, {}
    for name, expected in present.items():
        row = evidence.get(name, {})
        group = row.get("genotype_status", "UNREPORTED")
        bucket = evidence_classes.setdefault(group, dict(total=0, numeric=0, correct=0, absolute_errors=[]))
        bucket["total"] += 1
        if name in called:
            error = min(abs(called[name] - allele) for allele in expected["alleles"])
            bucket["numeric"] += 1
            bucket["correct"] += error == 0
            bucket["absolute_errors"].append(error)
    for bucket in evidence_classes.values():
        bucket["accuracy_all_loci"] = bucket["correct"] / bucket["total"]
        bucket["mae_among_numeric"] = statistics.mean(bucket["absolute_errors"]) if bucket["absolute_errors"] else None
    for field in ("model_support", "data_evidence_confidence"):
        points = [(float(evidence.get(name, {}).get(field) or row.get(field) or
                         (row.get("allele_confidence") if field == "model_support" else 0) or 0),
                   int(called[name] in present[name]["alleles"]))
                  for (_, name), row in calls.items() if name in called
                  and (field == "model_support" or field in row or field in evidence.get(name, {}))]
        bins = []
        for lower, upper in ((0, .25), (.25, .5), (.5, .75), (.75, 1.000001)):
            selected = [(support, correct) for support, correct in points if lower <= support < upper]
            bins.append(dict(lower=lower, upper=min(1, upper), n=len(selected),
                mean_support=statistics.mean(support for support, _ in selected) if selected else None,
                observed_accuracy=statistics.mean(correct for _, correct in selected) if selected else None))
        calibration[field] = dict(bins=bins, brier=statistics.mean((support - correct)**2 for support, correct in points) if points else None)
    strata = {}
    for field in ("coverage", "motif_length", "repeat_bases", "error_rate"):
        strata[field] = {}
        for value in sorted({row[field] for row in present.values()}):
            names = [name for name, row in present.items() if row[field] == value]
            observed = [name for name in names if name in called]
            strata[field][str(value)] = {"total": len(names), "called": len(observed),
                "correct": sum(called[name] in present[name]["alleles"] for name in observed)}
    sequence_path = Path(output) / "sequence_alleles.tsv"
    sequence_accuracy = None
    if sequence_path.exists():
        with sequence_path.open() as handle:
            sequences = list(csv.DictReader(handle, delimiter="\t"))
        sequence_accuracy = (sum(row["product_sequence"] in truth[row["locus"]]["products"] for row in sequences)
                             / len(sequences) if sequences else None)
    read_truth_path = Path(output).parent.parent / "inputs/read_truth.json"
    assignment_accuracy, assignment_rate = None, None
    evidence_path = Path(output) / "molecule_candidate_evidence.tsv"
    if read_truth_path.exists() and evidence_path.exists():
        read_truth = json.loads(read_truth_path.read_text())
        with evidence_path.open() as handle:
            assignments = {(row["molecule"], row["locus"]) for row in csv.DictReader(handle, delimiter="\t")
                           if "ambiguous_locus" not in row.get("quality_flags", "")}
        assignment_accuracy = sum(read_truth.get(molecule) == locus for molecule, locus in assignments) / len(assignments) if assignments else None
        assignment_rate = len({molecule for molecule, locus in assignments}) / max(1, len(read_truth))
    curve = []
    for threshold in (0, .5, .8, .95, .99):
        selected = [row for (_, locus), row in calls.items() if locus in called
                    and float(row.get("allele_confidence") or 0) >= threshold]
        curve.append({"model_support_threshold": threshold, "call_rate": len(selected) / max(1, len(present)),
            "accuracy": sum(row["repeat_count"] in present[row["locus_id"]]["alleles"] for row in selected) / len(selected) if selected else None})
    return {"call_rate": len(called) / max(1, len(present)),
        "configured_loci": len(truth), "truth_known_loci": len(present),
        "numerical_completion_rate": sum(values.get(name) is not None for name in truth) / max(1, len(truth)),
        "exact_accuracy_all_truth_known": sum(error == 0 for error in errors) / max(1, len(present)),
        "missing_truth_known": len(present) - len(called),
        "mae_all_truth_known": statistics.mean(errors) if len(called) == len(present) and errors else None,
        "evidence_classes": evidence_classes, "calibration": calibration,
        "accuracy_among_calls": sum(error == 0 for error in errors) / len(errors) if errors else None,
        "mae": statistics.mean(errors) if errors else None,
        "absent_numeric_prediction_rate": sum(values.get(name) is not None for name in absent) / max(1, len(absent)),
        "false_positive_genotype_rate": sum(row.get("present") == "yes" for (_, name), row in calls.items() if name in absent) / max(1, len(absent)),
        "sequence_accuracy_among_reported": sequence_accuracy, "strata": strata,
        "locus_assignment_accuracy": assignment_accuracy, "locus_assignment_rate": assignment_rate,
        "accuracy_call_rate_curve": curve,
        "note": "Mixture counts are scored against either known allele. Missing calls fail all-locus exact accuracy; their absolute error is undefined. "
                "Absent loci have no truth allele: numerical completion is separate from observed presence. Calibration is descriptive on this small synthetic set, not validation."}


def compare_versions(args, fixture):
    """Run both sources on the same retained files; agreement is not truth."""
    directory = args.workdir.resolve()
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("Comparison workdir must be empty; previous evidence is never overwritten")
    directory.mkdir(parents=True, exist_ok=True)
    sources = {"baseline": str(Path(args.baseline_source).resolve()),
               "updated": str(Path(__file__).resolve().parents[1])}
    report = {"sources": sources, "copies": args.copies, "sample_mode": args.sample_mode,
              "workflows": {}, "method": "Alternating serial fresh processes, identical input bytes. "
              "Wall time excludes imports. RSS is max(parent, largest child), not concurrent total. "
              "Synthetic expected alleles are specified independently of caller agreement. "
              "Preparation is excluded from timing. All outputs and disagreements are retained."}
    report["source_files_sha256"] = {version: {str(path.relative_to(source)): file_sha256(path)
        for path in sorted((Path(source) / "mlvamaps").glob("*.py"))} for version, source in sources.items()}
    test = fixture["test_short_read_fingerprint_matches_assembly"]

    def run(source, argv, target, ablation=""):
        target.mkdir(parents=True)
        arguments = list(argv)
        arguments[arguments.index("-o") + 1] = str(target / "calls")
        arguments.extend(["--keep-intermediates", "--debug-disagreements"])
        if args.max_candidates:
            arguments.extend(["--fastq-max-candidates", str(args.max_candidates)])
        if args.max_likelihood_cells:
            arguments.extend(["--fastq-max-likelihood-cells", str(args.max_likelihood_cells)])
        spec, measurement = target / "arguments.json", target / "measurement.json"
        spec.write_text(json.dumps(arguments))
        environment = dict(os.environ, PYTHONPATH=source)
        with (target / "stdout.txt").open("w") as stdout, (target / "stderr.txt").open("w") as stderr:
            completed = subprocess.run(
                [sys.executable, str(Path(__file__).resolve()), "--worker", source, str(spec), str(measurement), ablation],
                cwd=source, env=environment, stdout=stdout, stderr=stderr,
            )
        result = json.loads(measurement.read_text()) if measurement.exists() else {"exit_code": completed.returncode}
        result["output"] = str(target / "calls")
        return result

    for layout in args.layouts + (["challenge_paired", "challenge_long"] if args.challenge else []):
        prepared = directory / layout / "inputs"
        if not layout.startswith("challenge_"):
            prepared.mkdir(parents=True)
        captured = []
        original = test.__globals__["main"]

        def prepare(argv):
            argv[argv.index("--sample-mode") + 1] = args.sample_mode
            paths = [argv[argv.index("--fq1") + 1] if "--fq1" in argv else argv[argv.index("-i") + 1]]
            if "--fq2" in argv:
                paths.append(argv[argv.index("--fq2") + 1])
            for name in paths:
                path = Path(name)
                lines = path.read_text().splitlines()
                with path.open("w") as handle:
                    for copy in range(args.copies):
                        for start in range(0, len(lines), 4):
                            header, sequence, plus, quality = lines[start:start + 4]
                            handle.write(f"@copy{copy}_{header[1:]}\n{sequence}\n{plus}\n{quality}\n")
            captured.append((list(argv), paths))
            return original(argv)

        synthetic_truth = None
        if layout.startswith("challenge_"):
            from fastq_challenge_fixture import prepare_challenge
            argv, paths, synthetic_truth = prepare_challenge(prepared, layout.removeprefix("challenge_"), args.copies)
            argv[argv.index("--sample-mode") + 1] = args.sample_mode
            truth = {("sample", name): {"repeat_count": row["repeat_count"], "status": "KNOWN_SYNTHETIC_TRUTH"}
                     for name, row in synthetic_truth.items()}
        else:
            test.__globals__["main"] = prepare
            try:
                with pytest.MonkeyPatch.context() as patch, contextlib.redirect_stdout(io.StringIO()):
                    test(prepared, layout, .25, patch)
            finally:
                test.__globals__["main"] = original
            argv, paths = captured[0]
            truth = read_calls(prepared / "assembly", with_metadata=True)
        trials = []
        with Path(paths[0]).open() as handle:
            input_fragments = sum(1 for _ in handle) // 4
        for trial in range(args.repetitions):
            measurements = {}
            order = ("baseline", "updated") if trial % 2 == 0 else ("updated", "baseline")
            for version in order:
                result = run(sources[version], argv, directory / layout / f"{version}_{trial}")
                result["input_fragments_per_second"] = input_fragments / result["seconds"] if result.get("seconds") else None
                if result["exit_code"] == 0:
                    result["truth_audit"] = audit_calls(truth, read_calls(result["output"], with_metadata=True))
                    if synthetic_truth:
                        result["challenge_metrics"] = challenge_metrics(result["output"], synthetic_truth)
                measurements[version] = result
            for ablation in args.ablations:
                result = run(sources["updated"], argv, directory / layout / f"{ablation}_{trial}", ablation)
                if result["exit_code"] == 0:
                    result["truth_audit"] = audit_calls(truth, read_calls(result["output"], with_metadata=True))
                    if synthetic_truth:
                        result["challenge_metrics"] = challenge_metrics(result["output"], synthetic_truth)
                measurements[ablation] = result
            if all(result["exit_code"] == 0 for result in measurements.values()):
                before = read_calls(measurements["baseline"]["output"], with_metadata=True)
                after = read_calls(measurements["updated"]["output"], with_metadata=True)
                measurements["comparison"] = audit_calls(
                    before, after)
                measurements["comparison"]["version_disagreements"] = version_disagreements(before, after)
            trials.append(measurements)
        report["workflows"][layout] = {"inputs_sha256": {
            name: file_sha256(name) for name in paths},
            "trials": trials, "input_fragments": input_fragments,
            "deterministic_artifacts": {version: all(
                trial[version].get("artifact_sha256") == trials[0][version].get("artifact_sha256") for trial in trials)
                for version in (*sources, *args.ablations)}}
        (directory / "comparison.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(str(directory / "comparison.json"))
    if any(trial[version]["exit_code"] for workflow in report["workflows"].values()
           for trial in workflow["trials"] for version in (*sources, *args.ablations)):
        raise SystemExit(1)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        version_worker(sys.argv[2], json.loads(Path(sys.argv[3]).read_text()), sys.argv[4], sys.argv[5] if len(sys.argv) > 5 else "")
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--copies", type=int, default=25)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--baseline-source", type=Path)
    parser.add_argument("--workdir", type=Path, help="Retain identical FASTQs, both versions' outputs and disagreements")
    parser.add_argument("--challenge", action="store_true", help="Include deterministic difficult synthetic loci")
    parser.add_argument("--ablations", nargs="*", default=[],
                        choices=("no_error_model", "no_motif_refinement", "no_partial_bounds"))
    parser.add_argument("--layouts", nargs="+", choices=("paired", "shotgun", "long"),
                        default=["paired", "shotgun", "long"])
    parser.add_argument("--sample-mode", choices=("isolate", "metagenome"), default="isolate")
    parser.add_argument("--max-candidates", type=int, help="Apply the same candidate budget to both sources")
    parser.add_argument("--max-likelihood-cells", type=int, help="Apply the same likelihood cell budget to both sources")
    args = parser.parse_args()
    if args.repetitions < 1 or args.copies < 1:
        parser.error("repetitions and copies must be positive")
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root))
    fixture = runpy.run_path(str(root / "tests/test_fastq_assembly_concordance.py"))
    if args.baseline_source:
        if args.workdir is None:
            parser.error("--baseline-source requires --workdir")
        compare_versions(args, fixture)
        return
    test = fixture["test_short_read_fingerprint_matches_assembly"]
    original = test.__globals__["main"]
    timings = {}
    for layout in args.layouts:
        elapsed = []

        def timed_call(argv):
            argv[argv.index("--sample-mode") + 1] = args.sample_mode
            first = argv[argv.index("--fq1") + 1] if "--fq1" in argv else argv[argv.index("-i") + 1]
            paths = [first]
            if "--fq2" in argv:
                paths.append(argv[argv.index("--fq2") + 1])
            for name in paths:
                path = Path(name)
                lines = path.read_text().splitlines()
                with path.open("w") as handle:
                    for copy in range(args.copies):
                        for start in range(0, len(lines), 4):
                            header, sequence, plus, quality = lines[start:start + 4]
                            handle.write(f"@copy{copy}_{header[1:]}\n{sequence}\n{plus}\n{quality}\n")
            started = time.perf_counter()
            result = original(argv)
            elapsed.append(time.perf_counter() - started)
            return result

        test.__globals__["main"] = timed_call
        for _ in range(args.repetitions):
            with tempfile.TemporaryDirectory(prefix="mlvamaps-benchmark-") as directory:
                with pytest.MonkeyPatch.context() as patch, contextlib.redirect_stdout(io.StringIO()):
                    test(Path(directory), layout, .25, patch)
        timings[layout] = {"seconds": elapsed, "median_seconds": statistics.median(elapsed)}
    report = {"copies": args.copies, "repetitions": args.repetitions,
              "sample_mode": args.sample_mode,
              "threads": 2, "fingerprints_match_assembly": True, "workflows": timings}
    text = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
