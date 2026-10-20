#!/usr/bin/env python3
"""Audit matched FASTQ/assembly calls, or time two checkouts on paired FASTQs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys


def read_calls(path, *, with_metadata=False):
    path = Path(path)
    files = [path] if path.is_file() else sorted(
        file for file in path.rglob("calls.tsv") if "batch_summary" not in file.relative_to(path).parts
    )
    if not files:
        raise ValueError(f"No calls.tsv files found in {path}")
    calls = {}
    for file in files:
        with file.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = (row["sample_id"], row["locus_id"])
                if key in calls:
                    raise ValueError(f"Duplicate sample/locus {key} in {file}")
                value = row.get("repeat_count", "")
                count = float(value) if value not in ("", "NA", "N/A") else None
                if count is not None and not math.isfinite(count):
                    raise ValueError(f"Non-finite count for {key}")
                calls[key] = {**row, "repeat_count": count, "call_file": str(file.resolve())} if with_metadata else count
    return calls


def _ranks(values):
    ordered = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and values[ordered[end]] == values[ordered[start]]:
            end += 1
        for index in ordered[start:end]:
            ranks[index] = (start + end - 1) / 2
        start = end
    return ranks


def metrics(pairs):
    errors = [abs(a - b) for a, b in pairs]
    correlation = None
    if len(pairs) > 1:
        first, second = (_ranks(list(values)) for values in zip(*pairs))
        mean_first, mean_second = statistics.mean(first), statistics.mean(second)
        first = [value - mean_first for value in first]
        second = [value - mean_second for value in second]
        denominator = math.sqrt(sum(value * value for value in first) * sum(value * value for value in second))
        if denominator:
            correlation = sum(a * b for a, b in zip(first, second)) / denominator
    return {
        "n": len(pairs),
        "exact_agreement": sum(error == 0 for error in errors) / len(errors) if errors else None,
        "within_1": sum(error <= 1 for error in errors) / len(errors) if errors else None,
        "within_2": sum(error <= 2 for error in errors) / len(errors) if errors else None,
        "mae": statistics.mean(errors) if errors else None,
        "median_absolute_error": statistics.median(errors) if errors else None,
        "mean_signed_error": statistics.mean(b - a for a, b in pairs) if pairs else None,
        "spearman": correlation,
    }


def compare(assembly, short_reads):
    keys = set(assembly) | set(short_reads)

    def summarize(selected):
        pairs = [(assembly[key], short_reads[key]) for key in selected
                 if assembly.get(key) is not None and short_reads.get(key) is not None]
        return {**metrics(pairs), "total_sample_loci": len(selected),
                "missing_assembly": sum(assembly.get(key) is None for key in selected),
                "missing_short_read": sum(short_reads.get(key) is None for key in selected)}

    return {"overall": summarize(keys), "per_locus": {
        locus: summarize({key for key in keys if key[1] == locus})
        for locus in sorted({key[1] for key in keys})
    }}


def audit_calls(assembly, reads):
    """Retain uncertainty and sample identity while auditing numeric agreement."""
    assembly_counts = {key: row["repeat_count"] for key, row in assembly.items()}
    read_counts = {key: row["repeat_count"] for key, row in reads.items()}
    keys = set(assembly) | set(reads)

    def subset(selected):
        return compare({key: assembly_counts.get(key) for key in selected},
                       {key: read_counts.get(key) for key in selected})

    statuses = sorted({row.get("status") or "UNKNOWN" for row in reads.values()})
    methods = sorted({row.get("estimation_method") or "unspecified" for row in reads.values()
                      if row.get("status") == "ESTIMATED"})
    assembly_samples = {sample for sample, _locus in assembly}
    read_samples = {sample for sample, _locus in reads}
    disagreements = []
    for key in sorted(keys):
        expected, observed = assembly_counts.get(key), read_counts.get(key)
        if expected is not None and observed is not None and expected == observed:
            continue
        row = reads.get(key, {})
        disagreements.append({
            "sample_id": key[0], "locus_id": key[1], "assembly_repeat": expected,
            "read_repeat": observed, "signed_error": observed - expected
            if expected is not None and observed is not None else None,
            "assembly_status": assembly.get(key, {}).get("status", "UNKNOWN"),
            "read_status": row.get("status", "UNKNOWN"),
            **{field: row.get(field, "") for field in (
                "inference_method", "estimation_method", "estimation_reason",
                "estimate_support", "read_depth", "primary_read_depth", "call_file",
            )},
            "comparison": "disagreement" if expected is not None and observed is not None else "missing_call",
        })
    disagreements.sort(key=lambda row: (row["signed_error"] is None,
                       -abs(row["signed_error"] or 0), row["sample_id"], row["locus_id"]))
    return {
        **compare(assembly_counts, read_counts),
        "sample_matching": {"shared": sorted(assembly_samples & read_samples),
                            "assembly_only": sorted(assembly_samples - read_samples),
                            "reads_only": sorted(read_samples - assembly_samples)},
        "by_read_status": {status: subset({key for key, row in reads.items()
                                           if (row.get("status") or "UNKNOWN") == status})
                           for status in statuses},
        "by_estimation_method": {method: subset({key for key, row in reads.items()
                                  if row.get("status") == "ESTIMATED"
                                  and (row.get("estimation_method") or "unspecified") == method})
                                 for method in methods},
        "per_sample": {sample: subset({key for key in keys if key[0] == sample})["overall"]
                       for sample in sorted(assembly_samples | read_samples)},
        "disagreements": disagreements,
        "interpretation": "Paired locus-call agreement, not a Mantel test. Status strata retain missingness; estimates remain in the all-call result. Positive signed error means FASTQ overcalls repeats. Assembly is a comparator, not guaranteed truth.",
    }


def worker(args):
    import resource
    import time

    sys.path.insert(0, str(Path(args.source).resolve()))
    from mlvamaps.short_reads import run_short_read_call

    start = time.perf_counter()
    run_short_read_call(
        args.reads1, args.reads2, args.panel, args.outdir, args.sample_id,
        database_path=args.database, threads=args.threads, show_progress=False,
    )
    elapsed = time.perf_counter() - start
    rss = max(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    diagnostics = Path(args.outdir) / "short_read_diagnostics.tsv"
    fraction = None
    if diagnostics.exists():
        with diagnostics.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        fraction = sum(row["microassembly_attempted"] == "yes" for row in rows) / len(rows) if rows else 0
    print(json.dumps({"seconds": elapsed, "peak_rss_mib": rss / (1024 ** 2 if sys.platform == "darwin" else 1024),
                      "microassembly_fraction": fraction}))


def benchmark(args):
    output = Path(args.outdir).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Benchmark output must be empty to avoid mixing trials")
    output.mkdir(parents=True, exist_ok=True)
    results = {"baseline": [], "updated": []}
    for trial in range(args.trials):
        for version, source in (("baseline", args.baseline_source), ("updated", args.source)):
            directory = output / f"{version}_{trial}"
            command = [sys.executable, str(Path(__file__).resolve()), "worker",
                       "--source", source, "--reads1", args.reads1, "--reads2", args.reads2,
                       "--outdir", str(directory), "--threads", str(args.threads), "--sample-id", args.sample_id]
            for name in ("panel", "database"):
                if getattr(args, name):
                    command.extend(["--" + name, getattr(args, name)])
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            measurement = json.loads(result.stdout.splitlines()[-1])
            if args.assembly:
                measurement["accuracy"] = compare(read_calls(args.assembly), read_calls(directory))
            results[version].append(measurement)
    medians = {version: statistics.median(trial["seconds"] for trial in trials) for version, trials in results.items()}
    return {"trials": results, "median_seconds": medians,
            "runtime_percent_difference": (medians["updated"] / medians["baseline"] - 1) * 100,
            "method": "Alternating serial fresh processes; wall time excludes imports. RSS is max(parent, largest child), not concurrent total. OS caches are not cleared."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    comparison = commands.add_parser("compare")
    comparison.add_argument("--assembly", required=True)
    comparison.add_argument("--reads", "--short-reads", dest="short_reads", required=True,
                            help="Short- or long-read calls.tsv or result directory")
    comparison.add_argument("--baseline")
    for name in ("benchmark", "worker"):
        command = commands.add_parser(name)
        command.add_argument("--source", default=str(Path(__file__).resolve().parents[1]))
        command.add_argument("--reads1", required=True)
        command.add_argument("--reads2", required=True)
        command.add_argument("--panel", required=True)
        command.add_argument("--database")
        command.add_argument("--outdir", required=True)
        command.add_argument("--sample-id", default="sample")
        command.add_argument("--threads", type=int, default=2)
        if name == "benchmark":
            command.add_argument("--baseline-source", required=True)
            command.add_argument("--assembly")
            command.add_argument("--trials", type=int, default=3)
    args = parser.parse_args()
    if args.command == "worker":
        worker(args)
        return
    if args.command == "benchmark":
        if args.trials < 1:
            parser.error("--trials must be positive")
        result = benchmark(args)
    else:
        assembly = read_calls(args.assembly, with_metadata=True)
        result = {"updated": audit_calls(assembly, read_calls(args.short_reads, with_metadata=True))}
        if args.baseline:
            result["baseline"] = audit_calls(assembly, read_calls(args.baseline, with_metadata=True))
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
