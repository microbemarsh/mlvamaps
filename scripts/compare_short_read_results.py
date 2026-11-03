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
        diagnostics = {}
        if with_metadata:
            for name, filename, sample_field, locus_field in (
                ("genotype_evidence", "genotype_evidence.tsv", "sample", "locus"),
                ("short_read_diagnostics", "short_read_diagnostics.tsv", "sample_id", "locus_id"),
                ("reconstruction", "query_locus_reconstruction.tsv", "sample_id", "locus_id"),
                ("insert_calibration", "library_insert_distribution.tsv", "sample", None),
                ("run_summary", "sample_summary.tsv", "sample_id", None),
            ):
                source = file.parent / filename
                if not source.exists():
                    continue
                rows = {}
                with source.open(newline="") as handle:
                    for row in csv.DictReader(handle, delimiter="\t"):
                        key = (row[sample_field], row[locus_field] if locus_field else None)
                        if key in rows:
                            raise ValueError(f"Duplicate diagnostic sample/locus {key} in {source}")
                        rows[key] = {**row, "source_file": str(source.resolve())}
                diagnostics[name] = (rows, locus_field)
        with file.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = (row["sample_id"], row["locus_id"])
                if key in calls:
                    raise ValueError(f"Duplicate sample/locus {key} in {file}")
                value = row.get("repeat_count", "")
                count = float(value) if value not in ("", "NA", "N/A") else None
                if count is not None and not math.isfinite(count):
                    raise ValueError(f"Non-finite count for {key}")
                metadata = {name: rows[(key[0], key[1] if locus_field else None)]
                            for name, (rows, locus_field) in diagnostics.items()
                            if (key[0], key[1] if locus_field else None) in rows}
                calls[key] = {**row, "repeat_count": count, "call_file": str(file.resolve()),
                              "diagnostics": metadata} if with_metadata else count
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
        "greater_than_2": sum(error > 2 for error in errors) / len(errors) if errors else None,
        "max_absolute_error": max(errors) if errors else None,
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
        known = sum(assembly.get(key) is not None for key in selected)
        return {**metrics(pairs), "total_sample_loci": len(selected),
                "assembly_numeric_loci": known,
                "exact_agreement_all_assembly_numeric": sum(a == b for a, b in pairs) / known if known else None,
                "missing_assembly": sum(assembly.get(key) is None for key in selected),
                "missing_short_read": sum(short_reads.get(key) is None for key in selected)}

    return {"overall": summarize(keys), "per_locus": {
        locus: summarize({key for key in keys if key[1] == locus})
        for locus in sorted({key[1] for key in keys})
    }}


def evidence_observations(row):
    """Describe recorded evidence, without assigning a cause to a discrepancy."""
    diagnostics = row.get("diagnostics", {})
    evidence = diagnostics.get("genotype_evidence", {})
    fields = (
        ("genotype", evidence.get("genotype_status") or row.get("genotype_status")),
        ("estimate", row.get("estimation_method") or evidence.get("estimate_source")),
        ("microassembly", diagnostics.get("short_read_diagnostics", {}).get("microassembly_reason")),
        ("insert_calibration", diagnostics.get("insert_calibration", {}).get("source")),
        ("insert_calibration_status", diagnostics.get("insert_calibration", {}).get("calibration_status")),
        ("reconstruction", diagnostics.get("reconstruction", {}).get("status")),
        ("computational_fallback", evidence.get("computational_fallback") or row.get("computational_fallback")),
    )
    return [f"{name}:{value}" for name, value in fields if value]


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
    observations = {key: evidence_observations(row) for key, row in reads.items()}
    presence = {key: (assembly.get(key, {}).get("present"), reads.get(key, {}).get("present")) for key in keys}
    known_presence = {key: values for key, values in presence.items() if all(value in {"yes", "no"} for value in values)}
    disagreements = []
    for key in sorted(keys):
        expected, observed = assembly_counts.get(key), read_counts.get(key)
        presence_disagrees = key in known_presence and known_presence[key][0] != known_presence[key][1]
        if expected is not None and observed is not None and expected == observed and not presence_disagrees:
            continue
        row = reads.get(key, {})
        disagreements.append({
            "sample_id": key[0], "locus_id": key[1], "assembly_repeat": expected,
            "read_repeat": observed, "signed_error": observed - expected
            if expected is not None and observed is not None else None,
            "assembly_status": assembly.get(key, {}).get("status", "UNKNOWN"),
            "read_status": row.get("status", "UNKNOWN"),
            "assembly_present": presence[key][0], "read_present": presence[key][1],
            "presence_comparison": "disagreement" if presence_disagrees else "agreement" if key in known_presence else "unknown",
            **{field: row.get(field, "") for field in (
                "inference_method", "estimation_method", "estimation_reason",
                "estimate_support", "read_depth", "primary_read_depth", "call_file",
            )},
            "comparison": ("presence_disagreement" if expected == observed else "disagreement")
                          if expected is not None and observed is not None else "missing_call",
            "missing_side": "both" if expected is None and observed is None else "assembly" if expected is None else "reads" if observed is None else "",
            "evidence_observations": observations.get(key, []),
            "diagnostics": row.get("diagnostics", {}),
        })
    disagreements.sort(key=lambda row: (row["signed_error"] is None,
                       -abs(row["signed_error"] or 0), row["sample_id"], row["locus_id"]))
    return {
        **compare(assembly_counts, read_counts),
        "sample_matching": {"shared": sorted(assembly_samples & read_samples),
                            "assembly_only": sorted(assembly_samples - read_samples),
                            "reads_only": sorted(read_samples - assembly_samples)},
        "presence": {"compared": len(known_presence), "unknown": len(keys) - len(known_presence),
                     "agreements": sum(a == b for a, b in known_presence.values()),
                     "disagreements": sum(a != b for a, b in known_presence.values())},
        "by_read_status": {status: subset({key for key, row in reads.items()
                                           if (row.get("status") or "UNKNOWN") == status})
                           for status in statuses},
        "by_estimation_method": {method: subset({key for key, row in reads.items()
                                  if row.get("status") == "ESTIMATED"
                                  and (row.get("estimation_method") or "unspecified") == method})
                                 for method in methods},
        "by_evidence_observation": {label: subset({key for key, labels in observations.items() if label in labels})
                                    for label in sorted({label for labels in observations.values() for label in labels})},
        "per_sample": {sample: subset({key for key in keys if key[0] == sample})["overall"]
                       for sample in sorted(assembly_samples | read_samples)},
        "disagreements": disagreements,
        "interpretation": "Paired locus-call agreement, not a Mantel test. Exact agreement among numeric pairs is separate from agreement across all assembly-numeric loci, which counts missing FASTQ calls as failures. Estimates remain included. Presence compares recorded yes/no fields; NOT_FOUND does not establish biological absence. Evidence observations describe outputs, not proven causes. Positive signed error means FASTQ overcalls repeats. Assembly is a comparator, not guaranteed truth.",
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
        sample_mode=args.sample_mode, assembly_round_tolerance=args.assembly_round_tolerance,
        short_insert_mean=args.short_insert_mean, short_insert_sd=args.short_insert_sd,
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
    performance = Path(args.outdir) / "fastq_performance.tsv"
    stage_metrics = {}
    if performance.exists():
        with performance.open() as handle:
            stage_metrics = {row["metric"]: float(row["value"]) for row in csv.DictReader(handle, delimiter="\t")}
    print(json.dumps({"seconds": elapsed, "peak_rss_mib": rss / (1024 ** 2 if sys.platform == "darwin" else 1024),
                      "microassembly_fraction": fraction, "stage_metrics": stage_metrics}))


def benchmark(args):
    output = Path(args.outdir).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Benchmark output must be empty to avoid mixing trials")
    output.mkdir(parents=True, exist_ok=True)
    results = {"baseline": [], "updated": []}
    for trial in range(args.trials):
        sources = [("baseline", args.baseline_source), ("updated", args.source)]
        for version, source in sources if trial % 2 == 0 else reversed(sources):
            directory = output / f"{version}_{trial}"
            command = [sys.executable, str(Path(__file__).resolve()), "worker",
                       "--source", source, "--reads1", args.reads1, "--reads2", args.reads2,
                       "--outdir", str(directory), "--threads", str(args.threads), "--sample-id", args.sample_id]
            for name in ("panel", "database", "sample_mode", "assembly_round_tolerance", "short_insert_mean", "short_insert_sd"):
                if getattr(args, name) is not None:
                    command.extend(["--" + name.replace("_", "-"), str(getattr(args, name))])
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            measurement = json.loads(result.stdout.splitlines()[-1])
            if args.assembly:
                measurement["accuracy"] = audit_calls(read_calls(args.assembly, with_metadata=True),
                                                      read_calls(directory, with_metadata=True))
            results[version].append(measurement)
    medians = {version: statistics.median(trial["seconds"] for trial in trials) for version, trials in results.items()}
    return {"trials": results, "median_seconds": medians,
            "parameters": {name: getattr(args, name) for name in (
                "sample_mode", "assembly_round_tolerance", "short_insert_mean", "short_insert_sd", "threads")},
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
        command.add_argument("--sample-mode", choices=("isolate", "metagenome"), default="isolate")
        command.add_argument("--assembly-round-tolerance", type=float, default=.25)
        command.add_argument("--short-insert-mean", type=float)
        command.add_argument("--short-insert-sd", type=float)
        if name == "benchmark":
            command.add_argument("--baseline-source", required=True)
            command.add_argument("--assembly")
            command.add_argument("--trials", type=int, default=3)
    args = parser.parse_args()
    if args.command in {"worker", "benchmark"}:
        if not 0 <= args.assembly_round_tolerance <= .5:
            parser.error("--assembly-round-tolerance must be between 0 and 0.5")
        inserts = (args.short_insert_mean, args.short_insert_sd)
        if (inserts[0] is None) != (inserts[1] is None):
            parser.error("--short-insert-mean and --short-insert-sd must be supplied together")
        if any(value is not None and (not math.isfinite(value) or value <= 0) for value in inserts):
            parser.error("Insert mean and standard deviation must be finite and positive")
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
