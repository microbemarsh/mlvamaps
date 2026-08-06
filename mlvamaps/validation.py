from __future__ import annotations

import csv
from pathlib import Path

from .io import open_text, write_tsv


VALIDATION_FIELDS = [
    "sample_id",
    "locus_id",
    "workflow",
    "truth_repeat_count",
    "observed_repeat_count",
    "repeat_count_min",
    "repeat_count_max",
    "classification",
]
VALIDATION_SUMMARY_FIELDS = [
    "workflow",
    "total_truth_loci",
    "exact_calls",
    "exact_matches",
    "interval_calls",
    "interval_matches",
    "incorrect_calls",
    "unresolved_loci",
    "false_positives",
    "false_negatives",
    "exact_call_accuracy",
    "callable_locus_fraction",
    "interval_coverage",
    "false_exact_call_rate",
    "profile_recovery_rate",
]


def _read_calls(path: str | Path) -> list[dict[str, str]]:
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _value(row: dict[str, str], *names: str) -> str:
    for name in names:
        if row.get(name, "") not in ("", None):
            return str(row[name])
    return ""


def _number(value: str) -> float | None:
    try:
        return float(value) if value != "" else None
    except ValueError:
        return None


def compare_call_sets(
    truth_rows: list[dict[str, str]],
    observed_rows: list[dict[str, str]],
    workflow: str,
) -> list[dict]:
    truth = {(row.get("sample_id", ""), row["locus_id"]): row for row in truth_rows}
    observed = {(row.get("sample_id", ""), row["locus_id"]): row for row in observed_rows}
    output = []
    for key in sorted(set(truth) | set(observed)):
        truth_row = truth.get(key, {})
        observed_row = observed.get(key, {})
        truth_value = _value(truth_row, "repeat_count", "called_repeat_count")
        observed_value = _value(observed_row, "repeat_count", "called_repeat_count")
        lower = _value(observed_row, "repeat_count_min")
        upper = _value(observed_row, "repeat_count_max")
        truth_number = _number(truth_value)
        observed_number = _number(observed_value)
        lower_number = _number(lower)
        upper_number = _number(upper)
        if truth_number is None and observed_number is not None:
            classification = "false_positive"
        elif truth_number is not None and observed_number is not None:
            classification = "exact_match" if truth_number == observed_number else "incorrect_call"
        elif truth_number is not None and lower_number is not None and upper_number is not None:
            classification = "within_reported_interval" if lower_number <= truth_number <= upper_number else "incorrect_call"
        elif truth_number is not None and observed_row:
            classification = "unresolved"
        elif truth_number is not None:
            classification = "false_negative"
        else:
            classification = "unresolved"
        output.append(
            {
                "sample_id": key[0],
                "locus_id": key[1],
                "workflow": workflow,
                "truth_repeat_count": truth_value,
                "observed_repeat_count": observed_value,
                "repeat_count_min": lower,
                "repeat_count_max": upper,
                "classification": classification,
            }
        )
    return output


def summarize_validation(rows: list[dict]) -> list[dict]:
    output = []
    for workflow in sorted({str(row["workflow"]) for row in rows}):
        subset = [row for row in rows if row["workflow"] == workflow]
        counts = {classification: sum(row["classification"] == classification for row in subset) for classification in {str(row["classification"]) for row in subset}}
        truth_loci = sum(row["truth_repeat_count"] != "" for row in subset)
        exact_calls = sum(row["observed_repeat_count"] != "" for row in subset)
        exact_matches = counts.get("exact_match", 0)
        interval_calls = sum(row["observed_repeat_count"] == "" and row["repeat_count_min"] != "" for row in subset)
        interval_matches = counts.get("within_reported_interval", 0)
        incorrect = counts.get("incorrect_call", 0)
        unresolved = counts.get("unresolved", 0)
        false_positive = counts.get("false_positive", 0)
        false_negative = counts.get("false_negative", 0)
        by_sample: dict[str, list[dict]] = {}
        for row in subset:
            by_sample.setdefault(str(row["sample_id"]), []).append(row)
        recovered_profiles = sum(
            bool(sample_rows)
            and all(
                row["truth_repeat_count"] == ""
                or row["classification"] == "exact_match"
                for row in sample_rows
            )
            for sample_rows in by_sample.values()
        )
        output.append(
            {
                "workflow": workflow,
                "total_truth_loci": truth_loci,
                "exact_calls": exact_calls,
                "exact_matches": exact_matches,
                "interval_calls": interval_calls,
                "interval_matches": interval_matches,
                "incorrect_calls": incorrect,
                "unresolved_loci": unresolved,
                "false_positives": false_positive,
                "false_negatives": false_negative,
                "exact_call_accuracy": round(exact_matches / max(exact_calls, 1), 6),
                "callable_locus_fraction": round((exact_calls + interval_calls) / max(truth_loci, 1), 6),
                "interval_coverage": round(interval_matches / max(interval_calls, 1), 6),
                "false_exact_call_rate": round((incorrect + false_positive) / max(exact_calls, 1), 6),
                "profile_recovery_rate": round(recovered_profiles / max(len(by_sample), 1), 6),
            }
        )
    return output


def run_validation(
    truth_path: str,
    outdir: str,
    long_read_path: str | None = None,
    illumina_path: str | None = None,
) -> dict[str, Path]:
    if not long_read_path and not illumina_path:
        raise ValueError("validation requires --long-read and/or --illumina calls")
    truth = _read_calls(truth_path)
    rows: list[dict] = []
    if long_read_path:
        rows.extend(compare_call_sets(truth, _read_calls(long_read_path), "accurate-long"))
    if illumina_path:
        rows.extend(compare_call_sets(truth, _read_calls(illumina_path), "illumina"))
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    details = output / "validation_details.tsv"
    summary = output / "validation_summary.tsv"
    write_tsv(rows, details, VALIDATION_FIELDS)
    write_tsv(summarize_validation(rows), summary, VALIDATION_SUMMARY_FIELDS)
    return {"details": details, "summary": summary}
