from __future__ import annotations

import csv
import math
import re
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

from .combined_marker_export import COMBINED_OUTPUT_NAMES, export_combined_markers
from .io import open_text, write_tsv
from .phylogeny import neighbor_joining_tree_from_matrix
from .primers import read_loci_or_primers
from .sample_metadata import METADATA_ALIASES


PAIRWISE_FIELDS = [
    "sample_1",
    "sample_2",
    "loci_compared",
    "fraction_loci_compared",
    "categorical_differences",
    "categorical_distance",
    "repeat_distance_raw",
    "repeat_distance",
    "comparison_status",
]

EXCLUDED_FIELDS = [
    "sample_id",
    "reason",
    "scope",
    "path",
    "callable_loci",
    "total_loci",
    "callable_fraction",
    "details",
]

USED_FIELDS = [
    "sample_id",
    "path",
    "callable_loci",
    "total_loci",
    "callable_fraction",
    "metadata_found",
    "coordinates_valid",
    "coordinate_status",
]

OUTPUT_NAMES = (
    "myoga_metadata.tsv",
    "mlva_profiles.tsv",
    "mlva_calls_long.tsv",
    "mlva_pairwise_distances.tsv",
    "mlva_distance_matrix.tsv",
    "mlva_nj.tree",
    "samples_used.tsv",
    "samples_excluded.tsv",
    "export_summary.tsv",
    "export_summary.txt",
    *COMBINED_OUTPUT_NAMES,
)

_MISSING = {"", ".", "na", "nan", "none", "null"}
_NATURAL_PART = re.compile(r"(\d+)")


@dataclass
class SampleCalls:
    sample_id: str
    path: Path
    rows: list[dict[str, str]]
    locus_order: list[str]
    aggregate_source: bool = False


@dataclass
class MetadataTable:
    fields: list[str]
    rows_by_id: dict[str, dict[str, str]]
    id_field: str
    latitude_field: str | None
    longitude_field: str | None


def _natural_key(value: str) -> tuple:
    return tuple(
        int(part) if part.isdigit() else part.casefold()
        for part in _NATURAL_PART.split(str(value))
    )


def _read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with open_text(path, "rt") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            fields = [str(field).strip() for field in (reader.fieldnames or [])]
            if not fields:
                raise ValueError("table has no header")
            rows = []
            for number, raw in enumerate(reader, start=2):
                if None in raw:
                    raise ValueError(f"row {number} has more values than the header")
                rows.append(
                    {
                        str(key).strip(): "" if value is None else str(value).strip()
                        for key, value in raw.items()
                    }
                )
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValueError(f"could not read {path}: {exc}") from exc
    return fields, rows


def _numeric_repeat(value: object) -> float | None:
    text = "" if value is None else str(value).strip()
    if text.casefold() in _MISSING:
        return None
    try:
        parsed = float(text)
    except ValueError as exc:
        raise ValueError(f"non-numeric repeat_count {text!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"non-finite repeat_count {text!r}")
    return parsed


def _format_number(value: float | None, digits: int = 8) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.{digits}f}"


def _format_allele(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return ""
    return f"{value:.15g}"


def _status_records(results: Path) -> tuple[dict[str, list[tuple[str, Path, str]]], list[dict]]:
    statuses: dict[str, list[tuple[str, Path, str]]] = defaultdict(list)
    malformed: list[dict] = []
    for path in sorted(results.rglob("batch_status.tsv")):
        try:
            fields, rows = _read_table(path)
            if not {"sample_id", "status"}.issubset(fields):
                raise ValueError("requires sample_id and status columns")
            for row in rows:
                sample_id = row.get("sample_id", "").strip()
                if not sample_id:
                    continue
                statuses[sample_id].append(
                    (row.get("status", "").casefold(), path.parent, row.get("message", ""))
                )
        except ValueError as exc:
            malformed.append(
                _excluded(path.parent.name, "MALFORMED_RESULTS", "tree", path, details=str(exc))
            )
    return statuses, malformed


def _excluded(
    sample_id: str,
    reason: str,
    scope: str,
    path: str | Path,
    callable_loci: int | str = "",
    total_loci: int | str = "",
    callable_fraction: float | str = "",
    details: str = "",
) -> dict:
    # csv permits embedded newlines, but these audit files are routinely read by
    # line-oriented TSV tools (awk, cut, R read.delim).  Batch failures can carry
    # multi-line exception messages, so keep every exclusion to exactly one row.
    details = re.sub(r"\s+", " ", str(details)).strip()
    fraction = (
        callable_fraction
        if callable_fraction == ""
        else _format_number(float(callable_fraction), 6)
    )
    return {
        "sample_id": sample_id,
        "reason": reason,
        "scope": scope,
        "path": str(path),
        "callable_loci": callable_loci,
        "total_loci": total_loci,
        "callable_fraction": fraction,
        "details": details,
    }


def discover_sample_calls(results_path: str | Path) -> tuple[list[SampleCalls], list[dict], int]:
    """Discover valid per-sample calls, preferring leaf files over batch aggregates."""
    results = Path(results_path)
    if not results.exists() or not results.is_dir():
        raise ValueError(f"Results directory does not exist or is not a directory: {results}")

    statuses, excluded = _status_records(results)
    candidates: dict[str, list[SampleCalls]] = defaultdict(list)
    calls_paths = sorted(results.rglob("calls.tsv"))
    discovered_ids = set(statuses)
    for path in calls_paths:
        try:
            fields, raw_rows = _read_table(path)
            required = {"sample_id", "locus_id", "repeat_count"}
            if not required.issubset(fields):
                raise ValueError(
                    "missing required calls columns: " + ", ".join(sorted(required - set(fields)))
                )
            if not raw_rows:
                raise ValueError("calls table contains no locus rows")
            grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
            for row in raw_rows:
                sample_id = row.get("sample_id", "").strip()
                if not sample_id:
                    raise ValueError("calls row has an empty sample_id")
                grouped[sample_id].append(row)
            aggregate = len(grouped) > 1 or (path.parent / "batch_status.tsv").exists()
            for sample_id, rows in grouped.items():
                discovered_ids.add(sample_id)
                loci = [row.get("locus_id", "").strip() for row in rows]
                if any(not locus for locus in loci):
                    raise ValueError(f"sample {sample_id!r} has an empty locus_id")
                duplicates = sorted(
                    {locus for locus in loci if loci.count(locus) > 1}, key=_natural_key
                )
                if duplicates:
                    raise ValueError(
                        f"sample {sample_id!r} has duplicate loci: {', '.join(duplicates)}"
                    )
                for row in rows:
                    _numeric_repeat(row.get("repeat_count"))
                candidates[sample_id].append(
                    SampleCalls(sample_id, path, rows, loci, aggregate_source=aggregate)
                )
        except ValueError as exc:
            discovered_ids.add(path.parent.name)
            excluded.append(
                _excluded(path.parent.name, "MALFORMED_RESULTS", "tree", path, details=str(exc))
            )

    selected: list[SampleCalls] = []
    for sample_id in sorted(candidates, key=_natural_key):
        records = candidates[sample_id]
        leaf = [record for record in records if not record.aggregate_source]
        usable = leaf or records
        if len(usable) > 1:
            excluded.append(
                _excluded(
                    sample_id,
                    "DUPLICATE_SAMPLE_ID",
                    "tree",
                    ";".join(str(record.path) for record in usable),
                    details="multiple independent calls tables record this sample_id",
                )
            )
            continue
        record = usable[0]
        sample_statuses = statuses.get(sample_id, [])
        failed = [entry for entry in sample_statuses if entry[0] == "failed"]
        if failed:
            excluded.append(
                _excluded(
                    sample_id,
                    "FAILED_BATCH_SAMPLE",
                    "tree",
                    record.path.parent,
                    details=failed[-1][2],
                )
            )
            continue
        summary_path = record.path.parent / "sample_summary.tsv"
        if summary_path.exists():
            try:
                summary_fields, summary_rows = _read_table(summary_path)
                if "run_status" in summary_fields and summary_rows:
                    run_status = summary_rows[0].get("run_status", "").casefold()
                    if run_status and run_status != "success":
                        excluded.append(
                            _excluded(
                                sample_id,
                                "FAILED_BATCH_SAMPLE",
                                "tree",
                                record.path.parent,
                                details=f"sample_summary run_status={run_status}",
                            )
                        )
                        continue
            except ValueError as exc:
                excluded.append(
                    _excluded(sample_id, "MALFORMED_RESULTS", "tree", summary_path, details=str(exc))
                )
                continue
        selected.append(record)

    present_ids = set(candidates)
    for sample_id in sorted(statuses, key=_natural_key):
        if sample_id in present_ids:
            continue
        entries = statuses[sample_id]
        failed = [entry for entry in entries if entry[0] == "failed"]
        status, root, message = failed[-1] if failed else entries[-1]
        reason = "FAILED_BATCH_SAMPLE" if status == "failed" else "MISSING_CALLS_FILE"
        excluded.append(
            _excluded(sample_id, reason, "tree", root / sample_id, details=message)
        )

    marker_names = (
        "locus_repeat_counts.tsv",
        "mlva_fingerprint.tsv",
        "sample_summary.tsv",
        "profile_matches.tsv",
    )
    marker_directories = {
        path.parent
        for name in marker_names
        for path in results.rglob(name)
        if not (path.parent / "calls.tsv").exists()
        and not (path.parent / "batch_status.tsv").exists()
    }
    for directory in sorted(marker_directories):
        sample_id = directory.name
        marker = next((directory / name for name in marker_names if (directory / name).exists()), None)
        if marker is not None:
            try:
                fields, rows = _read_table(marker)
                recorded = {
                    row.get("sample_id", "").strip()
                    for row in rows
                    if row.get("sample_id", "").strip()
                }
                if len(recorded) == 1:
                    sample_id = recorded.pop()
            except ValueError:
                pass
        if sample_id not in present_ids and sample_id not in statuses:
            excluded.append(
                _excluded(sample_id, "MISSING_CALLS_FILE", "tree", directory)
            )
        discovered_ids.add(sample_id)

    return selected, excluded, len(discovered_ids)


def _resolve_field(fields: list[str], requested: str, aliases: Iterable[str] = ()) -> str | None:
    by_lower = {field.casefold(): field for field in fields}
    if requested.casefold() in by_lower:
        return by_lower[requested.casefold()]
    for alias in aliases:
        if alias.casefold() in by_lower:
            return by_lower[alias.casefold()]
    return None


def read_export_metadata(
    path: str | Path,
    id_column: str = "shared_identifier",
    latitude_column: str = "latitude",
    longitude_column: str = "longitude",
) -> MetadataTable:
    metadata_path = Path(path)
    try:
        with open_text(metadata_path, "rt") as handle:
            sample = handle.read(4096)
            handle.seek(0)
            first = sample.splitlines()[0] if sample.splitlines() else ""
            if "\t" in first:
                delimiter = "\t"
            elif "," in first:
                delimiter = ","
            else:
                raise ValueError("metadata must be comma- or tab-delimited")
            reader = csv.DictReader(handle, delimiter=delimiter)
            fields = [str(field).strip() for field in (reader.fieldnames or [])]
            if not fields:
                raise ValueError("metadata has no header")
            rows = []
            for number, raw in enumerate(reader, start=2):
                if None in raw:
                    raise ValueError(f"metadata row {number} has more values than the header")
                rows.append(
                    {
                        str(key).strip(): "" if value is None else str(value).strip()
                        for key, value in raw.items()
                    }
                )
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValueError(f"Could not read metadata table {metadata_path}: {exc}") from exc

    id_field = _resolve_field(fields, id_column)
    if id_field is None:
        raise ValueError(f"Metadata identifier column {id_column!r} was not found")
    latitude_field = _resolve_field(
        fields,
        latitude_column,
        METADATA_ALIASES["latitude"] if latitude_column == "latitude" else (),
    )
    longitude_field = _resolve_field(
        fields,
        longitude_column,
        METADATA_ALIASES["longitude"] if longitude_column == "longitude" else (),
    )
    rows_by_id: dict[str, dict[str, str]] = {}
    for number, row in enumerate(rows, start=2):
        identifier = row.get(id_field, "").strip()
        if not identifier:
            raise ValueError(f"Metadata row {number} has an empty {id_field!r} value")
        if identifier in rows_by_id:
            raise ValueError(f"Duplicate metadata identifier {identifier!r}")
        rows_by_id[identifier] = row
    return MetadataTable(fields, rows_by_id, id_field, latitude_field, longitude_field)


def _coordinate_status(
    row: dict[str, str] | None,
    metadata: MetadataTable,
) -> tuple[str, str, str]:
    if row is None:
        return "METADATA_NOT_FOUND", "", ""
    latitude = row.get(metadata.latitude_field, "") if metadata.latitude_field else ""
    longitude = row.get(metadata.longitude_field, "") if metadata.longitude_field else ""
    if not latitude.strip() or not longitude.strip():
        return "MISSING_COORDINATES", "", ""
    try:
        lat = float(latitude)
        lon = float(longitude)
    except ValueError:
        return "INVALID_COORDINATES", "", ""
    if not math.isfinite(lat) or not math.isfinite(lon) or not (-90 <= lat <= 90) or not (-180 <= lon <= 180):
        return "INVALID_COORDINATES", "", ""
    return "VALID", _format_allele(lat), _format_allele(lon)


def build_call_matrix(
    samples: list[SampleCalls], locus_order: list[str]
) -> tuple[np.ndarray, dict[str, dict[str, dict[str, str]]]]:
    matrix = np.full((len(samples), len(locus_order)), np.nan, dtype=np.float64)
    rows_by_sample: dict[str, dict[str, dict[str, str]]] = {}
    locus_index = {locus: index for index, locus in enumerate(locus_order)}
    for sample_index, sample in enumerate(samples):
        rows_by_locus = {row["locus_id"]: row for row in sample.rows}
        rows_by_sample[sample.sample_id] = rows_by_locus
        for locus, row in rows_by_locus.items():
            value = _numeric_repeat(row.get("repeat_count"))
            if value is not None:
                matrix[sample_index, locus_index[locus]] = value
    return matrix, rows_by_sample


def calculate_pairwise_distances(
    sample_ids: list[str],
    calls: np.ndarray,
    min_pairwise_loci: int = 1,
    min_pairwise_fraction: float = 0.5,
    *,
    applicable: np.ndarray | None = None,
    row_sink: Callable[[dict], object] | None = None,
    retain: str = "both",
) -> tuple[list[dict], np.ndarray, np.ndarray, np.ndarray]:
    """Calculate shared-call categorical and repeat distances in bounded-memory blocks."""
    if retain not in {"both", "categorical", "repeat"}:
        raise ValueError("retain must be 'both', 'categorical', or 'repeat'")
    sample_count, locus_count = calls.shape
    if applicable is None:
        applicable = np.ones(calls.shape, dtype=bool)
    else:
        applicable = np.asarray(applicable, dtype=bool)
        if applicable.shape != calls.shape:
            raise ValueError("applicable must have the same shape as calls")
    categorical = (
        np.full((sample_count, sample_count), np.nan, dtype=np.float64)
        if retain in {"both", "categorical"}
        else np.empty((0, 0), dtype=np.float64)
    )
    repeat = (
        np.full((sample_count, sample_count), np.nan, dtype=np.float64)
        if retain in {"both", "repeat"}
        else np.empty((0, 0), dtype=np.float64)
    )
    overlap = np.zeros((sample_count, sample_count), dtype=np.int32)
    if categorical.size:
        np.fill_diagonal(categorical, 0.0)
    if repeat.size:
        np.fill_diagonal(repeat, 0.0)
    np.fill_diagonal(overlap, applicable.sum(axis=1))
    rows: list[dict] = []
    for left in range(sample_count):
        if left + 1 >= sample_count:
            continue
        right_calls = calls[left + 1 :]
        valid = np.isfinite(right_calls) & np.isfinite(calls[left])[None, :]
        compared = valid.sum(axis=1)
        shared_applicable = (
            applicable[left + 1 :] & applicable[left][None, :]
        ).sum(axis=1)
        differences = np.abs(right_calls - calls[left][None, :])
        categorical_raw = ((differences > 1e-12) & valid).sum(axis=1)
        repeat_raw = np.where(valid, differences, 0.0).sum(axis=1)
        fractions = np.divide(
            compared,
            shared_applicable,
            out=np.zeros_like(compared, dtype=float),
            where=shared_applicable > 0,
        )
        supported = (compared >= min_pairwise_loci) & (fractions >= min_pairwise_fraction)
        categorical_values = np.divide(
            categorical_raw,
            compared,
            out=np.full(len(compared), np.nan, dtype=float),
            where=supported,
        )
        repeat_values = np.divide(
            repeat_raw,
            compared,
            out=np.full(len(compared), np.nan, dtype=float),
            where=supported,
        )
        for offset, right in enumerate(range(left + 1, sample_count)):
            overlap[left, right] = overlap[right, left] = int(compared[offset])
            if supported[offset]:
                if categorical.size:
                    categorical[left, right] = categorical[right, left] = categorical_values[offset]
                if repeat.size:
                    repeat[left, right] = repeat[right, left] = repeat_values[offset]
            row = {
                "sample_1": sample_ids[left],
                "sample_2": sample_ids[right],
                "loci_compared": int(compared[offset]),
                "fraction_loci_compared": _format_number(float(fractions[offset])),
                "categorical_differences": int(categorical_raw[offset]),
                "categorical_distance": _format_number(float(categorical_values[offset])),
                "repeat_distance_raw": _format_number(float(repeat_raw[offset])),
                "repeat_distance": _format_number(float(repeat_values[offset])),
                "comparison_status": "sufficient" if supported[offset] else "insufficient_overlap",
            }
            if row_sink is None:
                rows.append(row)
            else:
                row_sink(row)
    return rows, categorical, repeat, overlap


def _complete_matrix_subset(
    distances: np.ndarray,
    callable_counts: np.ndarray,
    sample_ids: list[str],
) -> tuple[list[int], list[int]]:
    active = list(range(len(sample_ids)))
    removed: list[int] = []
    while len(active) > 1:
        submatrix = distances[np.ix_(active, active)]
        unsupported = ~np.isfinite(submatrix)
        np.fill_diagonal(unsupported, False)
        counts = unsupported.sum(axis=1)
        if int(counts.max(initial=0)) == 0:
            break
        worst = int(counts.max())
        choices = [index for index, count in enumerate(counts) if int(count) == worst]
        choice = min(
            choices,
            key=lambda position: (
                int(callable_counts[active[position]]),
                _natural_key(sample_ids[active[position]]),
            ),
        )
        removed.append(active.pop(choice))
    return active, removed


def _ordered_loci(samples: list[SampleCalls]) -> list[str]:
    if not samples:
        return []
    anchor = min(
        samples,
        key=lambda sample: (-len(sample.locus_order), _natural_key(sample.sample_id)),
    )
    ordered = list(anchor.locus_order)
    seen = set(ordered)
    extras = sorted(
        {locus for sample in samples for locus in sample.locus_order if locus not in seen},
        key=_natural_key,
    )
    return ordered + extras


def _write_metadata(
    sample_ids: list[str],
    metadata: MetadataTable,
    output: Path,
) -> tuple[dict[str, tuple[bool, str]], list[dict]]:
    renamed_fields: list[tuple[str, str]] = []
    reserved_counts: dict[str, int] = defaultdict(int)
    for field in metadata.fields:
        lower = field.casefold()
        target = field
        if lower in {"sample_id", "latitude", "longitude"}:
            target = f"metadata_{lower}"
        while target in {existing for _, existing in renamed_fields} or target in {
            "sample_id", "latitude", "longitude"
        }:
            reserved_counts[target] += 1
            target = f"{target}_{reserved_counts[target] + 1}"
        renamed_fields.append((field, target))
    fields = ["sample_id", "latitude", "longitude", *[target for _, target in renamed_fields]]
    statuses: dict[str, tuple[bool, str]] = {}
    rows = []
    issues = []
    for sample_id in sample_ids:
        source = metadata.rows_by_id.get(sample_id)
        coordinate_status, latitude, longitude = _coordinate_status(source, metadata)
        statuses[sample_id] = (source is not None, coordinate_status)
        row = {field: "" for field in fields}
        row.update({"sample_id": sample_id, "latitude": latitude, "longitude": longitude})
        if source is not None:
            for original, target in renamed_fields:
                row[target] = source.get(original, "")
        rows.append(row)
        if coordinate_status != "VALID":
            issues.append(
                _excluded(
                    sample_id,
                    coordinate_status,
                    "geography",
                    output,
                    details="sample remains in the MLVA tree; canonical coordinates are blank",
                )
            )
    write_tsv(rows, output, fields)
    return statuses, issues


def _write_text_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def export_myoga(
    results_path: str | Path,
    metadata_path: str | Path,
    outdir: str | Path,
    *,
    metadata_id: str = "shared_identifier",
    latitude: str = "latitude",
    longitude: str = "longitude",
    min_callable_fraction: float = 0.0,
    min_callable_loci: int = 0,
    min_pairwise_loci: int = 1,
    min_pairwise_fraction: float = 0.5,
    distance: str = "repeat",
    combined_markers: bool = False,
    loci_path: str | Path | None = None,
    snp_weight: float = 1.0,
    repeat_weight: float = 1.0,
    threads: int = 1,
    mafft_bin: str = "mafft",
    raxml_ng_bin: str = "raxml-ng",
    raxml_model: str = "DNA",
    force: bool = False,
) -> dict[str, Path | int | str]:
    """Export completed MLVAmaps results as a relatedness dataset for MYOGA."""
    if not 0 <= min_callable_fraction <= 1:
        raise ValueError("min_callable_fraction must be between 0 and 1")
    if min_callable_loci < 0:
        raise ValueError("min_callable_loci must be non-negative")
    if min_pairwise_loci < 1:
        raise ValueError("min_pairwise_loci must be at least 1")
    if not 0 <= min_pairwise_fraction <= 1:
        raise ValueError("min_pairwise_fraction must be between 0 and 1")
    if distance not in {"repeat", "categorical"}:
        raise ValueError("distance must be 'repeat' or 'categorical'")
    if snp_weight < 0 or repeat_weight < 0 or snp_weight + repeat_weight <= 0:
        raise ValueError("SNP and repeat weights must be non-negative with a positive total")

    output = Path(outdir)
    existing = [output / name for name in OUTPUT_NAMES if (output / name).exists()]
    if combined_markers and (output / "locus_trees").exists():
        existing.append(output / "locus_trees")
    if existing and not force:
        raise ValueError(
            f"Output files already exist in {output}; use --force to replace this export"
        )
    output.mkdir(parents=True, exist_ok=True)

    samples, exclusions, discovered = discover_sample_calls(results_path)
    metadata = read_export_metadata(metadata_path, metadata_id, latitude, longitude)
    locus_order = _ordered_loci(samples)
    total_loci = len(locus_order)
    matrix, rows_by_sample = build_call_matrix(samples, locus_order)
    locus_index = {locus: index for index, locus in enumerate(locus_order)}
    applicable_matrix = np.zeros((len(samples), total_loci), dtype=bool)
    for sample_index, sample in enumerate(samples):
        for locus in sample.locus_order:
            applicable_matrix[sample_index, locus_index[locus]] = True
    assayed_counts = applicable_matrix.sum(axis=1)
    callable_counts = np.isfinite(matrix).sum(axis=1) if total_loci else np.zeros(len(samples), dtype=int)
    callable_fractions = np.divide(
        callable_counts,
        assayed_counts,
        out=np.zeros(len(samples), dtype=float),
        where=assayed_counts > 0,
    )
    required_counts = np.maximum(
        min_callable_loci,
        np.ceil(min_callable_fraction * assayed_counts).astype(int),
    )

    threshold_indices: list[int] = []
    for index, sample in enumerate(samples):
        callable_loci = int(callable_counts[index])
        assayed_loci = int(assayed_counts[index])
        required_loci = int(required_counts[index])
        fraction = float(callable_fractions[index])
        if callable_loci == 0:
            reason = "NO_CALLABLE_LOCI"
        elif callable_loci < required_loci:
            reason = "TOO_FEW_CALLABLE_LOCI"
        else:
            threshold_indices.append(index)
            continue
        exclusions.append(
            _excluded(
                sample.sample_id,
                reason,
                "tree",
                sample.path.parent,
                callable_loci,
                assayed_loci,
                fraction,
                details=f"requires at least {required_loci} callable loci",
            )
        )

    threshold_samples = [samples[index] for index in threshold_indices]
    threshold_matrix = matrix[threshold_indices, :] if threshold_indices else np.empty((0, total_loci))
    threshold_applicable = (
        applicable_matrix[threshold_indices, :]
        if threshold_indices
        else np.empty((0, total_loci), dtype=bool)
    )
    threshold_counts = callable_counts[threshold_indices] if threshold_indices else np.empty(0, dtype=int)
    threshold_ids = [sample.sample_id for sample in threshold_samples]
    pairwise_path = output / "mlva_pairwise_distances.tsv"
    pairwise_temporary = pairwise_path.with_name(pairwise_path.name + ".tmp")
    with pairwise_temporary.open("w", newline="", encoding="utf-8") as handle:
        pairwise_writer = csv.DictWriter(
            handle, fieldnames=PAIRWISE_FIELDS, delimiter="\t", extrasaction="ignore"
        )
        pairwise_writer.writeheader()
        pairwise_rows, categorical_matrix, repeat_matrix, overlap_matrix = calculate_pairwise_distances(
            threshold_ids,
            threshold_matrix,
            min_pairwise_loci=min_pairwise_loci,
            min_pairwise_fraction=min_pairwise_fraction,
            applicable=threshold_applicable,
            row_sink=pairwise_writer.writerow,
            retain=distance,
        )
    pairwise_temporary.replace(pairwise_path)
    selected_matrix = repeat_matrix if distance == "repeat" else categorical_matrix
    retained_local, removed_local = _complete_matrix_subset(
        selected_matrix, threshold_counts, threshold_ids
    )
    for local_index in removed_local:
        sample = threshold_samples[local_index]
        exclusions.append(
            _excluded(
                sample.sample_id,
                "INSUFFICIENT_PAIRWISE_OVERLAP",
                "tree",
                sample.path.parent,
                int(threshold_counts[local_index]),
                int(assayed_counts[threshold_indices[local_index]]),
                float(callable_fractions[threshold_indices[local_index]]),
                details="removed deterministically so every retained pair meets overlap thresholds",
            )
        )

    final_samples = [threshold_samples[index] for index in retained_local]
    final_ids = [sample.sample_id for sample in final_samples]
    final_matrix = selected_matrix[np.ix_(retained_local, retained_local)]
    final_source_indexes = [threshold_indices[index] for index in retained_local]

    metadata_status, metadata_issues = _write_metadata(
        final_ids, metadata, output / "myoga_metadata.tsv"
    )
    final_source_by_id = {
        sample.sample_id: source_index
        for sample, source_index in zip(final_samples, final_source_indexes)
    }
    final_sample_by_id = {sample.sample_id: sample for sample in final_samples}
    for issue in metadata_issues:
        source_index = final_source_by_id[str(issue["sample_id"])]
        issue.update(
            {
                "path": str(final_sample_by_id[str(issue["sample_id"])].path.parent),
                "callable_loci": int(callable_counts[source_index]),
                "total_loci": int(assayed_counts[source_index]),
                "callable_fraction": _format_number(
                    float(callable_fractions[source_index]), 6
                ),
            }
        )
    exclusions.extend(metadata_issues)

    call_fields = [
        "sample_id", "locus_id", "present", "repeat_count", "repeat_count_raw",
        "product_size_bp", "read_depth", "primary_read_depth", "mean_coverage",
        "allele_confidence", "status", "evidence",
    ]
    available_fields = {
        field for sample in final_samples for row in sample.rows for field in row
    }
    call_fields = [field for field in call_fields if field in available_fields or field in {"sample_id", "locus_id"}]
    call_fields.extend(
        sorted(available_fields - set(call_fields), key=_natural_key)
    )
    long_rows = []
    profile_rows = []
    for sample, source_index in zip(final_samples, final_source_indexes):
        by_locus = rows_by_sample[sample.sample_id]
        for locus in locus_order:
            source = by_locus.get(locus, {})
            long_rows.append(
                {**{field: source.get(field, "") for field in call_fields}, "sample_id": sample.sample_id, "locus_id": locus}
            )
        profile_rows.append(
            {
                "sample_id": sample.sample_id,
                **{
                    locus: _format_allele(matrix[source_index, locus_index])
                    for locus_index, locus in enumerate(locus_order)
                },
            }
        )
    write_tsv(long_rows, output / "mlva_calls_long.tsv", call_fields)
    write_tsv(profile_rows, output / "mlva_profiles.tsv", ["sample_id", *locus_order])

    distance_rows = []
    for row_index, sample_id in enumerate(final_ids):
        distance_rows.append(
            {
                "sample_id": sample_id,
                **{
                    other: _format_number(float(final_matrix[row_index, column_index]))
                    for column_index, other in enumerate(final_ids)
                },
            }
        )
    write_tsv(
        distance_rows,
        output / "mlva_distance_matrix.tsv",
        ["sample_id", *final_ids],
    )

    tree_path = output / "mlva_nj.tree"
    if final_ids:
        _write_text_atomic(
            tree_path,
            neighbor_joining_tree_from_matrix(final_ids, final_matrix),
        )
    elif tree_path.exists():
        tree_path.unlink()

    combined_result: dict[str, Path | int | str] = {}
    if combined_markers:
        panel_loci = (
            read_loci_or_primers(loci_path, None)
            if loci_path is not None
            else []
        )
        final_repeat_matrix = (
            matrix[final_source_indexes, :]
            if final_source_indexes
            else np.empty((0, total_loci), dtype=np.float64)
        )
        combined_result = export_combined_markers(
            final_samples,
            locus_order,
            final_repeat_matrix,
            output / "myoga_metadata.tsv",
            output,
            loci=panel_loci,
            min_pairwise_loci=min_pairwise_loci,
            min_pairwise_fraction=min_pairwise_fraction,
            snp_weight=snp_weight,
            repeat_weight=repeat_weight,
            threads=threads,
            mafft_bin=mafft_bin,
            raxml_ng_bin=raxml_ng_bin,
            raxml_model=raxml_model,
            force=force,
        )

    used_rows = []
    for sample, source_index in zip(final_samples, final_source_indexes):
        metadata_found, coordinate_status = metadata_status[sample.sample_id]
        used_rows.append(
            {
                "sample_id": sample.sample_id,
                "path": str(sample.path.parent),
                "callable_loci": int(callable_counts[source_index]),
                "total_loci": int(assayed_counts[source_index]),
                "callable_fraction": _format_number(float(callable_fractions[source_index]), 6),
                "metadata_found": "yes" if metadata_found else "no",
                "coordinates_valid": "yes" if coordinate_status == "VALID" else "no",
                "coordinate_status": coordinate_status,
            }
        )
    write_tsv(used_rows, output / "samples_used.tsv", USED_FIELDS)
    exclusions.sort(key=lambda row: (_natural_key(str(row["sample_id"])), row["scope"], row["reason"], row["path"]))
    write_tsv(exclusions, output / "samples_excluded.tsv", EXCLUDED_FIELDS)

    upper = overlap_matrix[np.triu_indices(len(threshold_ids), k=1)]
    callable_values = [int(value) for value in callable_counts]
    assayed_values = [int(value) for value in assayed_counts]
    effective_values = [int(value) for value in required_counts]
    discovery_failure_reasons = {
        "FAILED_BATCH_SAMPLE",
        "MISSING_CALLS_FILE",
        "MALFORMED_RESULTS",
        "DUPLICATE_SAMPLE_ID",
    }
    failed_discovery_ids = {
        str(row["sample_id"])
        for row in exclusions
        if row["scope"] == "tree" and row["reason"] in discovery_failure_reasons
    }
    tree_excluded_ids = {
        str(row["sample_id"]) for row in exclusions if row["scope"] == "tree"
    }
    summary_values: list[tuple[str, object]] = [
        ("result_directories_discovered", discovered),
        ("successful_mlvamaps_samples", len(samples)),
        ("failed_or_incomplete_samples", len(failed_discovery_ids)),
        ("samples_with_metadata", sum(row["metadata_found"] == "yes" for row in used_rows)),
        ("samples_without_metadata", sum(row["metadata_found"] == "no" for row in used_rows)),
        ("samples_with_coordinates", sum(row["coordinates_valid"] == "yes" for row in used_rows)),
        ("samples_without_coordinates", sum(row["coordinates_valid"] == "no" for row in used_rows)),
        ("total_mlva_loci", total_loci),
        ("callable_fraction_denominator", "sample_assayed_loci"),
        ("pairwise_fraction_denominator", "shared_assayed_loci"),
        ("samples_passing_callable_threshold", len(threshold_samples)),
        ("samples_excluded_from_tree", len(tree_excluded_ids)),
        ("final_tree_samples", len(final_ids)),
        ("pairwise_comparisons", len(threshold_ids) * (len(threshold_ids) - 1) // 2),
        ("supported_pairwise_comparisons", int(np.isfinite(selected_matrix[np.triu_indices(len(threshold_ids), k=1)]).sum())),
        ("chosen_distance_metric", distance),
        ("combined_marker_export", "yes" if combined_markers else "no"),
        ("combined_marker_loci_built", combined_result.get("loci_built", 0)),
        ("combined_marker_tree_samples", combined_result.get("tree_samples", 0)),
        ("minimum_callable_fraction", min_callable_fraction),
        ("minimum_callable_loci", min_callable_loci),
        (
            "effective_minimum_callable_loci",
            effective_values[0]
            if effective_values and len(set(effective_values)) == 1
            else "sample-specific",
        ),
        ("minimum_effective_callable_loci", min(effective_values, default=0)),
        ("maximum_effective_callable_loci", max(effective_values, default=0)),
        ("minimum_pairwise_loci", min_pairwise_loci),
        ("minimum_pairwise_fraction", min_pairwise_fraction),
        ("minimum_callable_loci_observed", min(callable_values, default=0)),
        ("median_callable_loci_observed", statistics.median(callable_values) if callable_values else 0),
        ("maximum_callable_loci_observed", max(callable_values, default=0)),
        ("minimum_assayed_loci_observed", min(assayed_values, default=0)),
        ("median_assayed_loci_observed", statistics.median(assayed_values) if assayed_values else 0),
        ("maximum_assayed_loci_observed", max(assayed_values, default=0)),
        ("minimum_pairwise_loci_observed", int(upper.min()) if upper.size else 0),
        ("median_pairwise_loci_observed", float(np.median(upper)) if upper.size else 0),
        ("maximum_pairwise_loci_observed", int(upper.max()) if upper.size else 0),
    ]
    write_tsv(
        [{"metric": metric, "value": value} for metric, value in summary_values],
        output / "export_summary.tsv",
        ["metric", "value"],
    )
    summary_text = [
        "MLVAmaps MYOGA export",
        "======================",
        "",
        *[f"{metric}: {value}" for metric, value in summary_values],
        "",
        "The tree is a neighbor-joining MLVA relatedness tree, not a whole-genome phylogeny.",
        "Distances use only exact repeat counts callable in both samples; missing loci are not imputed.",
    ]
    if combined_markers:
        summary_text.extend(
            [
                "",
                "The optional combined-marker tree is a repeat-aware multilocus marker tree, not a whole-genome phylogeny.",
                "Its SNP component comes from VNTR-masked accepted amplicons and is averaged over shared callable loci.",
            ]
        )
    _write_text_atomic(output / "export_summary.txt", "\n".join(summary_text) + "\n")

    return {
        "outdir": output,
        "metadata": output / "myoga_metadata.tsv",
        "profiles": output / "mlva_profiles.tsv",
        "calls_long": output / "mlva_calls_long.tsv",
        "pairwise": output / "mlva_pairwise_distances.tsv",
        "distance_matrix": output / "mlva_distance_matrix.tsv",
        "tree": tree_path if final_ids else "",
        "samples_used": output / "samples_used.tsv",
        "samples_excluded": output / "samples_excluded.tsv",
        "summary": output / "export_summary.tsv",
        "tree_samples": len(final_ids),
        "combined_marker_tree": combined_result.get("tree", ""),
        "combined_marker_distance_matrix": combined_result.get("distance_matrix", ""),
        "combined_marker_metadata": combined_result.get("metadata", ""),
        "combined_marker_locus_status": combined_result.get("locus_status", ""),
    }
