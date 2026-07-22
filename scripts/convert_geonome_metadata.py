#!/usr/bin/env python3
"""Convert Geonome reference-build metadata to MLVAMaps metadata TSV."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Iterable


OUTPUT_FIELDS = [
    "reference_id",
    "collection_date",
    "latitude",
    "longitude",
    "location",
    "source",
]


def _clean(value: object) -> str:
    return "" if value is None else str(value).strip()


def _delimiter(path: Path) -> str:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        first_line = handle.readline()
    return "\t" if "\t" in first_line else ","


def _read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise ValueError(f"metadata file does not exist: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=_delimiter(path))
        if not reader.fieldnames:
            raise ValueError(f"metadata has no header: {path}")
        fields = [_clean(field) for field in reader.fieldnames]
        rows = [
            {_clean(key): _clean(value) for key, value in row.items() if key is not None}
            for row in reader
        ]
    return fields, rows


def _metadata_from_manifest(manifest_path: Path) -> Path:
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid JSON in Geonome manifest {manifest_path}: {exc}") from exc
    relative_path = manifest.get("artifacts", {}).get("normalized_metadata")
    if not relative_path:
        raise ValueError(
            f"Geonome manifest does not declare artifacts.normalized_metadata: {manifest_path}"
        )
    metadata_path = manifest_path.parent / str(relative_path)
    if not metadata_path.is_file():
        raise ValueError(
            f"normalized metadata declared by the Geonome manifest does not exist: {metadata_path}"
        )
    return metadata_path


def resolve_input(input_path: Path) -> Path:
    """Resolve a reference directory, manifest, or table to its metadata table."""
    if input_path.is_dir():
        manifest_path = input_path / "reference_manifest.json"
        if manifest_path.is_file():
            return _metadata_from_manifest(manifest_path)
        metadata_path = input_path / "normalized_metadata.tsv"
        if metadata_path.is_file():
            return metadata_path
        raise ValueError(
            f"directory is not a Geonome reference build (no reference_manifest.json "
            f"or normalized_metadata.tsv): {input_path}"
        )
    if input_path.name == "reference_manifest.json" or input_path.suffix.lower() == ".json":
        return _metadata_from_manifest(input_path)
    return input_path


def _first(row: dict[str, str], names: Iterable[str]) -> str:
    return next((row[name] for name in names if row.get(name)), "")


def _location(row: dict[str, str]) -> str:
    explicit = _first(row, ("location_raw", "location"))
    if explicit:
        return explicit
    parts: list[str] = []
    for field in ("county", "state", "country"):
        value = row.get(field, "")
        if value and value not in parts:
            parts.append(value)
    return ", ".join(parts)


def convert_rows(fields: list[str], rows: list[dict[str, str]]) -> list[dict[str, str]]:
    id_field = next(
        (field for field in ("genome_id", "reference_id", "accession") if field in fields),
        None,
    )
    if id_field is None:
        raise ValueError("Geonome metadata needs a genome_id or accession column")

    converted: list[dict[str, str]] = []
    seen: set[str] = set()
    for line_number, row in enumerate(rows, 2):
        reference_id = row.get(id_field, "")
        if not reference_id:
            raise ValueError(f"empty {id_field} on metadata line {line_number}")
        if any(character.isspace() for character in reference_id):
            raise ValueError(
                f"reference identifier contains whitespace on metadata line {line_number}: "
                f"{reference_id!r}"
            )
        if reference_id in seen:
            raise ValueError(f"duplicate reference identifier: {reference_id}")
        seen.add(reference_id)

        latitude = row.get("latitude", "")
        longitude = row.get("longitude", "")
        if bool(latitude) != bool(longitude):
            raise ValueError(
                f"reference {reference_id!r} has only one of latitude and longitude"
            )
        converted.append(
            {
                "reference_id": reference_id,
                "collection_date": _first(row, ("collection_date", "date")),
                "latitude": latitude,
                "longitude": longitude,
                "location": _location(row),
                "source": _first(
                    row,
                    ("isolation_source", "sample_type_normalized", "sample_type_raw", "source"),
                ),
            }
        )
    if not converted:
        raise ValueError("Geonome metadata contains no data rows")
    return converted


def convert_metadata(input_path: Path, output_path: Path) -> int:
    metadata_path = resolve_input(input_path)
    fields, rows = _read_rows(metadata_path)
    converted = convert_rows(fields, rows)
    if output_path.resolve() == metadata_path.resolve():
        raise ValueError("output path must differ from the input metadata path")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(converted)
    return len(converted)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert normalized metadata from a Geonome reference build into the "
            "reference_metadata.tsv format used by MLVAMaps."
        )
    )
    parser.add_argument(
        "input",
        type=Path,
        help=(
            "Geonome reference directory, reference_manifest.json, or "
            "normalized_metadata.tsv"
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("reference_metadata.tsv"),
        help="output TSV (default: reference_metadata.tsv)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        count = convert_metadata(args.input, args.output)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(f"Wrote {count} MLVAMaps metadata rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
