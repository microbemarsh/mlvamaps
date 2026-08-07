from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable

from .io import open_text


METADATA_ALIASES = {
    "sample_id": ("sample_id", "run_accession", "sra_run", "accession"),
    "run_accession": ("run_accession", "sra_run", "accession", "sample_id"),
    "biosample": ("biosample", "biosample_accession", "metadata_id"),
    "collection_date": ("collection_date", "date", "sample_date", "isolation_date"),
    "latitude": ("latitude", "lat", "decimal_latitude", "decimal_lat"),
    "longitude": ("longitude", "lon", "lng", "decimal_longitude", "decimal_lon"),
    "location": ("location", "geo_loc_name", "geographic_location", "place"),
    "country": ("country",),
    "host": ("host",),
    "isolation_source": ("isolation_source", "source", "sample_type"),
    "study_accession": ("study_accession", "bioproject", "study"),
}

MYOGA_SAMPLE_FIELDS = [
    "genome_id",
    "sample_id",
    "run_accession",
    "biosample",
    "collection_date",
    "latitude",
    "longitude",
    "location",
    "country",
    "host",
    "isolation_source",
    "study_accession",
    "best_profile_id",
    "profile_distance",
    "profile_confidence",
    "complete_loci",
    "total_loci",
    "mlva_fingerprint",
    "read_technology",
]


def _delimiter(path: Path, sample: str) -> str:
    first = sample.splitlines()[0] if sample.splitlines() else ""
    if "\t" in first:
        return "\t"
    if "," in first:
        return ","
    raise ValueError(f"Metadata table {path} must be comma- or tab-delimited")


def read_sample_metadata(path: str | Path | None) -> list[dict[str, str]]:
    if path is None:
        return []
    metadata_path = Path(path)
    with open_text(metadata_path, "rt") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=_delimiter(metadata_path, sample))
        if not reader.fieldnames:
            raise ValueError(f"Metadata table {metadata_path} has no header")
        rows = [
            {
                str(key).strip(): "" if value is None else str(value).strip()
                for key, value in row.items()
                if key is not None
            }
            for row in reader
        ]
    normalized = [normalize_metadata_row(row) for row in rows]
    sample_ids = [row["sample_id"] for row in normalized]
    if any(not sample_id for sample_id in sample_ids):
        raise ValueError(
            f"Metadata table {metadata_path} has a row without a sample identifier"
        )
    duplicates = sorted(
        sample_id for sample_id in set(sample_ids) if sample_ids.count(sample_id) > 1
    )
    if duplicates:
        raise ValueError(
            "Metadata sample identifiers must be unique: " + ", ".join(duplicates)
        )
    return normalized


def normalize_metadata_row(row: dict[str, str]) -> dict[str, str]:
    """Append stable aliases while retaining every original metadata column."""
    lower_to_original = {key.strip().lower(): key for key in row}
    normalized = dict(row)
    for canonical, aliases in METADATA_ALIASES.items():
        value = ""
        for alias in aliases:
            original = lower_to_original.get(alias)
            if original is not None and row.get(original, "").strip():
                value = row[original].strip()
                break
        normalized[canonical] = value
    return normalized


def metadata_by_sample(rows: Iterable[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {str(row["sample_id"]): dict(row) for row in rows}


def myoga_sample_row(
    sample_id: str,
    metadata: dict[str, str] | None,
    summary: dict,
    fingerprint: dict[str, str],
    total_loci: int,
) -> dict[str, str | int]:
    """Create a MYOGA row whose ID exactly matches generated tree tip labels."""
    metadata = metadata or {}
    fingerprint_text = ";".join(
        f"{locus}={allele}"
        for locus, allele in fingerprint.items()
        if locus != "sample_id" and allele not in ("", None)
    )
    row: dict[str, str | int] = dict(metadata)
    row.update(
        {
            "genome_id": sample_id,
            "sample_id": sample_id,
            "run_accession": metadata.get("run_accession", sample_id),
            "biosample": metadata.get("biosample", ""),
            "collection_date": metadata.get("collection_date", ""),
            "latitude": metadata.get("latitude", ""),
            "longitude": metadata.get("longitude", ""),
            "location": metadata.get("location", ""),
            "country": metadata.get("country", ""),
            "host": metadata.get("host", ""),
            "isolation_source": metadata.get("isolation_source", ""),
            "study_accession": metadata.get("study_accession", ""),
            "best_profile_id": summary.get("best_profile_id", ""),
            "profile_distance": summary.get("best_profile_distance", ""),
            "profile_confidence": summary.get("profile_confidence", ""),
            "complete_loci": summary.get("complete_loci", 0),
            "total_loci": total_loci,
            "mlva_fingerprint": fingerprint_text,
            "read_technology": summary.get("read_technology", "illumina"),
        }
    )
    return row


def write_csv(rows: Iterable[dict], path: str | Path, preferred_fields: list[str]) -> Path:
    rows = list(rows)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    original_fields = sorted(
        {
            str(key)
            for row in rows
            for key in row
            if str(key) not in preferred_fields
        }
    )
    fields = preferred_fields + original_fields
    temporary = output.with_name(output.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output)
    return output
