"""Observed reference sequences, metadata, and assembly provenance."""
from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path

from .sequence import revcomp

_FASTA_SUFFIXES = (".fa", ".fas", ".fasta", ".fna", ".ffn")


def _is_fasta_path(path: Path) -> bool:
    name = path.name.lower()
    return any(
        name.endswith(suffix) or name.endswith(f"{suffix}.gz")
        for suffix in _FASTA_SUFFIXES
    )


def _fasta_stem(path: Path) -> str:
    name = path.name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in _FASTA_SUFFIXES:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


REFERENCE_ASSEMBLY_FIELDS = ["reference_id", "assembly_file", "assembly_sha256"]


def _read_fasta(path: str | Path) -> list[tuple[str, str]]:
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    records: list[tuple[str, str]] = []
    name: str | None = None
    sequence: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence).upper()))
                name = line[1:].split()[0]
                sequence = []
            elif name is None:
                raise ValueError(f"Sequence appeared before a FASTA header in {path}")
            else:
                sequence.append(line)
    if name is not None:
        records.append((name, "".join(sequence).upper()))
    return records


def canonical_assembly_digest(path: str | Path) -> str:
    """Hash assembly sequence content independent of headers/order/orientation."""
    canonical_contigs = sorted(
        min(sequence, revcomp(sequence))
        for _name, sequence in _read_fasta(path)
    )
    payload = "".join(f"{len(sequence)}:{sequence};" for sequence in canonical_contigs)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _read_database_fasta(path: Path, locus_ids: set[str]) -> dict[str, list[tuple[str, str]]]:
    records = _read_fasta(path)
    if path.parent != path and _fasta_stem(path) in locus_ids:
        return {_fasta_stem(path): records}
    by_locus: dict[str, list[tuple[str, str]]] = {}
    for header, sequence in records:
        parts = header.split("|")
        matching = [part for part in parts if part in locus_ids]
        if len(matching) != 1:
            raise ValueError(
                f"Could not identify one panel locus in FASTA header {header!r}. "
                "Use LOCUS.fasta files or headers such as reference_id|LOCUS."
            )
        locus_id = matching[0]
        reference_id = next((part for part in parts if part != locus_id), "")
        if not reference_id:
            raise ValueError(f"No reference id was present in FASTA header {header!r}")
        by_locus.setdefault(locus_id, []).append((reference_id, sequence))
    return by_locus


def read_sequence_database(
    database_path: str | Path, locus_ids: set[str]
) -> dict[str, list[tuple[str, str]]]:
    """Read per-locus references from a directory, FASTA, or long-form TSV."""
    path = Path(database_path)
    if not path.exists():
        raise ValueError(f"Sequence database path does not exist: {path}")
    by_locus: dict[str, list[tuple[str, str]]] = {}
    if path.is_dir():
        fasta_paths = sorted(
            item for item in path.iterdir() if item.is_file() and _is_fasta_path(item)
        )
        if not fasta_paths:
            raise ValueError(f"Sequence database directory contains no FASTA files: {path}")
        for fasta_path in fasta_paths:
            locus_name = _fasta_stem(fasta_path)
            if locus_name not in locus_ids:
                continue
            by_locus[locus_name] = _read_fasta(fasta_path)
    elif _is_fasta_path(path):
        by_locus = _read_database_fasta(path, locus_ids)
    else:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"reference_id", "locus_id", "sequence"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise ValueError(
                    "Sequence database TSV requires reference_id, locus_id, and sequence columns"
                )
            for row in reader:
                locus_id = row["locus_id"]
                if locus_id in locus_ids:
                    by_locus.setdefault(locus_id, []).append(
                        (row["reference_id"], row["sequence"].upper())
                    )
    for locus_id, records in by_locus.items():
        names = [name for name, _sequence in records]
        if not records:
            raise ValueError(f"No reference sequences found for locus {locus_id!r}")
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate reference ids found for locus {locus_id!r}")
        if any(not name or not sequence for name, sequence in records):
            raise ValueError(f"Blank reference id or sequence found for locus {locus_id!r}")
    if not by_locus:
        raise ValueError("Sequence database contains no loci matching the supplied panel")
    return by_locus


def read_reference_metadata(path: str | Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    metadata_path = Path(path)
    if not metadata_path.exists():
        raise ValueError(f"Reference metadata path does not exist: {metadata_path}")
    with metadata_path.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        first_line = sample.splitlines()[0] if sample.splitlines() else ""
        delimiter = "\t" if "\t" in first_line else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames or "reference_id" not in reader.fieldnames:
            raise ValueError("Reference metadata requires a reference_id column")
        rows = list(reader)
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        reference_id = str(row.get("reference_id", "")).strip()
        if not reference_id:
            continue
        if reference_id in result:
            raise ValueError(f"Duplicate reference metadata for {reference_id!r}")
        result[reference_id] = {str(key): str(value or "") for key, value in row.items()}
    return result
