from __future__ import annotations

import csv
import gzip
import shutil
from pathlib import Path
from itertools import zip_longest
from typing import Iterable, Iterator, Optional

import pysam

from .models import Locus, ReadPair, ReadRecord


def open_text(path: str | Path, mode: str = "rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode, newline="")


def read_fastq(path: str | Path) -> Iterator[ReadRecord]:
    try:
        with pysam.FastxFile(str(path), persist=False) as handle:
            for record in handle:
                if record.quality is None:
                    raise ValueError(f"Expected FASTQ qualities for record {record.name!r}")
                if len(record.sequence) != len(record.quality):
                    raise ValueError(
                        f"FASTQ sequence and quality lengths differ for {record.name!r}"
                    )
                yield ReadRecord(record.name, record.sequence.upper(), record.quality)
    except (OSError, ValueError) as exc:
        raise ValueError(f"Malformed FASTQ {str(path)!r}: {exc}") from exc


def normalize_read_id(read_id: str) -> tuple[str, int | None]:
    """Return a stable molecule ID and an optional mate number.

    Pysam exposes the first whitespace-delimited FASTQ token as ``name``.  The
    slash convention is normalized here; CASAVA's ``1:N:...``/``2:N:...``
    annotation is retained when callers pass a complete header.
    """
    tokens = read_id.strip().lstrip("@").split()
    if not tokens:
        raise ValueError("FASTQ record has an empty read ID")
    primary = tokens[0]
    mate: int | None = None
    if primary.endswith("/1") or primary.endswith("/2"):
        mate = int(primary[-1])
        primary = primary[:-2]
    elif len(tokens) > 1 and len(tokens[1]) >= 2 and tokens[1][0] in "12" and tokens[1][1] == ":":
        mate = int(tokens[1][0])
    return primary, mate


def read_fastq_pairs(
    reads1_path: str | Path,
    reads2_path: str | Path | None = None,
) -> Iterator[ReadPair]:
    """Stream separate paired FASTQ files with strict name/count validation."""
    if reads2_path is None:
        for read in read_fastq(reads1_path):
            molecule_id, mate = normalize_read_id(read.read_id)
            if mate == 2:
                raise ValueError(
                    f"Single-end FASTQ record {read.read_id!r} is labelled as mate 2"
                )
            yield ReadPair(molecule_id, read, None)
        return

    for record_number, (read1, read2) in enumerate(
        zip_longest(read_fastq(reads1_path), read_fastq(reads2_path)), start=1
    ):
        if read1 is None or read2 is None:
            shorter = reads1_path if read1 is None else reads2_path
            raise ValueError(
                f"Paired FASTQ files have different record counts; {shorter!s} "
                f"ended before pair {record_number}"
            )
        id1, mate1 = normalize_read_id(read1.read_id)
        id2, mate2 = normalize_read_id(read2.read_id)
        if id1 != id2:
            raise ValueError(
                f"Paired FASTQ IDs differ at pair {record_number}: "
                f"{read1.read_id!r} versus {read2.read_id!r}"
            )
        if mate1 not in (None, 1) or mate2 not in (None, 2):
            raise ValueError(
                f"Invalid mate labels at pair {record_number}: "
                f"{read1.read_id!r}, {read2.read_id!r}"
            )
        yield ReadPair(id1, read1, read2)


def read_fasta(path: str | Path) -> Iterator[tuple[str, str]]:
    with pysam.FastxFile(str(path), persist=False) as handle:
        for record in handle:
            yield record.name, record.sequence.upper()


def write_fastq(reads: Iterable[ReadRecord], path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        for read in reads:
            quality = read.quality if read.quality is not None else "I" * len(read.sequence)
            handle.write(f"@{read.read_id}\n{read.sequence}\n+\n{quality}\n")


def _int_or_none(value: str | None) -> Optional[int]:
    if value is None or value == "":
        return None
    return int(float(value))


def _int_or_default(value: str | None, default: int) -> int:
    parsed = _int_or_none(value)
    return default if parsed is None else parsed


def read_loci(path: str | Path) -> list[Locus]:
    loci = []
    with open_text(path, "rt") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        first_line = sample.splitlines()[0] if sample.splitlines() else ""
        delimiter = "," if "," in first_line and "\t" not in first_line else "\t"
        reader = csv.DictReader(handle, delimiter=delimiter)
        aliases = {
            "name": "locus_id",
            "locus": "locus_id",
            "id": "locus_id",
            "forward": "forward_primer",
            "reverse": "reverse_primer",
            "fwd": "forward_primer",
            "rev": "reverse_primer",
            "repeat_bp": "repeat_unit_length_bp",
            "amplicon_bp": "expected_product_size_bp",
            "product_bp": "expected_product_size_bp",
            "units": "nominal_repeat_units",
        }
        normalized_fields = {
            aliases.get(str(field).strip().lower(), str(field).strip().lower())
            for field in (reader.fieldnames or [])
        }
        if "locus_id" not in normalized_fields:
            fields = ", ".join(str(field) for field in (reader.fieldnames or []))
            raise ValueError(
                f"Loci table needs a locus_id, name, locus, or id column; "
                f"found: {fields or '(no header)'}"
            )
        for raw_row in reader:
            row = {
                aliases.get(str(key).strip().lower(), str(key).strip().lower()): value
                for key, value in raw_row.items()
                if key is not None
            }
            loci.append(
                Locus(
                    locus_id=row["locus_id"],
                    chrom_or_contig=row.get("chrom_or_contig", ""),
                    start=_int_or_none(row.get("start")),
                    end=_int_or_none(row.get("end")),
                    forward_primer=row.get("forward_primer", "").upper(),
                    reverse_primer=row.get("reverse_primer", "").upper(),
                    left_flank_sequence=row.get("left_flank_sequence", "").upper(),
                    right_flank_sequence=row.get("right_flank_sequence", "").upper(),
                    repeat_motif=row.get("repeat_motif", "").upper(),
                    expected_min_repeats=_int_or_default(row.get("expected_min_repeats"), 0),
                    expected_max_repeats=_int_or_default(row.get("expected_max_repeats"), 100),
                    expected_amplicon_min_bp=_int_or_default(row.get("expected_amplicon_min_bp"), 0),
                    expected_amplicon_max_bp=_int_or_default(row.get("expected_amplicon_max_bp"), 100000),
                    pool_id=row.get("pool_id", ""),
                    repeat_unit_length_bp=_int_or_default(row.get("repeat_unit_length_bp"), 0),
                    expected_product_size_bp=_int_or_default(row.get("expected_product_size_bp"), 0),
                    nominal_repeat_units=_int_or_default(row.get("nominal_repeat_units"), 0),
                )
            )
    return loci


def read_profiles(path: str | Path | None) -> list[dict[str, str]]:
    if not path:
        return []
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(rows: Iterable[dict], path: str | Path, fieldnames: list[str]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_fasta(records: Iterable[tuple[str, str]], path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for idx in range(0, len(sequence), 80):
                handle.write(sequence[idx : idx + 80] + "\n")


def gzip_output_file(path: str | Path) -> Path:
    """Replace one generated sequence file with a gzip-compressed copy."""
    source = Path(path)
    if source.suffix.lower() == ".gz":
        return source
    destination = Path(f"{source}.gz")
    temporary = destination.with_name(f".{destination.name}.tmp")
    with source.open("rb") as input_handle, gzip.open(temporary, "wb") as output_handle:
        shutil.copyfileobj(input_handle, output_handle)
    temporary.replace(destination)
    source.unlink()
    return destination
