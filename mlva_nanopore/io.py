from __future__ import annotations

import csv
import gzip
from pathlib import Path
from typing import Iterable, Iterator, Optional

from .models import Locus, ReadRecord


def open_text(path: str | Path, mode: str = "rt"):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, mode)
    return path.open(mode, newline="")


def read_fastq(path: str | Path) -> Iterator[ReadRecord]:
    with open_text(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            plus = handle.readline()
            quality = handle.readline().strip()
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"Malformed FASTQ record near {header.strip()!r}")
            read_id = header[1:].strip().split()[0]
            yield ReadRecord(read_id, sequence.upper(), quality)


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
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
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
