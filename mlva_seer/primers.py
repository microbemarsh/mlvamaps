from __future__ import annotations

import csv
import re
from pathlib import Path

from .io import open_text
from .models import Locus


LEGACY_NAME_RE = re.compile(r"^(?P<name>.+)_(?P<repeat_bp>\d+)bp_(?P<amplicon_bp>\d+)bp_(?P<units>\d+)[Uu]$")


def _sniff_delimiter(sample: str) -> str | None:
    if "\t" in sample:
        return "\t"
    if "," in sample:
        return ","
    return None


def _int_value(value: str | None) -> int | None:
    if value in (None, ""):
        return None
    return int(float(value))


def _locus_from_primer_row(
    locus_id: str,
    forward: str,
    reverse: str,
    repeat_unit_length_bp: str | None = None,
    expected_product_size_bp: str | None = None,
    nominal_repeat_units: str | None = None,
) -> Locus:
    repeat_motif = ""
    expected_min = 0
    expected_max = 100
    expected_min_bp = 0
    expected_max_bp = 100000
    repeat_bp = _int_value(repeat_unit_length_bp)
    amplicon_bp = _int_value(expected_product_size_bp)
    units = _int_value(nominal_repeat_units)
    match = LEGACY_NAME_RE.match(locus_id)
    if match and repeat_bp is None:
        repeat_bp = int(match.group("repeat_bp"))
    if match and amplicon_bp is None:
        amplicon_bp = int(match.group("amplicon_bp"))
    if match and units is None:
        units = int(match.group("units"))
    if repeat_bp:
        repeat_motif = "N" * repeat_bp
    if units is not None:
        expected_min = max(0, units - 10)
        expected_max = units + 10
    if repeat_bp and amplicon_bp:
        expected_min_bp = max(len(forward) + len(reverse), amplicon_bp - repeat_bp * 10)
        expected_max_bp = amplicon_bp + repeat_bp * 10
    return Locus(
        locus_id=locus_id,
        forward_primer=forward.upper(),
        reverse_primer=reverse.upper(),
        repeat_motif=repeat_motif,
        expected_min_repeats=expected_min,
        expected_max_repeats=expected_max,
        expected_amplicon_min_bp=expected_min_bp,
        expected_amplicon_max_bp=expected_max_bp,
        repeat_unit_length_bp=repeat_bp or 0,
        expected_product_size_bp=amplicon_bp or 0,
        nominal_repeat_units=units or 0,
    )


def read_primer_pairs(path: str | Path) -> list[Locus]:
    """Read primer-pair files with locus, forward, reverse columns.

    Supports comma-delimited CSV, tab-delimited TSV, and whitespace-delimited
    legacy MLVA_finder primer files such as:

    locus_id forward_primer reverse_primer
    """
    with open_text(path, "rt") as handle:
        content = handle.read()
    lines = [line.strip() for line in content.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        return []

    delimiter = _sniff_delimiter(lines[0])
    rows: list[list[str]]
    if delimiter:
        reader = csv.reader(lines, delimiter=delimiter)
        rows = [[cell.strip() for cell in row] for row in reader]
    else:
        rows = [line.split() for line in lines]

    header = [cell.lower() for cell in rows[0]]
    header_aliases = {
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
    normalized_header = [header_aliases.get(cell, cell) for cell in header]
    has_header = {"locus_id", "forward_primer", "reverse_primer"}.issubset(normalized_header)

    loci = []
    data_rows = rows[1:] if has_header else rows
    if has_header:
        idx = {name: normalized_header.index(name) for name in ("locus_id", "forward_primer", "reverse_primer")}
        optional_idx = {
            name: normalized_header.index(name)
            for name in ("repeat_unit_length_bp", "expected_product_size_bp", "nominal_repeat_units")
            if name in normalized_header
        }
        for row in data_rows:
            if len(row) <= max(idx.values()):
                continue
            optional = {
                name: row[position] if len(row) > position else None
                for name, position in optional_idx.items()
            }
            loci.append(
                _locus_from_primer_row(
                    row[idx["locus_id"]],
                    row[idx["forward_primer"]],
                    row[idx["reverse_primer"]],
                    optional.get("repeat_unit_length_bp"),
                    optional.get("expected_product_size_bp"),
                    optional.get("nominal_repeat_units"),
                )
            )
        return loci

    for row in data_rows:
        if len(row) < 3:
            continue
        loci.append(_locus_from_primer_row(row[0], row[1], row[2]))
    return loci


def read_loci_or_primers(loci_path: str | Path | None = None, primers_path: str | Path | None = None) -> list[Locus]:
    if loci_path:
        from .io import read_loci

        return read_loci(loci_path)
    if primers_path:
        return read_primer_pairs(primers_path)
    raise ValueError("Provide either a loci TSV with --loci or a primer-pair file with --primers.")
