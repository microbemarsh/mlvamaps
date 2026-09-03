"""Shared generation and persistence of candidate MLVA allele contexts."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Iterable

from .io import read_fasta, write_fasta, write_tsv
from .models import Locus
from .repeat_calibration import (
    assembly_equivalent_product_allele,
    repeat_unit_length,
)
from .sequence import revcomp


CONTEXT_SCHEMA_VERSION = "2.0"
CONTEXT_FIELDS = [
    "schema_version", "candidate_id", "locus_id", "repeat_count",
    "repeat_unit_length", "repeat_start", "repeat_end", "left_flank",
    "right_flank", "reference_accession", "taxon_id", "taxon_name",
    "reference_sequence_id", "reference_id", "sequence_sha256",
    "context_length_bp", "expected_product_size_bp", "repeat_motif",
    "reference_contig", "reference_start", "reference_end", "strand",
    "source", "provenance_count",
]


@dataclass(frozen=True)
class CandidateContext:
    candidate_id: str
    locus_id: str
    sequence: str
    repeat_count: int | float | None
    repeat_unit_length: int
    repeat_start: int
    repeat_end: int
    left_flank: str
    right_flank: str
    reference_accession: str = ""
    taxon_id: str = ""
    taxon_name: str = ""
    reference_sequence_id: str = ""
    reference_id: str = ""
    expected_product_size_bp: int = 0
    repeat_motif: str = ""
    reference_contig: str = ""
    reference_start: int | None = None
    reference_end: int | None = None
    strand: str = "+"
    source: str = ""
    provenance_count: int = 1

    # Compatibility names used by the original Illumina implementation.
    @property
    def context_id(self) -> str:
        return self.candidate_id

    @property
    def expected_repeat_count(self) -> int | float | None:
        return self.repeat_count

    @property
    def repeat_unit_length_bp(self) -> int:
        return self.repeat_unit_length

    @property
    def upstream_flank_bp(self) -> int:
        return self.repeat_start

    @property
    def downstream_flank_bp(self) -> int:
        return len(self.sequence) - self.repeat_end

    def row(self) -> dict[str, object]:
        return {
            "schema_version": CONTEXT_SCHEMA_VERSION,
            "candidate_id": self.candidate_id,
            "locus_id": self.locus_id,
            "repeat_count": "" if self.repeat_count is None else self.repeat_count,
            "repeat_unit_length": self.repeat_unit_length,
            "repeat_start": self.repeat_start,
            "repeat_end": self.repeat_end,
            "left_flank": self.left_flank,
            "right_flank": self.right_flank,
            "reference_accession": self.reference_accession,
            "taxon_id": self.taxon_id,
            "taxon_name": self.taxon_name,
            "reference_sequence_id": self.reference_sequence_id,
            "reference_id": self.reference_id or self.reference_accession,
            "sequence_sha256": hashlib.sha256(self.sequence.encode()).hexdigest(),
            "context_length_bp": len(self.sequence),
            "expected_product_size_bp": self.expected_product_size_bp or len(self.sequence),
            "repeat_motif": self.repeat_motif,
            "reference_contig": self.reference_contig,
            "reference_start": "" if self.reference_start is None else self.reference_start,
            "reference_end": "" if self.reference_end is None else self.reference_end,
            "strand": self.strand,
            "source": self.source,
            "provenance_count": self.provenance_count,
        }


# Compatibility alias while downstream integrations migrate terminology.
LocusContext = CandidateContext


def _database_root(database: str | Path | None) -> Path | None:
    if not database:
        return None
    path = Path(database)
    for candidate in (path, path / "database"):
        if candidate.is_dir():
            return candidate
    return None


def _metadata(database: Path | None) -> dict[str, dict[str, str]]:
    if database is None or not (database / "reference_metadata.tsv").is_file():
        return {}
    with (database / "reference_metadata.tsv").open(newline="") as handle:
        return {
            str(row.get("reference_id", "")): row
            for row in csv.DictReader(handle, delimiter="\t")
        }


def _synthetic_product(locus: Locus) -> str:
    unit = repeat_unit_length(locus)
    motif = locus.repeat_motif
    if (
        not locus.forward_primer or not locus.reverse_primer or not motif or not unit
        or set(motif) - set("ACGT")
        or not locus.left_flank_sequence or not locus.right_flank_sequence
    ):
        return ""
    nominal = locus.nominal_repeat_units or max(
        locus.expected_min_repeats,
        min(locus.expected_max_repeats, round((locus.expected_min_repeats + locus.expected_max_repeats) / 2)),
    )
    repeated = (motif * math.ceil(max(nominal * unit, 1) / len(motif)))[: nominal * unit]
    return (
        locus.forward_primer + locus.left_flank_sequence + repeated
        + locus.right_flank_sequence + revcomp(locus.reverse_primer)
    )


def repeat_interval(sequence: str, locus: Locus) -> tuple[int, int]:
    left = sequence.find(locus.left_flank_sequence) if locus.left_flank_sequence else -1
    if left >= 0:
        start = left + len(locus.left_flank_sequence)
        end = sequence.find(locus.right_flank_sequence, start) if locus.right_flank_sequence else -1
        if end >= start:
            return start, end
    unit = repeat_unit_length(locus)
    nominal_bp = unit * locus.nominal_repeat_units
    if nominal_bp and locus.expected_product_size_bp:
        start = max(0, (len(sequence) - nominal_bp) // 2)
        return start, min(len(sequence), start + nominal_bp)
    start = len(locus.forward_primer) + len(locus.left_flank_sequence)
    end = len(sequence) - len(locus.right_flank_sequence) - len(locus.reverse_primer)
    return max(0, start), max(start, end)


def candidate_repeat_counts(
    locus: Locus,
    observed_states: Iterable[int | float | None],
    maximum: int = 100,
    expansion: int = 1,
) -> list[int]:
    if maximum < 1:
        raise ValueError("maximum candidate repeat count must be at least 1")
    if expansion < 0:
        raise ValueError("candidate expansion cannot be negative")
    expected_lower = max(0, int(locus.expected_min_repeats))
    expected_upper = min(maximum, int(locus.expected_max_repeats))
    states = set(range(expected_lower, expected_upper + 1))
    for observed in observed_states:
        if observed is None or not float(observed).is_integer():
            continue
        center = int(observed)
        for state in range(max(0, center - expansion), min(maximum, center + expansion) + 1):
            # Database observations are retained, but expansion remains local
            # and the explicit maximum prevents unbounded state generation.
            states.add(state)
    return sorted(states)


def _load_v2(database: Path, loci: list[Locus]) -> list[CandidateContext] | None:
    metadata_path = database / "candidate_metadata.tsv"
    fasta_path = database / "candidate_contexts.fasta"
    if not fasta_path.is_file():
        fasta_path = database / "candidate_contexts.fasta.gz"
    if not metadata_path.is_file() or not fasta_path.is_file():
        return None
    sequences = dict(read_fasta(fasta_path))
    with metadata_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    missing = set(CONTEXT_FIELDS) - set(rows[0] if rows else [])
    if missing:
        raise ValueError("Invalid candidate metadata; missing columns: " + ", ".join(sorted(missing)))
    locus_ids = {locus.locus_id for locus in loci}
    contexts = []
    for row in rows:
        candidate_id = row["candidate_id"]
        if row["locus_id"] not in locus_ids or candidate_id not in sequences:
            continue
        sequence = sequences[candidate_id]
        if hashlib.sha256(sequence.encode()).hexdigest() != row["sequence_sha256"]:
            raise ValueError(f"Candidate context hash mismatch for {candidate_id!r}")
        contexts.append(CandidateContext(
            candidate_id, row["locus_id"], sequence,
            float(row["repeat_count"]) if "." in row["repeat_count"] else int(row["repeat_count"]),
            int(row["repeat_unit_length"]), int(row["repeat_start"]), int(row["repeat_end"]),
            row["left_flank"], row["right_flank"], row["reference_accession"],
            row["taxon_id"], row["taxon_name"], row["reference_sequence_id"],
            row["reference_id"], int(row["expected_product_size_bp"] or 0),
            row["repeat_motif"], row["reference_contig"],
            int(row["reference_start"]) if row["reference_start"] else None,
            int(row["reference_end"]) if row["reference_end"] else None,
            row["strand"], row["source"], int(row["provenance_count"] or 1),
        ))
    return contexts


def _base_contexts(loci: list[Locus], database_path: str | Path | None) -> list[CandidateContext]:
    database = _database_root(database_path)
    if database:
        stored = _load_v2(database, loci)
        if stored is not None:
            return stored
    metadata = _metadata(database)
    contexts: list[CandidateContext] = []
    for locus in loci:
        records: list[tuple[str, str]] = []
        if database:
            for suffix in (".fasta.gz", ".fasta", ".fa.gz", ".fa"):
                path = database / f"{locus.locus_id}{suffix}"
                if path.is_file():
                    records = list(read_fasta(path))
                    break
        if not records:
            sequence = _synthetic_product(locus)
            if sequence:
                records = [("panel", sequence)]
        for index, (reference_id, sequence) in enumerate(records, 1):
            sequence = sequence.upper()
            start, end = repeat_interval(sequence, locus)
            repeat = assembly_equivalent_product_allele(locus, len(sequence))[1]
            if repeat is None and repeat_unit_length(locus):
                repeat = (end - start) / repeat_unit_length(locus)
            meta = metadata.get(reference_id.split()[0], {})
            contexts.append(CandidateContext(
                f"base{len(contexts)+1:07d}", locus.locus_id, sequence, repeat,
                repeat_unit_length(locus), start, end, locus.left_flank_sequence,
                locus.right_flank_sequence, reference_id.split()[0],
                str(meta.get("taxon_id") or meta.get("taxid") or ""),
                str(meta.get("taxon_name") or meta.get("organism_name") or ""),
                reference_id, reference_id, len(sequence), locus.repeat_motif,
                source="database_amplicon" if database else "panel_synthetic",
            ))
    if not contexts:
        raise ValueError(
            "Candidate contexts require a current database or a rich panel with "
            "flanks, a concrete repeat motif, and repeat bounds."
        )
    return contexts


def generate_candidate_contexts(
    loci: list[Locus],
    database_path: str | Path | None = None,
    maximum: int = 100,
    expansion: int = 1,
) -> list[CandidateContext]:
    bases = _base_contexts(loci, database_path)
    if bases and all(context.candidate_id.startswith("candidate") for context in bases):
        return bases
    loci_by_id = {locus.locus_id: locus for locus in loci}
    observed: dict[str, list[int | float | None]] = defaultdict(list)
    for context in bases:
        observed[context.locus_id].append(context.repeat_count)
    collapsed: dict[tuple[str, str, int], dict[str, object]] = {}
    for base in bases:
        locus = loci_by_id[base.locus_id]
        motif = locus.repeat_motif
        if not motif or set(motif) - set("ACGT") or base.repeat_end <= base.repeat_start:
            continue
        for count in candidate_repeat_counts(locus, observed[base.locus_id], maximum, expansion):
            repeat_sequence = (motif * math.ceil(max(count * base.repeat_unit_length, 1) / len(motif)))[
                : count * base.repeat_unit_length
            ]
            sequence = base.sequence[:base.repeat_start] + repeat_sequence + base.sequence[base.repeat_end:]
            key = (base.locus_id, sequence, count)
            item = collapsed.setdefault(key, {"base": base, "references": set(), "taxa": set(), "names": set(), "sequence_ids": set()})
            item["references"].add(base.reference_accession)
            item["taxa"].add(base.taxon_id)
            item["names"].add(base.taxon_name)
            item["sequence_ids"].add(base.reference_sequence_id)
    contexts = []
    ordered = sorted(
        collapsed.items(), key=lambda item: (item[0][0], float(item[0][2]), item[0][1])
    )
    for index, ((locus_id, sequence, count), item) in enumerate(ordered, 1):
        base = item["base"]
        contexts.append(replace(
            base,
            candidate_id=f"candidate{index:07d}", sequence=sequence, repeat_count=count,
            repeat_end=base.repeat_start + count * base.repeat_unit_length,
            reference_accession=";".join(sorted(item["references"] - {""})),
            reference_id=";".join(sorted(item["references"] - {""})),
            taxon_id=";".join(sorted(item["taxa"] - {""})),
            taxon_name=";".join(sorted(item["names"] - {""})),
            reference_sequence_id=";".join(sorted(item["sequence_ids"] - {""})),
            expected_product_size_bp=len(sequence), provenance_count=max(len(item["references"] - {""}), 1),
        ))
    if not contexts:
        raise ValueError("Candidate contexts do not encode usable discrete repeat-count states")
    return contexts


def write_candidate_contexts(contexts: list[CandidateContext], directory: str | Path) -> dict[str, Path]:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / "candidate_contexts.fasta"
    metadata = directory / "candidate_metadata.tsv"
    provenance = directory / "candidate_provenance.json"
    write_fasta(((context.candidate_id, context.sequence) for context in contexts), fasta)
    write_tsv((context.row() for context in contexts), metadata, CONTEXT_FIELDS)
    locus_ranges: dict[str, list[int | float]] = defaultdict(list)
    for context in contexts:
        if context.repeat_count is not None:
            locus_ranges[context.locus_id].append(context.repeat_count)
    provenance.write_text(json.dumps({
        "schema_version": CONTEXT_SCHEMA_VERSION,
        "candidate_count": len(contexts),
        "generation": {
            "states": "panel expected range plus observed database states and bounded local expansion",
            "duplicate_normalization": "identical locus/sequence/repeat contexts collapsed with provenance_count",
            "locus_repeat_ranges": {
                locus_id: [min(states, key=float), max(states, key=float)]
                for locus_id, states in sorted(locus_ranges.items()) if states
            },
        },
        "candidate_signature": hashlib.sha256(
            "\n".join(f"{c.candidate_id}\t{c.locus_id}\t{c.repeat_count}\t{c.sequence}" for c in contexts).encode()
        ).hexdigest(),
    }, indent=2, sort_keys=True) + "\n")
    return {"fasta": fasta, "metadata": metadata, "provenance": provenance}


# Compatibility API for the former short-read-only context module.
def load_locus_contexts(loci: list[Locus], database_path: str | Path | None = None) -> list[CandidateContext]:
    return _base_contexts(loci, database_path)


def expand_candidate_contexts(
    contexts: list[CandidateContext], loci: list[Locus], maximum: int = 100
) -> list[CandidateContext]:
    # Expand supplied contexts using the same bounded logic without a database.
    loci_by_id = {locus.locus_id: locus for locus in loci}
    observed: dict[str, list[int | float | None]] = defaultdict(list)
    for context in contexts:
        observed[context.locus_id].append(context.repeat_count)
    expanded = []
    seen = set()
    for context in contexts:
        locus = loci_by_id[context.locus_id]
        for count in candidate_repeat_counts(locus, observed[context.locus_id], maximum, 0):
            sequence = context.sequence[:context.repeat_start] + locus.repeat_motif * count + context.sequence[context.repeat_end:]
            key = (context.locus_id, sequence, count)
            if key in seen:
                continue
            seen.add(key)
            expanded.append(replace(
                context, candidate_id=f"candidate{len(expanded)+1:07d}", sequence=sequence,
                repeat_count=count, repeat_end=context.repeat_start + len(locus.repeat_motif) * count,
                expected_product_size_bp=len(sequence),
            ))
    if not expanded:
        raise ValueError("Candidate contexts do not encode usable discrete repeat-count states")
    return expanded