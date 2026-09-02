"""Competitive Bowtie2 mapping and conservative read-level VNTR inference.

The module deliberately separates locus detection from genotype resolution.  A
few flank-anchored molecules can establish presence without supplying enough
information to choose between adjacent repeat counts.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import statistics
import subprocess
from collections import defaultdict
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import pysam

from . import __version__
from .calling import assembly_equivalent_product_allele, repeat_unit_length
from .io import open_text, read_fasta, read_profiles, write_fasta, write_tsv
from .locus_measurement import measure_locus_product
from .models import Locus
from .profile_matching import (
    PROFILE_MATCH_LOCUS_FIELDS,
    build_fingerprint,
    match_profiles,
    profile_match_locus_rows,
)
from .sequence import revcomp


CONTEXT_SCHEMA_VERSION = "1.0"
CONTEXT_FIELDS = [
    "schema_version", "context_id", "locus_id", "reference_id", "taxon_id", "taxon_name",
    "sequence_sha256", "context_length_bp", "expected_product_size_bp",
    "repeat_motif", "repeat_unit_length_bp", "expected_repeat_count",
    "repeat_start", "repeat_end", "reference_contig", "reference_start",
    "reference_end", "strand", "upstream_flank_bp", "downstream_flank_bp",
    "source",
]

MAPPING_EVIDENCE_FIELDS = [
    "sample_id", "locus_id", "state", "repeat_count", "locus_length_bp",
    "confidence", "supporting_fragments", "proper_spanning_pairs",
    "junction_reads", "full_spanning_reads", "cigar_indel_reads",
    "mean_mapq", "median_mapq", "candidate_scores", "best_contexts",
    "context_taxa", "reason",
]


@dataclass(frozen=True)
class LocusContext:
    context_id: str
    locus_id: str
    sequence: str
    reference_id: str
    taxon_id: str = ""
    taxon_name: str = ""
    expected_product_size_bp: int = 0
    repeat_motif: str = ""
    repeat_unit_length_bp: int = 0
    expected_repeat_count: int | float | None = None
    repeat_start: int = 0
    repeat_end: int = 0
    reference_contig: str = ""
    reference_start: int | None = None
    reference_end: int | None = None
    strand: str = "+"
    upstream_flank_bp: int = 0
    downstream_flank_bp: int = 0
    source: str = ""

    def row(self) -> dict:
        row = asdict(self)
        row.pop("sequence")
        row["schema_version"] = CONTEXT_SCHEMA_VERSION
        row["sequence_sha256"] = hashlib.sha256(self.sequence.encode()).hexdigest()
        row["context_length_bp"] = len(self.sequence)
        return row


@dataclass(frozen=True)
class InsertSizeEstimate:
    median: float | None
    mean: float | None
    standard_deviation: float | None
    mad: float | None
    pairs_used: int


@dataclass
class _CandidateEvidence:
    locus_id: str
    repeat_count: int | float
    context_ids: set[str]
    molecules: set[str]
    spanning_pairs: set[str]
    junction_reads: set[str]
    full_reads: set[str]
    indel_reads: set[str]
    mapqs: list[int]
    alignment_scores: list[float]
    geometry_log_likelihoods: list[float]


def _metadata(database: Path | None) -> dict[str, dict[str, str]]:
    if database is None:
        return {}
    for root in (database, database / "database"):
        path = root / "reference_metadata.tsv"
        if path.exists():
            with path.open(newline="") as handle:
                return {
                    str(row.get("reference_id", "")): row
                    for row in csv.DictReader(handle, delimiter="\t")
                }
    return {}


def _database_root(database: str | Path | None, loci: list[Locus]) -> Path | None:
    if not database:
        return None
    path = Path(database)
    roots = (path / "database", path) if (path / "database").is_dir() else (path,)
    for root in roots:
        if (root / "mlva_contexts.tsv").exists() or (root / "mlva_contexts.fasta.gz").exists():
            return root
        if any((root / f"{locus.locus_id}.fasta.gz").exists() or
               (root / f"{locus.locus_id}.fasta").exists() for locus in loci):
            return root
    raise ValueError(f"Reference database has no MLVA data: {path}")


def _synthetic_product(locus: Locus) -> str:
    count = locus.nominal_repeat_units
    if not count:
        count = max(locus.expected_min_repeats, min(locus.expected_max_repeats, 1))
    motif = locus.repeat_motif
    if not motif or set(motif) == {"N"}:
        return ""
    return (
        locus.forward_primer + locus.left_flank_sequence + motif * count
        + locus.right_flank_sequence + revcomp(locus.reverse_primer)
    )


def _repeat_interval(sequence: str, locus: Locus) -> tuple[int, int]:
    measured = measure_locus_product(sequence, locus, source="reference_context")
    if measured.repeat_start is not None and measured.repeat_end is not None:
        return measured.repeat_start, measured.repeat_end
    unit = repeat_unit_length(locus)
    count = locus.nominal_repeat_units
    repeat_bp = unit * count
    if repeat_bp and locus.expected_product_size_bp:
        left = max(0, (len(sequence) - repeat_bp) // 2)
        return left, min(len(sequence), left + repeat_bp)
    left = len(locus.forward_primer) + len(locus.left_flank_sequence)
    right = len(sequence) - len(locus.right_flank_sequence) - len(locus.reverse_primer)
    return max(0, left), max(left, right)


def load_locus_contexts(
    loci: list[Locus], database_path: str | Path | None = None
) -> list[LocusContext]:
    """Load all candidate contexts without pre-classifying the sample taxon."""
    database = _database_root(database_path, loci)
    metadata = _metadata(database)
    stored_manifest = database / "mlva_contexts.tsv" if database else None
    stored_fasta = database / "mlva_contexts.fasta.gz" if database else None
    if database and not (
        stored_manifest and stored_manifest.is_file()
        and stored_fasta and stored_fasta.is_file()
    ):
        raise ValueError(
            "This reference database predates the Illumina context schema. "
            "Rebuild it with the current mlvamaps build-reference command."
        )
    if stored_manifest and stored_manifest.exists() and stored_fasta and stored_fasta.exists():
        sequences = dict(read_fasta(stored_fasta))
        with stored_manifest.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            missing_fields = set(CONTEXT_FIELDS) - set(reader.fieldnames or [])
            if missing_fields:
                raise ValueError(
                    "Invalid MLVA context manifest; missing columns: "
                    + ", ".join(sorted(missing_fields))
                )
            rows = list(reader)
        loci_by_id = {locus.locus_id: locus for locus in loci}
        seen: set[str] = set()
        for row in rows:
            context_id = row["context_id"]
            locus_id = row["locus_id"]
            if row["schema_version"] != CONTEXT_SCHEMA_VERSION:
                raise ValueError(
                    f"Unsupported MLVA context schema {row['schema_version']!r}; "
                    "rebuild the reference database"
                )
            if context_id in seen:
                raise ValueError(f"Duplicate MLVA context id {context_id!r}")
            seen.add(context_id)
            if context_id not in sequences:
                raise ValueError(f"MLVA context sequence is missing for {context_id!r}")
            if locus_id not in loci_by_id:
                raise ValueError(f"MLVA context refers to unknown locus {locus_id!r}")
            sequence = sequences[context_id]
            if hashlib.sha256(sequence.encode()).hexdigest() != row["sequence_sha256"]:
                raise ValueError(f"MLVA context hash mismatch for {context_id!r}")
            start = int(row.get("repeat_start") or 0)
            end = int(row.get("repeat_end") or 0)
            if not 0 <= start < end <= len(sequence):
                raise ValueError(f"Invalid repeat boundaries for MLVA context {context_id!r}")
        extra_sequences = set(sequences) - seen
        if extra_sequences:
            raise ValueError(
                "MLVA context FASTA has records absent from the manifest: "
                + ", ".join(sorted(extra_sequences))
            )
        return [
            LocusContext(
                context_id=row["context_id"], locus_id=row["locus_id"],
                sequence=sequences[row["context_id"]], reference_id=row["reference_id"],
                taxon_id=row.get("taxon_id", ""), taxon_name=row.get("taxon_name", ""),
                expected_product_size_bp=int(row.get("expected_product_size_bp") or 0),
                repeat_motif=row.get("repeat_motif", ""),
                repeat_unit_length_bp=int(row.get("repeat_unit_length_bp") or 0),
                expected_repeat_count=(
                    float(row["expected_repeat_count"])
                    if row.get("expected_repeat_count")
                    else assembly_equivalent_product_allele(
                        loci_by_id[row["locus_id"]], len(sequences[row["context_id"]])
                    )[1]
                ),
                repeat_start=int(row.get("repeat_start") or 0),
                repeat_end=int(row.get("repeat_end") or 0),
                reference_contig=row.get("reference_contig", ""),
                reference_start=int(row["reference_start"]) if row.get("reference_start") else None,
                reference_end=int(row["reference_end"]) if row.get("reference_end") else None,
                strand=row.get("strand", "+"),
                upstream_flank_bp=int(row.get("upstream_flank_bp") or 0),
                downstream_flank_bp=int(row.get("downstream_flank_bp") or 0),
                source=row.get("source", "database_context"),
            ) for row in rows
        ]

    contexts: list[LocusContext] = []
    for locus in loci:
        records: list[tuple[str, str]] = []
        if not records:
            sequence = _synthetic_product(locus)
            if sequence:
                records = [("panel", sequence)]
        # Identical sequences are represented once per locus while provenance is
        # retained as a semicolon-separated reference set.
        by_sequence: dict[str, list[str]] = defaultdict(list)
        for reference_id, sequence in records:
            by_sequence[sequence.upper()].append(reference_id.split()[0])
        for index, (sequence, reference_ids) in enumerate(sorted(by_sequence.items()), 1):
            repeat_start, repeat_end = _repeat_interval(sequence, locus)
            expected = assembly_equivalent_product_allele(locus, len(sequence))[1]
            references = sorted(set(reference_ids))
            meta_rows = [metadata.get(reference, {}) for reference in references]
            contexts.append(LocusContext(
                context_id=f"ctx{len(contexts)+1:06d}", locus_id=locus.locus_id,
                sequence=sequence, reference_id=";".join(references),
                taxon_id=";".join(sorted({str(row.get("taxon_id") or row.get("taxid") or "") for row in meta_rows} - {""})),
                taxon_name=";".join(sorted({str(row.get("taxon_name") or row.get("organism_name") or "") for row in meta_rows} - {""})),
                expected_product_size_bp=len(sequence), repeat_motif=locus.repeat_motif,
                repeat_unit_length_bp=repeat_unit_length(locus),
                expected_repeat_count=expected, repeat_start=repeat_start,
                repeat_end=repeat_end, upstream_flank_bp=repeat_start,
                downstream_flank_bp=len(sequence) - repeat_end,
                source="database_amplicon" if database else "panel_synthetic",
            ))
    if not contexts:
        raise ValueError(
            "No MLVA locus contexts can be constructed. Rebuild the reference "
            "database or provide a rich panel with flanks and a concrete repeat motif."
        )
    return contexts


def candidate_repeat_counts(locus: Locus, reference_counts: Iterable[int | float | None],
                            maximum: int = 100) -> list[int]:
    if maximum < 1:
        raise ValueError("maximum candidate repeat count must be at least 1")
    lower = max(0, locus.expected_min_repeats)
    upper = min(maximum, locus.expected_max_repeats)
    values = set(range(lower, upper + 1))
    for value in reference_counts:
        if value is not None and float(value).is_integer() and 0 <= int(value) <= maximum:
            values.add(int(value))
    return sorted(values)


def expand_candidate_contexts(contexts: list[LocusContext], loci: list[Locus],
                              maximum: int = 100) -> list[LocusContext]:
    loci_by_id = {locus.locus_id: locus for locus in loci}
    reference_counts: dict[str, list[int | float | None]] = defaultdict(list)
    for context in contexts:
        reference_counts[context.locus_id].append(context.expected_repeat_count)
    expanded: list[LocusContext] = []
    seen: set[tuple[str, str, int]] = set()
    for context in contexts:
        locus = loci_by_id[context.locus_id]
        motif = locus.repeat_motif
        if not motif or set(motif) == {"N"} or context.repeat_end <= context.repeat_start:
            counts = [context.expected_repeat_count] if context.expected_repeat_count is not None else []
        else:
            counts = candidate_repeat_counts(locus, reference_counts[locus.locus_id], maximum)
        for count_value in counts:
            if count_value is None or not float(count_value).is_integer():
                continue
            count = int(count_value)
            sequence = (context.sequence[:context.repeat_start] + motif * count
                        + context.sequence[context.repeat_end:])
            key = (context.locus_id, sequence, count)
            if key in seen:
                continue
            seen.add(key)
            expanded.append(LocusContext(
                **{**asdict(context), "context_id": f"allele{len(expanded)+1:07d}",
                   "sequence": sequence, "expected_repeat_count": count,
                   "repeat_end": context.repeat_start + len(motif) * count,
                   "expected_product_size_bp": len(sequence)}
            ))
    if not expanded:
        raise ValueError("Locus contexts do not encode usable discrete repeat-count candidates")
    return expanded


@lru_cache(maxsize=None)
def check_bowtie2(executable: str = "bowtie2") -> str:
    resolved = shutil.which(executable)
    if resolved is None:
        raise RuntimeError(
            f"Bowtie2 is required for Illumina calling but {executable!r} "
            "is not available on PATH. Install bowtie2 or pass --bowtie2-bin."
        )
    return resolved


@lru_cache(maxsize=None)
def check_bowtie2_build(executable: str = "bowtie2-build") -> str:
    resolved = shutil.which(executable)
    if resolved is None:
        raise RuntimeError(
            f"bowtie2-build is required to index MLVA contexts but {executable!r} "
            "is not available on PATH. Install bowtie2 or pass --bowtie2-build-bin."
        )
    return resolved


def bowtie2_version(executable: str) -> str:
    result = subprocess.run([executable, "--version"], text=True, capture_output=True, check=False)
    return (result.stdout or result.stderr).splitlines()[0].strip() if (result.stdout or result.stderr) else "unknown"


def build_context_index(contexts: list[LocusContext], directory: str | Path,
                        bowtie2_build_bin: str = "bowtie2-build") -> dict[str, Path]:
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    fasta = root / "mlva_context_candidates.fasta"
    manifest = root / "mlva_context_candidates.tsv"
    metadata = root / "mlva_context_index.json"
    prefix = root / "mlva_contexts"
    digest = hashlib.sha256()
    for context in contexts:
        digest.update(json.dumps(context.row(), sort_keys=True).encode())
        digest.update(context.sequence.encode())
    signature = digest.hexdigest()
    small_indexes = [Path(f"{prefix}.{suffix}.bt2") for suffix in (1, 2, 3, 4, "rev.1", "rev.2")]
    large_indexes = [Path(f"{prefix}.{suffix}.bt2l") for suffix in (1, 2, 3, 4, "rev.1", "rev.2")]
    if metadata.exists() and (all(path.exists() for path in small_indexes) or all(path.exists() for path in large_indexes)):
        document = json.loads(metadata.read_text())
        if document.get("context_signature") == signature and document.get("schema_version") == CONTEXT_SCHEMA_VERSION:
            return {"fasta": fasta, "manifest": manifest, "metadata": metadata, "prefix": prefix}
    write_fasta(((context.context_id, context.sequence) for context in contexts), fasta)
    write_tsv((context.row() for context in contexts), manifest, CONTEXT_FIELDS)
    executable = check_bowtie2_build(bowtie2_build_bin)
    result = subprocess.run([executable, "--quiet", str(fasta), str(prefix)], text=True,
                            capture_output=True, check=False)
    if result.returncode:
        raise RuntimeError(f"bowtie2-build failed ({result.returncode}): {(result.stderr or result.stdout).strip()}")
    metadata.write_text(json.dumps({"schema_version": CONTEXT_SCHEMA_VERSION,
                                    "context_signature": signature,
                                    "candidate_contexts": len(contexts)},
                                   indent=2, sort_keys=True) + "\n")
    return {"fasta": fasta, "manifest": manifest, "metadata": metadata, "prefix": prefix}


def bowtie2_mapping_command(executable: str, index_prefix: Path, reads1: Path,
                            reads2: Path | None, output_sam: Path, threads: int,
                            include_secondary: bool = True, orphans: Path | None = None) -> list[str]:
    command = [executable, "--very-sensitive", "--no-unal", "--no-mixed", "-p", str(threads)]
    command += ["-a"] if include_secondary else ["-k", "1"]
    command += ["-x", str(index_prefix)]
    if reads2 is not None:
        command += ["-1", str(reads1), "-2", str(reads2)]
    else:
        command += ["-U", str(reads1)]
    if orphans is not None and orphans.exists() and orphans.stat().st_size:
        command += ["-U", str(orphans)]
    command += ["-S", str(output_sam)]
    return command


def run_bowtie2(command: list[str], log_path: Path) -> None:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    log_path.write_text((result.stderr or "") + (result.stdout or ""))
    if result.returncode:
        raise RuntimeError(f"Bowtie2 MLVA context mapping failed ({result.returncode}); see {log_path}")


def estimate_insert_size_distribution(values: Iterable[int], minimum_pairs: int = 3) -> InsertSizeEstimate:
    sizes = [abs(int(value)) for value in values if int(value)]
    if len(sizes) < minimum_pairs:
        return InsertSizeEstimate(None, None, None, None, len(sizes))
    median = statistics.median(sizes)
    deviations = [abs(value - median) for value in sizes]
    mad = statistics.median(deviations)
    robust_limit = max(6 * mad, 50.0)
    retained = [value for value in sizes if abs(value - median) <= robust_limit]
    return InsertSizeEstimate(
        statistics.median(retained), statistics.mean(retained),
        statistics.stdev(retained) if len(retained) > 1 else 0.0,
        statistics.median(abs(value - statistics.median(retained)) for value in retained),
        len(retained),
    )


def _alignment_score(alignment: pysam.AlignedSegment) -> float:
    return float(alignment.get_tag("AS")) if alignment.has_tag("AS") else float(alignment.query_alignment_length)


def _crosses(alignment: pysam.AlignedSegment, position: int) -> bool:
    return alignment.reference_start < position < (alignment.reference_end or 0)


def _full_spanning(alignment: pysam.AlignedSegment, context: LocusContext) -> bool:
    return alignment.reference_start <= context.repeat_start and (alignment.reference_end or 0) >= context.repeat_end


def _pair_spans(left: pysam.AlignedSegment, right: pysam.AlignedSegment,
                context: LocusContext) -> bool:
    starts = sorted((left.reference_start, right.reference_start))
    ends = sorted(((left.reference_end or 0), (right.reference_end or 0)))
    return starts[0] < context.repeat_start and ends[-1] > context.repeat_end and starts[1] >= context.repeat_end


def _normal_logpdf(value: float, center: float, spread: float) -> float:
    spread = max(spread, 10.0)
    return -0.5 * ((value - center) / spread) ** 2 - math.log(spread)


def _softmax(scores: dict[int | float, float]) -> dict[int | float, float]:
    if not scores:
        return {}
    maximum = max(scores.values())
    weights = {key: math.exp(min(0.0, value - maximum)) for key, value in scores.items()}
    total = sum(weights.values())
    return {key: value / total for key, value in weights.items()}


def infer_mapping_calls(
    bam_path: str | Path, contexts: list[LocusContext], loci: list[Locus], sample_id: str,
    min_mapq: int = 10, min_supporting_fragments: int = 3, min_spanning_pairs: int = 2,
    confidence_threshold: float = 0.8,
) -> tuple[list[dict], InsertSizeEstimate]:
    by_context = {context.context_id: context for context in contexts}
    alignments: dict[tuple[str, str], list[pysam.AlignedSegment]] = defaultdict(list)
    insert_values: dict[str, int] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for alignment in bam.fetch(until_eof=True):
            if alignment.is_unmapped or alignment.mapping_quality < min_mapq:
                continue
            context = by_context.get(bam.get_reference_name(alignment.reference_id))
            if context is None:
                continue
            alignments[(alignment.query_name, context.context_id)].append(alignment)
            # Same-flank pairs estimate the library without VNTR-length dependence.
            if alignment.is_proper_pair and alignment.is_read1 and alignment.template_length:
                mate_start = alignment.next_reference_start
                same_left = max(alignment.reference_start, mate_start) < context.repeat_start
                same_right = min(alignment.reference_start, mate_start) >= context.repeat_end
                if same_left or same_right:
                    insert_values.setdefault(alignment.query_name, abs(alignment.template_length))
    insert = estimate_insert_size_distribution(insert_values.values())
    # Recover overlap-spanning molecules before candidate scoring. This is not
    # assembly: it is an exact two-read molecule observation and uses the same
    # product-length calibration as the assembly caller.
    from .models import ReadPair, ReadRecord
    from .short_reads import _merge_overlap_is_repeat_only, merge_read_pair

    direct_alleles: dict[tuple[str, int | float], set[str]] = defaultdict(set)
    measured_molecule_loci: set[tuple[str, str]] = set()
    loci_by_id = {locus.locus_id: locus for locus in loci}
    for (molecule, context_id), records in alignments.items():
        context = by_context[context_id]
        molecule_locus = molecule, context.locus_id
        if molecule_locus in measured_molecule_loci:
            continue
        measured_molecule_loci.add(molecule_locus)
        first = next((record for record in records if record.is_read1), None)
        second = next((record for record in records if record.is_read2), None)
        if first is None and records and not any(record.is_paired for record in records):
            first = records[0]
        if first is None:
            continue
        pair = ReadPair(
            molecule,
            ReadRecord(f"{molecule}/1", first.get_forward_sequence()),
            ReadRecord(f"{molecule}/2", second.get_forward_sequence()) if second is not None else None,
        )
        merged = merge_read_pair(pair)
        if merged is not None and _merge_overlap_is_repeat_only(
            pair, merged, loci_by_id[context.locus_id]
        ):
            merged = None
        sequences = [merged.sequence] if merged is not None else [first.query_sequence]
        for sequence in sequences:
            measurements = [
                measure_locus_product(sequence, loci_by_id[context.locus_id], source="paired_read_molecule"),
                measure_locus_product(revcomp(sequence), loci_by_id[context.locus_id], source="paired_read_molecule"),
            ]
            measured = max(measurements, key=lambda value: (value.called_allele is not None, value.confidence or 0))
            if measured.called_allele is not None and measured.status == "FULL_PRODUCT":
                direct_alleles[(context.locus_id, measured.called_allele)].add(molecule)
                break
    evidence: dict[tuple[str, int | float], _CandidateEvidence] = {}
    any_molecules: dict[str, set[str]] = defaultdict(set)
    locus_mapqs: dict[str, list[int]] = defaultdict(list)
    for (molecule, context_id), records in alignments.items():
        context = by_context[context_id]
        if context.expected_repeat_count is None:
            continue
        key = (context.locus_id, context.expected_repeat_count)
        item = evidence.setdefault(key, _CandidateEvidence(
            context.locus_id, context.expected_repeat_count, set(), set(), set(), set(),
            set(), set(), [], [], []))
        item.context_ids.add(context_id)
        item.molecules.add(molecule)
        any_molecules[context.locus_id].add(molecule)
        for record in records:
            item.mapqs.append(record.mapping_quality)
            locus_mapqs[context.locus_id].append(record.mapping_quality)
            item.alignment_scores.append(_alignment_score(record))
            if _crosses(record, context.repeat_start) or _crosses(record, context.repeat_end):
                item.junction_reads.add(f"{molecule}/{int(record.is_read2)+1}")
            if _full_spanning(record, context):
                item.full_reads.add(f"{molecule}/{int(record.is_read2)+1}")
            if any(operation in (1, 2) for operation, _length in (record.cigartuples or [])):
                item.indel_reads.add(f"{molecule}/{int(record.is_read2)+1}")
        first_records = [record for record in records if record.is_read1]
        second_records = [record for record in records if record.is_read2]
        if first_records and second_records:
            first, second = first_records[0], second_records[0]
            if _pair_spans(first, second, context):
                item.spanning_pairs.add(molecule)
                tlen = abs(first.template_length or second.template_length)
                if tlen and insert.median is not None:
                    spread = max(insert.standard_deviation or 0, (insert.mad or 0) * 1.4826, 10)
                    item.geometry_log_likelihoods.append(_normal_logpdf(tlen, insert.median, spread))

    rows: list[dict] = []
    for locus_id in [locus.locus_id for locus in loci]:
        candidates = {count: item for (candidate_locus, count), item in evidence.items()
                      if candidate_locus == locus_id}
        detected = any_molecules[locus_id]
        if not detected:
            state, reason = "no_evidence", "no reads passed locus-context mapping thresholds"
            probabilities: dict[int | float, float] = {}
        else:
            scores = {}
            for count, item in candidates.items():
                independent = len(item.molecules)
                favored_mapq = statistics.mean(item.mapqs) / 60 if item.mapqs else 0
                direct = 4 * len(item.full_reads) + 2.5 * len(item.junction_reads)
                pair = 3 * len(item.spanning_pairs)
                indel = 1.5 * len(item.indel_reads)
                alignment = statistics.mean(item.alignment_scores) / 100 if item.alignment_scores else 0
                geometry = sum(item.geometry_log_likelihoods) / max(len(item.geometry_log_likelihoods), 1)
                molecule_calls = len(direct_alleles.get((locus_id, count), set()))
                scores[count] = (math.log1p(independent) + favored_mapq + direct + pair
                                 + indel + alignment + geometry + 20 * molecule_calls)
            probabilities = _softmax(scores)
            ranked = sorted(probabilities, key=lambda value: (-probabilities[value], float(value)))
            best = ranked[0] if ranked else None
            best_item = candidates.get(best) if best is not None else None
            runner = probabilities[ranked[1]] if len(ranked) > 1 else 0.0
            directly_observed = {
                count: molecules
                for (direct_locus, count), molecules in direct_alleles.items()
                if direct_locus == locus_id and molecules
            }
            has_direct = bool(best is not None and directly_observed.get(best))
            if best is None:
                state, reason = "detected_unresolved", "mapped reads do not support a countable candidate"
            elif len(directly_observed) > 1:
                state, reason = "ambiguous", "multiple alleles are directly supported by independent molecules"
            elif len(detected) < min_supporting_fragments and has_direct:
                state, reason = "low_coverage", "insufficient independent supporting fragments"
            elif len(detected) < min_supporting_fragments:
                state, reason = "detected_unresolved", "insufficient independent supporting fragments"
            elif not (
                direct_alleles.get((locus_id, best))
                or (
                    insert.median is not None
                    and len(best_item.spanning_pairs) >= min_spanning_pairs
                )
            ):
                state, reason = "detected_unresolved", "flank evidence establishes presence but does not span the VNTR"
            elif probabilities[best] < confidence_threshold or probabilities[best] - runner < 0.2:
                state, reason = "ambiguous", "neighboring repeat counts are not decisively separated"
            else:
                state, reason = "called", "single discrete repeat count is supported"
        ranked = sorted(probabilities, key=lambda value: (-probabilities[value], float(value)))
        best = ranked[0] if ranked else None
        best_item = candidates.get(best) if best is not None else None
        direct_ranked = sorted(
            ((count, len(molecules)) for (direct_locus, count), molecules in direct_alleles.items()
             if direct_locus == locus_id and molecules),
            key=lambda item: (-item[1], float(item[0])),
        )
        dominant_direct = (
            direct_ranked[0][0]
            if direct_ranked and (len(direct_ranked) == 1 or direct_ranked[0][1] > direct_ranked[1][1])
            else None
        )
        called = best if state == "called" or (
            state == "low_coverage" and direct_alleles.get((locus_id, best))
        ) else dominant_direct if state == "ambiguous" else ""
        unit = repeat_unit_length(loci_by_id[locus_id])
        nonrepeat = None
        if best_item:
            representative = by_context[sorted(best_item.context_ids)[0]]
            nonrepeat = len(representative.sequence) - unit * float(representative.expected_repeat_count)
        rows.append({
            "sample_id": sample_id, "locus_id": locus_id, "state": state,
            "repeat_count": called,
            "locus_length_bp": int(nonrepeat + unit * float(best)) if nonrepeat is not None and best is not None else "",
            "confidence": round(probabilities.get(best, 0), 6) if best is not None else 0,
            "supporting_fragments": len(detected),
            "proper_spanning_pairs": len(best_item.spanning_pairs) if best_item else 0,
            "junction_reads": len(best_item.junction_reads) if best_item else 0,
            "full_spanning_reads": len(best_item.full_reads) if best_item else 0,
            "cigar_indel_reads": len(best_item.indel_reads) if best_item else 0,
            "mean_mapq": round(statistics.mean(locus_mapqs[locus_id]), 3) if locus_mapqs[locus_id] else "",
            "median_mapq": statistics.median(locus_mapqs[locus_id]) if locus_mapqs[locus_id] else "",
            "candidate_scores": ";".join(f"{count}:{probabilities[count]:.6f}" for count in ranked),
            "best_contexts": ";".join(sorted(best_item.context_ids)) if best_item else "",
            "context_taxa": ";".join(sorted({by_context[c].taxon_name or by_context[c].taxon_id for c in best_item.context_ids} - {""})) if best_item else "",
            "reason": reason,
        })
    return rows, insert


def compatibility_call_rows(evidence_rows: list[dict], contexts: list[LocusContext],
                            loci: list[Locus] | None = None) -> list[dict]:
    contexts_by_id = {context.context_id: context for context in contexts}
    loci_by_id = {locus.locus_id: locus for locus in (loci or [])}
    output = []
    for row in evidence_rows:
        distribution = row["candidate_scores"]
        ranked = [entry.rsplit(":", 1) for entry in distribution.split(";") if entry]
        second = ranked[1] if len(ranked) > 1 else ("", "")
        context_ids = [value for value in row["best_contexts"].split(";") if value]
        references = sorted({reference for context_id in context_ids
                             for reference in contexts_by_id[context_id].reference_id.split(";") if reference})
        status_map = {"called": "PASS", "detected_unresolved": "PRESENT_COUNT_UNKNOWN",
                      "low_coverage": "LOW_DEPTH", "ambiguous": "AMBIGUOUS",
                      "no_evidence": "NOT_FOUND", "mapping_conflict": "MAPPING_CONFLICT"}
        output.append({
            "sample_id": row["sample_id"], "locus_id": row["locus_id"],
            "present": "no" if row["state"] == "no_evidence" else "yes",
            "repeat_count": row["repeat_count"], "repeat_count_raw": row["repeat_count"],
            "product_size_bp": row["locus_length_bp"], "read_depth": row["supporting_fragments"],
            "primary_read_depth": row["supporting_fragments"], "mean_coverage": "",
            "allele_confidence": row["confidence"], "second_best_repeat_count": second[0],
            "second_best_probability": second[1], "inference_method": "bowtie2_context_likelihood",
            "dominant_variant_fraction": row["confidence"], "num_candidate_variants": len(ranked),
            "num_confirmed_secondary_variants": int(row["state"] == "ambiguous"),
            "secondary_alleles": ";".join(":".join(item) for item in ranked[1:]),
            "allele_distribution": distribution, "status": status_map[row["state"]],
            "evidence": row["reason"], "read_technology": "illumina",
            "evidence_class": (
                "MULTIPLE_ALLELES" if row["state"] == "ambiguous" and row["repeat_count"] != ""
                else "BOUNDARY_SPANNING_READ_PAIR" if row["state"] == "called"
                or (row["state"] == "detected_unresolved" and row["supporting_fragments"]
                    and row["junction_reads"] >= 2)
                else "PARTIAL_REPEAT_EVIDENCE" if row["state"] == "detected_unresolved"
                else "LOW_DEPTH" if row["state"] == "low_coverage"
                else row["state"].upper()
            ),
            "informative_molecule_count": row["supporting_fragments"],
            "boundary_1_support": row["junction_reads"], "boundary_2_support": row["junction_reads"],
            "both_boundary_support": row["proper_spanning_pairs"],
            "repeat_count_min": (
                loci_by_id[row["locus_id"]].expected_min_repeats
                if row["locus_id"] in loci_by_id else min(int(float(item[0])) for item in ranked)
                if row["state"] == "detected_unresolved" and ranked else ""
            ),
            "repeat_count_max": (
                loci_by_id[row["locus_id"]].expected_max_repeats
                if row["locus_id"] in loci_by_id else max(int(float(item[0])) for item in ranked)
                if row["state"] == "detected_unresolved" and ranked else ""
            ),
            "repeat_count_interval_reason": row["reason"], "confidence_reason": row["reason"],
            "short_read_warning": "" if row["state"] == "called" else row["reason"],
            "recruited_read_pairs": row["supporting_fragments"],
            "failure_reason": "" if row["state"] == "called" else row["reason"],
            "primary_allele_support": row["supporting_fragments"], "secondary_allele_support": 0,
            "informative_allele_reads": row["full_spanning_reads"] + row["junction_reads"],
            "uninformative_locus_reads": max(0, row["supporting_fragments"] - row["junction_reads"]),
            "estimated_primary_fraction": row["confidence"], "estimated_secondary_fraction": "",
            "mixture_status": "AMBIGUOUS" if row["state"] == "ambiguous" else "SINGLE",
            "evidence_sources": "bowtie2_mapping", "mapping_state": row["state"],
            "supporting_fragments": row["supporting_fragments"],
            "proper_spanning_pairs": row["proper_spanning_pairs"], "junction_read_count": row["junction_reads"],
            "full_spanning_read_count": row["full_spanning_reads"], "cigar_indel_read_count": row["cigar_indel_reads"],
            "mean_mapq": row["mean_mapq"], "median_mapq": row["median_mapq"],
            "candidate_allele_scores": distribution, "best_reference_contexts": row["best_contexts"],
            "reference_context_taxa": row["context_taxa"], "reference_context_provenance": ";".join(references),
            "mlva_method": "Bowtie2 short-read mapping", "query_sequence": "",
        })
    return output


def write_run_metadata(path: Path, bowtie2: str, database_path: str | None,
                       insert: InsertSizeEstimate, parameters: dict) -> None:
    path.write_text(json.dumps({
        "schema_version": "1.0", "mlvamaps_version": __version__,
        "method": "bowtie2_context_mapping",
        "bowtie2_version": bowtie2, "database": database_path or "panel-derived",
        "insert_size": asdict(insert), "parameters": parameters,
    }, indent=2, sort_keys=True) + "\n")


def run_mapping_short_read_call(
    *, reads1_path: str, reads2_path: str | None, loci_path: str | None,
    primers_path: str | None, profiles_path: str | None, database_path: str | None,
    outdir: str, sample_id: str, sample_metadata: dict[str, str] | None,
    short_min_read_length: int, short_min_mean_quality: float, short_trim_quality: int,
    short_min_pair_retention: float, min_depth: int, threads: int,
    keep_intermediates: bool, sample_mode: str, bowtie2_bin: str,
    bowtie2_build_bin: str, short_min_mapping_quality: int,
    short_min_spanning_pairs: int, short_confidence_threshold: float,
    short_max_candidate_repeat_count: int, short_consider_secondary: bool,
    show_progress: bool = True,
) -> dict[str, Path]:
    """Run one competitive map and emit the established MLVAmap result family."""
    # Lazy imports avoid a module cycle with the shared Illumina helpers.
    from .concurrency import resolve_threads
    from .io import read_fastq_pairs
    from .primers import read_loci_or_primers
    from .report import write_report
    from .sample_metadata import MYOGA_SAMPLE_FIELDS, myoga_sample_row, write_csv
    from .short_reads import (
        SAMPLE_SUMMARY_FIELDS, SHORT_CALL_FIELDS, SHORT_CALL_EXTRA_FIELDS,
        SHORT_QC_FIELDS, _allele_rows, qc_read_pairs,
    )
    from .pipeline import ALLELE_DISTRIBUTION_FIELDS, MATCH_FIELDS, REPEAT_COUNT_FIELDS, allele_distribution_rows

    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    loci = read_loci_or_primers(loci_path, primers_path)
    profiles = read_profiles(profiles_path)
    bowtie2 = check_bowtie2(bowtie2_bin)
    check_bowtie2_build(bowtie2_build_bin)
    thread_count = resolve_threads(threads)
    filtered1 = output / "filtered_reads_1.fastq.gz"
    filtered2 = output / "filtered_reads_2.fastq.gz"
    orphans = output / "filtered_orphan_reads.fastq.gz"
    counters: dict[str, int] = defaultdict(int)

    if show_progress:
        print(f"[{sample_id}] Filtering and validating Illumina read pairs")
    with open_text(filtered1, "wt") as first_handle, open_text(filtered2, "wt") as second_handle, open_text(orphans, "wt") as orphan_handle:
        chunk = []

        def consume() -> None:
            retained, metrics = qc_read_pairs(
                chunk, short_min_read_length, short_min_mean_quality,
                short_trim_quality, short_min_pair_retention,
            )
            for key, value in metrics.items():
                counters[key] += int(value)
            for pair in retained:
                if pair.read2 is None:
                    target = orphan_handle if reads2_path else first_handle
                    quality = pair.read1.quality or "I" * len(pair.read1.sequence)
                    target.write(f"@{pair.read1.read_id}\n{pair.read1.sequence}\n+\n{quality}\n")
                else:
                    for record, handle in ((pair.read1, first_handle), (pair.read2, second_handle)):
                        quality = record.quality or "I" * len(record.sequence)
                        handle.write(f"@{record.read_id}\n{record.sequence}\n+\n{quality}\n")
            chunk.clear()

        for pair in read_fastq_pairs(reads1_path, reads2_path):
            chunk.append(pair)
            if len(chunk) == 5000:
                consume()
        if chunk:
            consume()

    contexts = expand_candidate_contexts(
        load_locus_contexts(loci, database_path), loci, short_max_candidate_repeat_count
    )
    work = output / "bowtie2_mlva"
    if show_progress:
        print(f"[{sample_id}] Building/using Bowtie2 MLVA context index")
    index = build_context_index(contexts, work / "index", bowtie2_build_bin)
    sam = work / "locus_context_alignments.sam"
    bam = work / "locus_context_alignments.bam"
    if show_progress:
        print(f"[{sample_id}] Mapping paired-end reads with Bowtie2")
    command = bowtie2_mapping_command(
        bowtie2, index["prefix"], filtered1,
        filtered2 if reads2_path else None, sam, thread_count,
        include_secondary=short_consider_secondary,
        orphans=orphans if reads2_path else None,
    )
    run_bowtie2(command, work / "bowtie2.log")
    with pysam.AlignmentFile(str(sam), "r") as source, pysam.AlignmentFile(str(bam), "wb", template=source) as destination:
        for alignment in source.fetch(until_eof=True):
            destination.write(alignment)
    evidence_rows, insert = infer_mapping_calls(
        bam, contexts, loci, sample_id, short_min_mapping_quality, min_depth,
        short_min_spanning_pairs, short_confidence_threshold,
    )
    write_tsv(evidence_rows, output / "short_read_mapping_evidence.tsv", MAPPING_EVIDENCE_FIELDS)
    calls = compatibility_call_rows(evidence_rows, contexts, loci)
    mapping_fields = [
        "mapping_state", "supporting_fragments", "proper_spanning_pairs",
        "junction_read_count", "full_spanning_read_count", "cigar_indel_read_count",
        "mean_mapq", "median_mapq", "candidate_allele_scores",
        "best_reference_contexts", "reference_context_taxa",
        "reference_context_provenance", "mlva_method",
    ]
    call_fields = SHORT_CALL_FIELDS + [field for field in mapping_fields if field not in SHORT_CALL_FIELDS]
    calls_path = output / "calls.tsv"
    write_tsv(calls, calls_path, call_fields)
    write_tsv(calls, output / "locus_repeat_counts.tsv",
              REPEAT_COUNT_FIELDS + [field for field in SHORT_CALL_EXTRA_FIELDS + mapping_fields if field not in REPEAT_COUNT_FIELDS])
    allele_rows = _allele_rows(calls)
    write_tsv(allele_distribution_rows(sample_id, calls), output / "allele_probability_distribution.tsv", ALLELE_DISTRIBUTION_FIELDS)
    fingerprint, probabilistic = build_fingerprint(sample_id, allele_rows, loci)
    write_tsv(fingerprint, output / "mlva_fingerprint.tsv", ["sample_id"] + [locus.locus_id for locus in loci])
    write_tsv(probabilistic, output / "mlva_fingerprint_probabilistic.tsv",
              ["sample_id", "locus_id", "repeat_count", "posterior_probability"])
    matches = match_profiles(sample_id, fingerprint[0], profiles, allele_rows=allele_rows)
    write_tsv(matches, output / "profile_matches.tsv", MATCH_FIELDS)
    write_tsv(profile_match_locus_rows(sample_id, fingerprint[0], profiles, matches, allele_rows),
              output / "profile_match_loci.tsv", PROFILE_MATCH_LOCUS_FIELDS)

    qc_values = dict(counters)
    qc_values.update({
        "insert_size_median": "" if insert.median is None else round(insert.median, 3),
        "insert_size_mean": "" if insert.mean is None else round(insert.mean, 3),
        "insert_size_standard_deviation": "" if insert.standard_deviation is None else round(insert.standard_deviation, 3),
        "insert_size_mad": "" if insert.mad is None else round(insert.mad, 3),
        "insert_size_pairs": insert.pairs_used,
    })
    write_tsv(({"sample_id": sample_id, "metric": key, "value": value}
               for key, value in sorted(qc_values.items())),
              output / "short_read_qc_summary.tsv", SHORT_QC_FIELDS)
    recruitment = [{
        "sample_id": sample_id, "locus_id": row["locus_id"],
        "read_pairs_examined": counters.get("input_pairs", 0),
        "read_pairs_recruited": row["supporting_fragments"],
        "uniquely_recruited_pairs": row["supporting_fragments"], "ambiguous_pairs": 0,
        "discordant_pairs": int(row["state"] == "mapping_conflict"), "orphan_reads": 0,
        "mean_mapping_quality": row["mean_mapq"], "mean_alignment_identity": "",
        "presence_status": "NO_EVIDENCE" if row["state"] == "no_evidence" else
                           "PRESENT_GENOTYPED" if row["state"] == "called" else "PRESENT_UNTYPED",
        "mapped_reads": row["supporting_fragments"], "full_product_reads": row["full_spanning_reads"],
        "genotype_informative_reads": row["junction_reads"], "candidate_alleles": row["candidate_scores"],
        "reference_source": "bowtie2_competitive_contexts",
    } for row in evidence_rows]
    recruitment_fields = [
        "sample_id", "locus_id", "read_pairs_examined", "read_pairs_recruited",
        "uniquely_recruited_pairs", "ambiguous_pairs", "discordant_pairs", "orphan_reads",
        "mean_mapping_quality", "mean_alignment_identity", "presence_status", "mapped_reads",
        "full_product_reads", "genotype_informative_reads", "candidate_alleles", "reference_source",
    ]
    write_tsv(recruitment, output / "short_read_recruitment_summary.tsv", recruitment_fields)
    write_tsv(recruitment, output / "locus_presence.tsv", recruitment_fields)
    best = matches[0] if matches else {}
    states = defaultdict(int)
    for row in evidence_rows:
        states[row["state"]] += 1
    summary = {
        "sample_id": sample_id, "input_read_1": str(Path(reads1_path)),
        "input_read_2": "" if reads2_path is None else str(Path(reads2_path)),
        "read_technology": "illumina", "sample_mode": sample_mode,
        "total_reads": counters.get("input_reads", 0), "total_read_pairs": counters.get("input_pairs", 0),
        "retained_reads": counters.get("retained_reads", 0), "retained_pairs": counters.get("retained_pairs", 0),
        "callable_loci": states["called"], "complete_loci": states["called"],
        "partial_loci": states["detected_unresolved"] + states["low_coverage"] + states["ambiguous"],
        "presence_only_loci": states["detected_unresolved"], "mixed_loci": states["ambiguous"],
        "missing_loci": states["no_evidence"], "best_profile_id": best.get("best_profile_id", ""),
        "best_profile_distance": best.get("distance", ""), "profile_confidence": best.get("confidence", ""),
        "run_status": "success_partial" if states["called"] < len(loci) else "success",
        "warnings": ";".join(sorted({row["reason"] for row in evidence_rows if row["state"] != "called"})),
    }
    write_tsv([summary], output / "sample_summary.tsv", SAMPLE_SUMMARY_FIELDS)
    write_csv([myoga_sample_row(sample_id, sample_metadata, summary, fingerprint[0], len(loci))],
              output / "myoga_samples.csv", MYOGA_SAMPLE_FIELDS)
    write_csv([{"genome_id": sample_id, "sample_id": sample_id, "locus_id": row["locus_id"],
                "repeat_count": row["repeat_count"], "repeat_count_min": "", "repeat_count_max": "",
                "evidence_class": row["mapping_state"], "confidence": row["allele_confidence"]}
               for row in calls], output / "myoga_loci.csv",
              ["genome_id", "sample_id", "locus_id", "repeat_count", "repeat_count_min",
               "repeat_count_max", "evidence_class", "confidence"])
    write_run_metadata(output / "short_read_run_metadata.json", bowtie2_version(bowtie2),
                       database_path, insert, {"minimum_mapq": short_min_mapping_quality,
                       "minimum_supporting_fragments": min_depth,
                       "minimum_spanning_pairs": short_min_spanning_pairs,
                       "confidence_threshold": short_confidence_threshold,
                       "maximum_candidate_repeat_count": short_max_candidate_repeat_count,
                       "secondary_alignments": short_consider_secondary})
    write_report(output, sample_id, allele_rows, loci, matches, profiles,
                 presence_rows=recruitment, local_assembly_rows=[], short_read_rows=calls)
    report_path = output / "report.html"
    report_path.write_text(report_path.read_text().replace(
        "MLVA analysis report", "MLVA analysis report · Method: Bowtie2 short-read mapping", 1
    ))
    if show_progress:
        detected = len(loci) - states["no_evidence"]
        print(f"[{sample_id}] {detected}/{len(loci)} loci have mapping evidence")
        print(f"[{sample_id}] {states['called']} called, {states['ambiguous']} ambiguous, "
              f"{states['low_coverage']} low coverage, {states['no_evidence']} no evidence")
    if not keep_intermediates:
        sam.unlink(missing_ok=True)
        for path in (filtered1, filtered2, orphans):
            path.unlink(missing_ok=True)
        shutil.rmtree(work, ignore_errors=True)
    return {
        "outdir": output, "calls": calls_path, "repeat_counts": output / "locus_repeat_counts.tsv",
        "allele_distribution": output / "allele_probability_distribution.tsv",
        "fingerprint": output / "mlva_fingerprint.tsv", "profile_matches": output / "profile_matches.tsv",
        "profile_match_loci": output / "profile_match_loci.tsv", "report": report_path,
        "sample_summary": output / "sample_summary.tsv", "myoga_samples": output / "myoga_samples.csv",
        "myoga_loci": output / "myoga_loci.csv", "short_read_qc": output / "short_read_qc_summary.tsv",
        "short_read_recruitment": output / "short_read_recruitment_summary.tsv",
        "short_read_mapping": output / "short_read_mapping_evidence.tsv",
        "run_metadata": output / "short_read_run_metadata.json",
        **({"mapping_bam": bam} if keep_intermediates else {}),
    }