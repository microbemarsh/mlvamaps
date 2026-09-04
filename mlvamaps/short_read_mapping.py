"""Illumina orchestration and legacy benchmark mapping helpers.

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
from collections import defaultdict
from dataclasses import asdict, dataclass
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
    sequence_reference_match_rows,
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


def run_mapping_short_read_call(
    *, reads1_path: str, reads2_path: str | None, loci_path: str | None,
    primers_path: str | None, profiles_path: str | None, database_path: str | None,
    outdir: str, sample_id: str, sample_metadata: dict[str, str] | None,
    short_min_read_length: int, short_min_mean_quality: float, short_trim_quality: int,
    short_min_pair_retention: float, min_depth: int, threads: int,
    keep_intermediates: bool, sample_mode: str, minimap2_bin: str,
    short_min_mapping_quality: int,
    short_min_spanning_pairs: int, short_confidence_threshold: float,
    short_max_candidate_repeat_count: int, short_consider_secondary: bool,
    mafft_bin: str, raxml_ng_bin: str, epa_ng_bin: str, raxml_model: str,
    phylogeny_snp_weight: float, phylogeny_repeat_weight: float,
    reference_metadata_path: str | None, target_taxon_id: str | None,
    taxon_calibration_path: str | None, taxon_alpha: float | None,
    taxon_min_loci: int | None, taxon_min_locus_fraction: float,
    taxon_bootstrap_replicates: int, taxon_min_bootstrap_support: float,
    taxon_max_mean_placement_entropy: float | None,
    taxon_min_median_placement_lwr: float | None,
    taxon_identification: bool | None, taxon_k: int, taxon_minimum_margin: float,
    show_progress: bool = True,
) -> dict[str, Path]:
    """Run competitive minimap2 mapping and emit established output views."""
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

    from .unified_fastq import common_calls_to_compatibility, run_unified_fastq_inference

    if show_progress:
        print(f"[{sample_id}] Competitively mapping Illumina molecules with minimap2")
    common_calls, molecule_evidence, _molecule_calls, unified_paths = run_unified_fastq_inference(
        reads1=filtered1,
        reads2=filtered2 if reads2_path else None,
        loci=loci,
        database_path=database_path,
        outdir=output,
        sample_id=sample_id,
        technology="illumina",
        minimap2_bin=minimap2_bin,
        threads=thread_count,
        minimum_molecules=min_depth,
        minimum_probability=short_confidence_threshold,
        maximum_candidate_repeat_count=short_max_candidate_repeat_count,
        keep_alignments=keep_intermediates,
    )
    calls = common_calls_to_compatibility(common_calls)
    by_locus = {str(row["locus"]): row for row in common_calls}
    evidence_rows = []
    for locus in loci:
        row = by_locus[locus.locus_id]
        state = str(row["status"])
        evidence_rows.append({
            "sample_id": sample_id, "locus_id": locus.locus_id,
            "state": "no_evidence" if state == "not_found" else state,
            "repeat_count": row["repeat_count"], "locus_length_bp": "",
            "confidence": row["best_probability"],
            "supporting_fragments": row["molecule_support"],
            "proper_spanning_pairs": 0, "junction_reads": row["junction_support"],
            "full_spanning_reads": row["full_span_support"], "cigar_indel_reads": 0,
            "mean_mapq": "", "median_mapq": "",
            "candidate_scores": row["candidate_distribution"], "best_contexts": "",
            "context_taxa": "", "reason": calls[len(evidence_rows)]["evidence"],
        })
    write_tsv(evidence_rows, output / "short_read_mapping_evidence.tsv", MAPPING_EVIDENCE_FIELDS)
    insert = InsertSizeEstimate(None, None, None, None, 0)
    for call, common in zip(calls, common_calls):
        call.update({
            "read_technology": "illumina",
            "evidence_class": str(common["status"]).upper(),
            "informative_molecule_count": common["molecule_support"],
            "boundary_1_support": common["junction_support"],
            "boundary_2_support": common["junction_support"],
            "both_boundary_support": common["full_span_support"],
            "repeat_count_min": "", "repeat_count_max": "",
            "repeat_count_interval_reason": "", "confidence_reason": call["evidence"],
            "short_read_warning": "" if common["status"] == "called" else call["evidence"],
            "recruited_read_pairs": common["molecule_support"],
            "failure_reason": "" if common["status"] in {"called", "low_coverage"} else call["evidence"],
            "primary_allele_support": common["molecule_support"], "secondary_allele_support": "",
            "informative_allele_reads": common["full_span_support"],
            "uninformative_locus_reads": max(0, int(common["molecule_support"]) - int(common["full_span_support"])),
            "estimated_primary_fraction": common["dominant_fraction"],
            "estimated_secondary_fraction": common["secondary_fraction"],
            "mixture_status": "MIXED" if common["status"] == "mixed" else "SINGLE",
            "evidence_sources": "competitive_minimap2", "mapping_state": common["status"],
            "supporting_fragments": common["molecule_support"], "proper_spanning_pairs": 0,
            "junction_read_count": common["junction_support"],
            "full_spanning_read_count": common["full_span_support"], "cigar_indel_read_count": 0,
            "mean_mapq": "", "median_mapq": "", "candidate_allele_scores": common["candidate_distribution"],
            "best_reference_contexts": "", "reference_context_taxa": "",
            "reference_context_provenance": "", "mlva_method": "competitive minimap2 shared inference",
            "query_sequence": "", "direct_product_support": common["direct_product_support"],
        })
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
    write_tsv(profile_match_locus_rows(sample_id, fingerprint[0], profiles, matches, allele_rows),
              output / "profile_match_loci.tsv", PROFILE_MATCH_LOCUS_FIELDS)

    phylogeny_paths: dict[str, Path] = {}
    if database_path:
        from .phylogeny import run_phylogenetic_placement

        query_sequences = dict(read_fasta(unified_paths["taxonomic_query_sequences"]))
        if show_progress:
            print(
                f"[{sample_id}] Matching {len(query_sequences):,} Illumina marker "
                "sequence(s) against the reference database"
            )
        phylogeny_paths = run_phylogenetic_placement(
            query_sequences, database_path, output, sample_id, loci, thread_count,
            mafft_bin=mafft_bin, raxml_ng_bin=raxml_ng_bin, epa_ng_bin=epa_ng_bin,
            raxml_model=raxml_model, snp_weight=phylogeny_snp_weight,
            repeat_weight=phylogeny_repeat_weight,
            reference_metadata_path=reference_metadata_path,
            target_taxon_id=target_taxon_id,
            taxon_calibration_path=taxon_calibration_path, taxon_alpha=taxon_alpha,
            taxon_min_loci=taxon_min_loci,
            taxon_min_locus_fraction=taxon_min_locus_fraction,
            taxon_bootstrap_replicates=taxon_bootstrap_replicates,
            taxon_min_bootstrap_support=taxon_min_bootstrap_support,
            taxon_max_mean_placement_entropy=taxon_max_mean_placement_entropy,
            taxon_min_median_placement_lwr=taxon_min_median_placement_lwr,
            taxon_identification=taxon_identification, taxon_k=taxon_k,
            taxon_minimum_margin=taxon_minimum_margin, input_mode="illumina",
            locus_quality={
                str(row["locus_id"]): {
                    "depth": row["primary_read_depth"],
                    "consensus_strength": row["allele_confidence"],
                    "status": row["status"],
                }
                for row in calls
            },
        )
    phylogenetic_rows = (
        read_profiles(phylogeny_paths["combined_marker_matches"])
        if phylogeny_paths else []
    )
    closest_reference_bands = (
        read_profiles(phylogeny_paths["closest_reference_bands"])
        if phylogeny_paths else []
    )
    write_tsv(
        matches + sequence_reference_match_rows(phylogenetic_rows),
        output / "profile_matches.tsv",
        MATCH_FIELDS,
    )

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
        "reference_source": "minimap2_competitive_candidate_contexts",
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
    from .minimap_mapping import minimap2_version
    metadata_path = output / "short_read_run_metadata.json"
    metadata_path.write_text(json.dumps({
                       "schema_version": "2.0", "mlvamaps_version": __version__,
                       "method": "competitive_minimap2_shared_inference",
                       "minimap2_version": minimap2_version(minimap2_bin),
                       "database": database_path or "panel-derived",
                       "insert_size": asdict(insert), "parameters": {"minimum_mapq": short_min_mapping_quality,
                       "minimum_supporting_fragments": min_depth,
                       "minimum_spanning_pairs": short_min_spanning_pairs,
                       "confidence_threshold": short_confidence_threshold,
                       "maximum_candidate_repeat_count": short_max_candidate_repeat_count,
                       "secondary_alignments": short_consider_secondary}}, indent=2, sort_keys=True) + "\n")
    write_report(
        output, sample_id, allele_rows, loci, matches, profiles,
        phylogenetic_rows=phylogenetic_rows,
        closest_reference_bands=closest_reference_bands,
        presence_rows=recruitment, local_assembly_rows=[], short_read_rows=calls,
    )
    report_path = output / "report.html"
    report_path.write_text(report_path.read_text().replace(
        "MLVA analysis report", "MLVA analysis report · Method: competitive minimap2 shared inference", 1
    ))
    if show_progress:
        detected = len(loci) - states["no_evidence"]
        print(f"[{sample_id}] {detected}/{len(loci)} loci have mapping evidence")
        print(f"[{sample_id}] {states['called']} called, {states['ambiguous']} ambiguous, "
              f"{states['low_coverage']} low coverage, {states['no_evidence']} no evidence")
    if not keep_intermediates:
        for path in (filtered1, filtered2, orphans):
            path.unlink(missing_ok=True)
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
        **phylogeny_paths,
        **unified_paths,
    }
