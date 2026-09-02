"""Shared Illumina helpers and the canonical Bowtie2 short-read entry point."""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

from .concurrency import DEFAULT_THREADS
from .models import Locus, ReadPair, ReadRecord
from .pipeline import SIMPLE_CALL_FIELDS
from .sequence import mean_qscore, revcomp


SHORT_CALL_EXTRA_FIELDS = [
    "read_technology", "evidence_class", "informative_molecule_count",
    "boundary_1_support", "boundary_2_support", "both_boundary_support",
    "repeat_count_min", "repeat_count_max",
    "repeat_count_interval_reason", "confidence_reason", "short_read_warning",
    "recruited_read_pairs",
    "failure_reason", "primary_allele_support", "secondary_allele_support",
    "informative_allele_reads", "uninformative_locus_reads",
    "estimated_primary_fraction", "estimated_secondary_fraction", "mixture_status",
    "evidence_sources", "mapping_state", "supporting_fragments",
    "proper_spanning_pairs", "junction_read_count", "full_spanning_read_count",
    "cigar_indel_read_count", "mean_mapq", "median_mapq",
    "candidate_allele_scores", "best_reference_contexts",
    "reference_context_taxa", "reference_context_provenance", "mlva_method",
]
SHORT_CALL_FIELDS = SIMPLE_CALL_FIELDS + SHORT_CALL_EXTRA_FIELDS
SHORT_QC_FIELDS = ["sample_id", "metric", "value"]
SAMPLE_SUMMARY_FIELDS = [
    "sample_id", "input_read_1", "input_read_2", "read_technology", "sample_mode",
    "total_reads", "total_read_pairs", "retained_reads", "retained_pairs",
    "callable_loci", "complete_loci", "partial_loci", "presence_only_loci",
    "mixed_loci", "missing_loci", "best_profile_id", "best_profile_distance",
    "profile_confidence", "run_status", "warnings",
]


def _trim_read(read: ReadRecord, trim_quality: int) -> ReadRecord:
    if trim_quality <= 0 or read.quality is None:
        return read
    end = len(read.sequence)
    while end and ord(read.quality[end - 1]) - 33 < trim_quality:
        end -= 1
    return ReadRecord(read.read_id, read.sequence[:end], read.quality[:end])


def qc_read_pairs(
    pairs: list[ReadPair], min_length: int, min_mean_quality: float,
    trim_quality: int, min_pair_retention: float,
) -> tuple[list[ReadPair], dict[str, int]]:
    if min_length < 1 or min_mean_quality < 0 or trim_quality < 0:
        raise ValueError("short-read QC thresholds must be non-negative")
    if not 0 <= min_pair_retention <= 1:
        raise ValueError("short-min-pair-retention must be between 0 and 1")
    retained: list[ReadPair] = []
    metrics: Counter[str] = Counter()
    for pair in pairs:
        metrics["input_pairs"] += 1
        mates = [pair.read1] + ([pair.read2] if pair.read2 is not None else [])
        metrics["input_reads"] += len(mates)
        passing = [
            read for read in (_trim_read(item, trim_quality) for item in mates)
            if len(read.sequence) >= min_length
            and mean_qscore(read.quality) >= min_mean_quality
        ]
        if not passing or len(passing) / len(mates) < min_pair_retention:
            metrics["rejected_pairs"] += 1
            metrics["rejected_reads"] += len(mates)
            continue
        read1 = passing[0]
        read2 = passing[1] if len(passing) == 2 else None
        if len(mates) == 2 and len(passing) == 1:
            metrics["orphan_reads"] += 1
        retained.append(ReadPair(pair.molecule_id, read1, read2))
        metrics["retained_pairs"] += 1
        metrics["retained_reads"] += len(passing)
    return retained, dict(metrics)


def merge_read_pair(
    pair: ReadPair, min_overlap: int = 20, max_mismatch_fraction: float = 0.03,
) -> ReadRecord | None:
    """Merge a directly overlapping pair; this does not perform assembly."""
    if pair.read2 is None:
        return None
    left = pair.read1
    right_sequence = revcomp(pair.read2.sequence)
    right_quality = pair.read2.quality[::-1] if pair.read2.quality else "I" * len(right_sequence)
    left_quality = left.quality or "I" * len(left.sequence)
    overlap = next((
        size for size in range(min(len(left.sequence), len(right_sequence)), min_overlap - 1, -1)
        if sum(a != b for a, b in zip(left.sequence[-size:], right_sequence[:size])) / size
        <= max_mismatch_fraction
    ), None)
    if overlap is None:
        return None
    sequence = list(left.sequence)
    quality = list(left_quality)
    start = len(left.sequence) - overlap
    for offset in range(overlap):
        index = start + offset
        if sequence[index] != right_sequence[offset] and right_quality[offset] > quality[index]:
            sequence[index] = right_sequence[offset]
        quality[index] = max(quality[index], right_quality[offset])
    sequence.extend(right_sequence[overlap:])
    quality.extend(right_quality[overlap:])
    return ReadRecord(f"{pair.molecule_id}/merged", "".join(sequence), "".join(quality))


def _merge_overlap_is_repeat_only(pair: ReadPair, merged: ReadRecord, locus: Locus) -> bool:
    if pair.read2 is None or not locus.repeat_motif:
        return False
    overlap = len(pair.read1.sequence) + len(pair.read2.sequence) - len(merged.sequence)
    if overlap <= 0:
        return False
    motif = locus.repeat_motif.upper()
    if not motif or set(motif) - set("ACGT"):
        return False
    periodic = motif * math.ceil((overlap + len(motif) * 2) / len(motif))
    return pair.read1.sequence[-overlap:].upper() in periodic


def _allele_rows(call_rows: list[dict]) -> list[dict]:
    return [{
        "sample_id": row["sample_id"], "locus_id": row["locus_id"],
        "called_repeat_count": row["repeat_count"],
        "posterior_probability": row["allele_confidence"],
        "second_best_repeat_count": row["second_best_repeat_count"],
        "second_best_posterior": row["second_best_probability"],
        "read_depth": row["read_depth"], "primary_read_depth": row["primary_read_depth"],
        "num_vntr_asvs": row["num_candidate_variants"],
        "num_meaningful_variants": 1 + row["num_confirmed_secondary_variants"] if row["repeat_count"] != "" else 0,
        "num_candidate_variants": row["num_candidate_variants"],
        "num_confirmed_secondary_variants": row["num_confirmed_secondary_variants"],
        "dominant_vntr_asv": row["repeat_count"],
        "dominant_variant_fraction": row["estimated_primary_fraction"],
        "secondary_alleles": row["secondary_alleles"],
        "allele_distribution": row["allele_distribution"], "call_status": row["status"],
        "primary_product_size_bp": row["product_size_bp"],
        "primary_repeat_count_raw": row["repeat_count_raw"],
        "primary_measurement_source": row["evidence_sources"],
        "evidence_status": row["evidence_class"],
    } for row in call_rows]


def run_short_read_call(
    reads1_path: str, reads2_path: str | None, loci_path: str | None,
    outdir: str, sample_id: str, primers_path: str | None = None,
    profiles_path: str | None = None, database_path: str | None = None,
    sample_metadata: dict[str, str] | None = None, short_min_read_length: int = 40,
    short_min_mean_quality: float = 15.0, short_trim_quality: int = 0,
    short_min_pair_retention: float = 0.5, min_depth: int = 3,
    threads: int = DEFAULT_THREADS, keep_intermediates: bool = False,
    sample_mode: str = "isolate", bowtie2_bin: str = "bowtie2",
    bowtie2_build_bin: str = "bowtie2-build", short_min_mapping_quality: int = 0,
    short_min_spanning_pairs: int = 2, short_confidence_threshold: float = 0.8,
    short_max_candidate_repeat_count: int = 100, short_consider_secondary: bool = True,
    show_progress: bool = True,
) -> dict[str, Path]:
    """Call Illumina data with the sole supported Bowtie2 context algorithm."""
    from .short_read_mapping import run_mapping_short_read_call

    return run_mapping_short_read_call(
        reads1_path=reads1_path, reads2_path=reads2_path, loci_path=loci_path,
        primers_path=primers_path, profiles_path=profiles_path,
        database_path=database_path, outdir=outdir, sample_id=sample_id,
        sample_metadata=sample_metadata, short_min_read_length=short_min_read_length,
        short_min_mean_quality=short_min_mean_quality, short_trim_quality=short_trim_quality,
        short_min_pair_retention=short_min_pair_retention, min_depth=min_depth,
        threads=threads, keep_intermediates=keep_intermediates, sample_mode=sample_mode,
        bowtie2_bin=bowtie2_bin, bowtie2_build_bin=bowtie2_build_bin,
        short_min_mapping_quality=short_min_mapping_quality,
        short_min_spanning_pairs=short_min_spanning_pairs,
        short_confidence_threshold=short_confidence_threshold,
        short_max_candidate_repeat_count=short_max_candidate_repeat_count,
        short_consider_secondary=short_consider_secondary, show_progress=show_progress,
    )