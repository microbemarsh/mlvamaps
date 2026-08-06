from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional


@dataclass(frozen=True)
class ReadRecord:
    read_id: str
    sequence: str
    quality: Optional[str] = None


@dataclass(frozen=True)
class ReadPair:
    """One sequencing molecule while retaining mate identity and orphans."""

    molecule_id: str
    read1: ReadRecord
    read2: Optional[ReadRecord] = None


@dataclass(frozen=True)
class Locus:
    locus_id: str
    chrom_or_contig: str = ""
    start: Optional[int] = None
    end: Optional[int] = None
    forward_primer: str = ""
    reverse_primer: str = ""
    left_flank_sequence: str = ""
    right_flank_sequence: str = ""
    repeat_motif: str = ""
    expected_min_repeats: int = 0
    expected_max_repeats: int = 100
    expected_amplicon_min_bp: int = 0
    expected_amplicon_max_bp: int = 100000
    pool_id: str = ""
    repeat_unit_length_bp: int = 0
    expected_product_size_bp: int = 0
    nominal_repeat_units: int = 0


@dataclass(frozen=True)
class Assignment:
    read_id: str
    sample_id: str
    assigned_locus: str
    assignment_score: float
    orientation: str
    primer_forward_detected: bool
    primer_reverse_detected: bool
    passes_assignment_qc: bool
    oriented_sequence: str
    oriented_quality: Optional[str]
    forward_start: Optional[int]
    forward_end: Optional[int]
    reverse_start: Optional[int]
    reverse_end: Optional[int]
    forward_mismatches: Optional[int]
    reverse_mismatches: Optional[int]
    product_size_bp: Optional[int] = None
    measurement_status: str = ""
    failure_reason: Optional[str] = None
    recruitment_reference: str = ""
    recruitment_reference_allele: int | float | str | None = None
    recruitment_mapq: Optional[int] = None
    recruitment_alignment_score: Optional[int] = None
    recruitment_cigar: str = ""
    extracted_query_start: Optional[int] = None
    extracted_query_end: Optional[int] = None
    locus_measurement: Optional["LocusMeasurement"] = None


@dataclass(frozen=True)
class AnchorMeasurement:
    start: int
    end: int
    identity: float
    edit_distance: int
    mismatches: int
    insertions: int
    deletions: int
    complete: bool = True


@dataclass(frozen=True)
class LocusMeasurement:
    locus_id: str
    sequence_id: Optional[str]
    source: str
    product_start: Optional[int]
    product_end: Optional[int]
    forward_anchor_start: Optional[int]
    forward_anchor_end: Optional[int]
    reverse_anchor_start: Optional[int]
    reverse_anchor_end: Optional[int]
    repeat_start: Optional[int]
    repeat_end: Optional[int]
    repeat_length_bp: Optional[int]
    repeat_unit_length: Optional[int]
    raw_repeat_count: Optional[float]
    called_allele: int | float | str | None
    forward_anchor_identity: Optional[float]
    reverse_anchor_identity: Optional[float]
    repeat_motif_identity: Optional[float]
    allele_likelihoods: Optional[dict[int | float, float]]
    confidence: Optional[float]
    status: str
    failure_reason: Optional[str]
    forward_anchor: Optional[AnchorMeasurement] = None
    reverse_anchor: Optional[AnchorMeasurement] = None
    credible_alleles: tuple[int | float, ...] = ()
    second_allele: int | float | str | None = None
    second_allele_probability: Optional[float] = None
    product_sequence: str = ""
    product_quality: Optional[str] = None
    repeat_sequence: str = ""
    measurement_method: str = ""
    metadata: Optional[dict[str, Any]] = None


@dataclass(frozen=True)
class RepeatFeature:
    read_id: str
    locus_id: str
    repeat_region_start: int
    repeat_region_end: int
    repeat_region_length_bp: int
    repeat_motif: str
    raw_repeat_count_estimate: float
    nearest_integer_repeat_count: int
    flank_quality_score: float
    repeat_pattern: str
    repeat_sequence: str
    mean_qscore: float
    mismatch_count_in_repeat_region: int
    motif_kmer_count: int
    left_primer_score: float
    right_primer_score: float
    left_flank_score: float
    right_flank_score: float
    amplicon_sequence: str = ""
    amplicon_quality: Optional[str] = None
    product_size_bp: int = 0
    repeat_measurement_method: str = "repeat_region_length"


@dataclass(frozen=True)
class ReadPrediction:
    read_id: str
    locus_id: str
    predicted_repeat_count: int | float
    probability: float
    top_alt_repeat_count: Optional[int | float]
    top_alt_probability: float
    variant_id: str
    insertions_vs_representative: int
    deletions_vs_representative: int
    substitutions_vs_representative: int
    evidence_weight: float
    raw_repeat_count_estimate: Optional[float] = None
    measurement_sigma: Optional[float] = None
    measurement_repeat_count_estimate: Optional[float] = None
