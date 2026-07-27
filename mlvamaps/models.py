from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class ReadRecord:
    read_id: str
    sequence: str
    quality: Optional[str] = None


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
