"""Technology-neutral competitive alignment and molecule evidence models."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class CandidateAlignment:
    molecule_id: str
    read_id: str
    mate: int | None
    locus_id: str
    candidate_id: str
    repeat_count: int | float
    reference_id: str
    alignment_score: float
    mapping_quality: int
    alignment_identity: float
    query_coverage: float
    reference_coverage: float
    query_start: int
    query_end: int
    reference_start: int
    reference_end: int
    cigar: str
    cs: str
    primary: bool
    secondary: bool
    reverse: bool
    query_sequence: str = ""
    query_quality: str | None = None
    template_length: int = 0
    next_reference_start: int = -1


@dataclass(frozen=True)
class CandidateEvidence:
    locus_id: str
    repeat_count: int | float
    molecule_id: str
    alignment_score: float = 0.0
    alignment_identity: float = 0.0
    locus_confidence: float = 0.0
    direct_product_measurement: bool = False
    full_repeat_span: bool = False
    left_boundary_span: bool = False
    right_boundary_span: bool = False
    repeat_indel_support: float = 0.0
    pair_geometry_support: float = 0.0
    technology: str = ""
    measured_repeat_count: int | float | None = None
    candidate_id: str = ""
    reference_id: str = ""
    mapping_quality: int = 0
    query_coverage: float = 0.0
    reference_coverage: float = 0.0
    cigar: str = ""
    cs: str = ""
    evidence_tier: str = "generic_locus_mapping"
    metadata: dict[str, object] = field(default_factory=dict)


EVIDENCE_FIELDS = [
    "sample", "locus", "molecule", "technology", "candidate_repeat",
    "candidate_id", "reference_id", "alignment_score", "identity",
    "mapping_quality", "query_coverage", "reference_coverage",
    "direct_product", "full_repeat_span", "left_boundary", "right_boundary",
    "repeat_indel_delta", "pair_geometry_score", "evidence_tier",
    "inferred_molecule_repeat", "cigar", "cs",
]


def evidence_row(
    sample_id: str,
    evidence: CandidateEvidence,
    inferred_repeat: int | float | None = None,
) -> dict[str, object]:
    return {
        "sample": sample_id,
        "locus": evidence.locus_id,
        "molecule": evidence.molecule_id,
        "technology": evidence.technology,
        "candidate_repeat": evidence.repeat_count,
        "candidate_id": evidence.candidate_id,
        "reference_id": evidence.reference_id,
        "alignment_score": round(evidence.alignment_score, 6),
        "identity": round(evidence.alignment_identity, 6),
        "mapping_quality": evidence.mapping_quality,
        "query_coverage": round(evidence.query_coverage, 6),
        "reference_coverage": round(evidence.reference_coverage, 6),
        "direct_product": "yes" if evidence.direct_product_measurement else "no",
        "full_repeat_span": "yes" if evidence.full_repeat_span else "no",
        "left_boundary": "yes" if evidence.left_boundary_span else "no",
        "right_boundary": "yes" if evidence.right_boundary_span else "no",
        "repeat_indel_delta": round(evidence.repeat_indel_support, 6),
        "pair_geometry_score": round(evidence.pair_geometry_support, 6),
        "evidence_tier": evidence.evidence_tier,
        "inferred_molecule_repeat": "" if inferred_repeat is None else inferred_repeat,
        "cigar": evidence.cigar,
        "cs": evidence.cs,
    }