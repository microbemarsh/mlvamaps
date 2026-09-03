"""Orchestration helpers for the unified FASTQ allele-calling architecture."""

from __future__ import annotations

from pathlib import Path

from .alignment_evidence import CandidateEvidence, EVIDENCE_FIELDS, evidence_row
from .allele_inference import (
    COMMON_LOCUS_CALL_FIELDS,
    InferenceThresholds,
    infer_alleles,
)
from .candidate_contexts import generate_candidate_contexts, write_candidate_contexts
from .io import write_tsv
from .long_read_evidence import extract_long_read_evidence
from .minimap_mapping import map_reads_to_candidates
from .models import Locus
from .short_read_evidence import extract_short_read_evidence


def run_unified_fastq_inference(
    *,
    reads1: str | Path,
    reads2: str | Path | None,
    loci: list[Locus],
    database_path: str | Path | None,
    outdir: str | Path,
    sample_id: str,
    technology: str,
    minimap2_bin: str,
    threads: int,
    minimum_molecules: int,
    minimum_probability: float,
    maximum_candidate_repeat_count: int = 100,
    keep_alignments: bool = False,
) -> tuple[list[dict[str, object]], list[CandidateEvidence], dict[tuple[str, str], int | float], dict[str, Path]]:
    outdir = Path(outdir)
    work = outdir / "candidate_mapping"
    work.mkdir(parents=True, exist_ok=True)
    contexts = generate_candidate_contexts(
        loci, database_path, maximum=maximum_candidate_repeat_count
    )
    paths = write_candidate_contexts(contexts, work)
    sam = work / "candidate_alignments.sam"
    database = Path(database_path) if database_path else None
    if database is not None and (database / "database").is_dir():
        database = database / "database"
    resource = (
        database / "competitive_mapping"
        if database is not None and (database / "competitive_mapping").is_dir()
        else database
    )
    index_name = "short.mmi" if technology == "illumina" else "long.mmi"
    cached_index = resource / index_name if resource is not None else None
    mapping_reference = cached_index if cached_index is not None and cached_index.is_file() else paths["fasta"]
    alignments = map_reads_to_candidates(
        mapping_reference, reads1, reads2, contexts, sam, threads, technology,
        executable=minimap2_bin,
    )
    if technology == "illumina":
        evidence = extract_short_read_evidence(alignments, contexts, loci)
    else:
        evidence = extract_long_read_evidence(alignments, contexts, loci, technology)
    calls, molecule_calls = infer_alleles(
        evidence, loci, contexts, sample_id, technology,
        InferenceThresholds(
            minimum_molecules=minimum_molecules,
            minimum_probability=minimum_probability,
        ),
    )
    common_calls = outdir / "common_locus_calls.tsv"
    evidence_path = outdir / "molecule_candidate_evidence.tsv"
    write_tsv(calls, common_calls, COMMON_LOCUS_CALL_FIELDS)
    write_tsv(
        (evidence_row(sample_id, row, molecule_calls.get((row.locus_id, row.molecule_id))) for row in evidence),
        evidence_path,
        EVIDENCE_FIELDS,
    )
    if not keep_alignments:
        sam.unlink(missing_ok=True)
    return calls, evidence, molecule_calls, {
        "common_locus_calls": common_calls,
        "molecule_evidence": evidence_path,
        "candidate_contexts": paths["fasta"],
        "candidate_metadata": paths["metadata"],
        "candidate_provenance": paths["provenance"],
        **({"candidate_alignments": sam} if keep_alignments else {}),
    }


def common_calls_to_compatibility(calls: list[dict[str, object]]) -> list[dict[str, object]]:
    """Project shared calls into the established compact ``calls.tsv`` schema."""
    statuses = {
        "called": "PASS",
        "low_coverage": "LOW_DEPTH",
        "detected_unresolved": "PRESENT_COUNT_UNKNOWN",
        "ambiguous": "AMBIGUOUS",
        "not_found": "NOT_FOUND",
        "mixed": "MULTIPLE_VARIANTS",
    }
    output = []
    for row in calls:
        repeat = row["repeat_count"]
        output.append({
            "sample_id": row["sample"],
            "locus_id": row["locus"],
            "present": "no" if row["status"] == "not_found" else "yes",
            "repeat_count": repeat,
            "repeat_count_raw": repeat,
            "product_size_bp": "",
            "read_depth": row["molecule_support"],
            "primary_read_depth": row["molecule_support"],
            "mean_coverage": "",
            "allele_confidence": row["best_probability"],
            "second_best_repeat_count": "",
            "second_best_probability": row["second_best_probability"],
            "inference_method": "shared_competitive_minimap2_inference",
            "dominant_variant_fraction": row["dominant_fraction"],
            "num_candidate_variants": len(str(row["candidate_distribution"]).split(";")) if row["candidate_distribution"] else 0,
            "num_confirmed_secondary_variants": 1 if row["status"] == "mixed" else 0,
            "secondary_alleles": row["secondary_repeat"],
            "allele_distribution": row["candidate_distribution"],
            "status": statuses[str(row["status"])],
            "evidence": (
                f"{row['molecule_support']} informative molecule(s); "
                f"{row['direct_product_support']} direct product; "
                f"{row['full_span_support']} full repeat span"
            ),
        })
    return output