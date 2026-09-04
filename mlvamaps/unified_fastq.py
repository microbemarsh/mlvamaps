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
from .io import write_fasta, write_tsv
from .long_read_evidence import extract_long_read_evidence
from .minimap_mapping import map_reads_to_candidates_bam
from .models import Locus
from .short_read_evidence import extract_short_read_evidence


def taxonomic_query_sequences(
    calls: list[dict[str, object]],
    evidence: list[CandidateEvidence],
    molecule_calls: dict[tuple[str, str], int | float],
    contexts: list[CandidateContext],
) -> dict[str, str]:
    """Choose the best-supported reference-guided marker for each FASTQ call."""
    context_by_id = {context.candidate_id: context for context in contexts}
    sequences: dict[str, str] = {}
    for call in calls:
        locus_id = str(call["locus"])
        repeat_count = call.get("repeat_count", "")
        if call.get("status") not in {"called", "low_coverage", "mixed"} or repeat_count == "":
            continue
        scores: dict[str, tuple[set[str], float]] = {}
        for row in evidence:
            if (
                row.locus_id != locus_id
                or molecule_calls.get((locus_id, row.molecule_id)) != repeat_count
                or row.candidate_id not in context_by_id
                or (
                    row.technology == "illumina"
                    and float(row.metadata.get("background_alignment_margin", 0)) <= 0
                )
            ):
                continue
            molecules, score = scores.setdefault(row.candidate_id, (set(), 0.0))
            molecules.add(row.molecule_id)
            scores[row.candidate_id] = (molecules, score + row.alignment_score)
        if scores:
            ranked = sorted(
                scores,
                key=lambda item: (-len(scores[item][0]), -scores[item][1], item),
            )
            candidate_id = ranked[0]
            if len(ranked) > 1 and len(scores[ranked[0]][0]) == len(scores[ranked[1]][0]):
                continue
            sequences[locus_id] = context_by_id[candidate_id].sequence
    return sequences


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
    database = Path(database_path) if database_path else None
    if database is not None and (database / "database").is_dir():
        database = database / "database"
    resource = (
        database / "competitive_mapping"
        if database is not None and (database / "competitive_mapping").is_dir()
        else database
    )
    cached_fasta = resource / "candidate_contexts.fasta" if resource is not None else None
    if cached_fasta is not None and cached_fasta.is_file():
        paths = {
            "fasta": cached_fasta,
            "metadata": resource / "candidate_metadata.tsv",
            "provenance": resource / "candidate_provenance.json",
        }
    else:
        paths = write_candidate_contexts(contexts, work)
    bam = work / "candidate_alignments.bam"
    index_name = "short.mmi" if technology == "illumina" else "long.mmi"
    cached_index = resource / index_name if resource is not None else None
    mapping_reference = cached_index if cached_index is not None and cached_index.is_file() else paths["fasta"]
    alignments = map_reads_to_candidates_bam(
        mapping_reference, reads1, reads2, contexts, bam, threads, technology,
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
    taxonomic_queries = outdir / "taxonomic_query_sequences.fasta"
    if database_path:
        write_fasta(
            taxonomic_query_sequences(calls, evidence, molecule_calls, contexts).items(),
            taxonomic_queries,
        )
    if not keep_alignments:
        bam.unlink(missing_ok=True)
    return calls, evidence, molecule_calls, {
        "common_locus_calls": common_calls,
        "molecule_evidence": evidence_path,
        "candidate_contexts": paths["fasta"],
        "candidate_metadata": paths["metadata"],
        "candidate_provenance": paths["provenance"],
        **({"taxonomic_query_sequences": taxonomic_queries} if database_path else {}),
        **({"candidate_alignments": bam} if keep_alignments else {}),
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
