"""Orchestration helpers for the unified FASTQ allele-calling architecture."""

from __future__ import annotations

from pathlib import Path
from dataclasses import replace
from time import perf_counter

from .alignment_evidence import CandidateAlignment, CandidateEvidence, EVIDENCE_FIELDS, evidence_row
from .allele_inference import (
    COMMON_LOCUS_CALL_FIELDS,
    InferenceThresholds,
    ESTIMATION_FIELDS,
    CONFIDENCE_FIELDS,
    estimate_unresolved_calls,
    finalize_repeat_count_calls,
    infer_alleles,
)
from .candidate_contexts import CandidateContext, generate_candidate_contexts, write_candidate_contexts
from .io import ReadSpool, normalize_read_id, write_fasta, write_tsv
from .long_read_evidence import extract_long_read_evidence, flank_error_model, long_read_measurement_metadata
from .minimap_mapping import map_reads_to_candidates_bam
from .models import Assignment, Locus
from .short_read_evidence import extract_short_read_evidence


def taxonomic_query_sequences(calls, products):
    """Export only reconstructed sample products consistent with the final call."""
    return {
        str(call["locus"]): products[str(call["locus"])][0]
        for call in calls
        if call["status"] in {"called", "low_coverage"}
        and str(call["locus"]) in products
        and call["repeat_count"] == products[str(call["locus"])][1]
    }


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
    round_tolerance: float = 0.25,
    contexts: list[CandidateContext] | None = None,
    recruited_fragments: dict | None = None,
    primer_assignments: list[Assignment] | None = None,
    rescue_reads: list | None = None,
    max_primer_mismatches: int = 2,
    reference_alignments: list[CandidateAlignment] | None = None,
    sample_mode: str = "metagenome",
    insert_mean: float | None = None,
    insert_sd: float | None = None,
    mixture_min_fraction: float = 0.2,
    mixture_min_molecules: int = 2,
    inference_limits: InferenceThresholds | None = None,
) -> tuple[list[dict[str, object]], list[CandidateEvidence], dict[tuple[str, str], int | float], dict[str, Path]]:
    """Map observed loci, measure molecules, and optionally retain mappings for typing.

    ``reference_alignments`` is filled in place so callers can reuse this exact
    observed-reference mapping for classification without mapping FASTQ twice.
    """
    outdir = Path(outdir)
    thresholds = replace(inference_limits or InferenceThresholds(),
        minimum_molecules=minimum_molecules, minimum_probability=minimum_probability,
        mixture_min_fraction=mixture_min_fraction, mixture_min_molecules=mixture_min_molecules)
    work = outdir / "candidate_mapping"
    work.mkdir(parents=True, exist_ok=True)
    if contexts is None:
        if database_path:
            from .mapping_classification import observed_reference_contexts

            _root, contexts, _members = observed_reference_contexts(database_path, loci, round_tolerance)
        else:
            contexts = generate_candidate_contexts(loci, maximum=maximum_candidate_repeat_count)
    paths = write_candidate_contexts(contexts, work)
    bam = work / "candidate_alignments.bam"
    unmapped_mates = {}
    started = perf_counter()
    alignments = map_reads_to_candidates_bam(
        paths["fasta"], reads1, reads2, contexts, bam, threads, technology,
        executable=minimap2_bin,
        max_secondary=max(100, len(contexts)),
        include_unknown_repeats=bool(database_path),
        retain_all_competitors=bool(database_path),
        unmapped_mates=unmapped_mates if technology == "illumina" else None,
        query_reads=rescue_reads if technology != "illumina" and hasattr(rescue_reads, "get") else None,
    )
    mapped_at = perf_counter()
    if database_path and reference_alignments is not None:
        reference_alignments.extend(alignments)
    from .mapped_reconstruction import apply_reconstructed_calls, recruit_mapping_fragments, resolve_query_products

    ambiguous_assignments = {}
    if technology == "illumina":
        recruited_fragments = recruit_mapping_fragments(
            alignments, contexts, recruited_fragments, unmapped_mates,
            ambiguous_assignments=ambiguous_assignments,
        )
        owners = {name: locus for locus, pool in recruited_fragments.items() for name in pool}
        inference_alignments = [row for row in alignments if owners.get(row.molecule_id) == row.locus_id]
        insert_calibration = None
        if insert_mean is None and insert_sd is None:
            from .short_read_evidence import estimate_insert_distribution
            insert_calibration = estimate_insert_distribution(inference_alignments, contexts)
            if insert_calibration is not None:
                insert_mean, insert_sd, _ = insert_calibration
        evidence = extract_short_read_evidence(
            inference_alignments, contexts, loci, round_tolerance=round_tolerance,
            recruited_fragments=recruited_fragments,
            insert_mean=insert_mean, insert_sd=insert_sd,
        )
        if insert_mean is not None:
            for row in evidence:
                if "pair_repeat_count" in row.metadata:
                    row.metadata["insert_calibration"] = "mapped_unique_flanks" if insert_calibration else "user_supplied"
                    row.metadata["insert_calibration_fragments"] = insert_calibration[2] if insert_calibration else ""
    else:
        pools = recruit_mapping_fragments(alignments, contexts, ambiguous_assignments=ambiguous_assignments)
        owners = {name: locus for locus, pool in pools.items() for name in pool}
        inference_alignments = [row for row in alignments if owners.get(row.molecule_id) == row.locus_id]
        evidence = extract_long_read_evidence(inference_alignments, contexts, loci, technology,
            round_tolerance=round_tolerance, limits=thresholds)
        if rescue_reads is not None:
            from .in_silico_pcr import read_pcr_results, run_in_silico_pcr_loci
            from .locus_assignment import assignments_from_pcr

            complete = {row.molecule_id for row in evidence if row.direct_product_measurement}
            rescue_fasta = outdir / "primer_rescue.fasta.gz"
            with ReadSpool() as unresolved:
                def rescue_records():
                    for read in rescue_reads:
                        if normalize_read_id(read.read_id)[0] not in complete:
                            unresolved.append(read)
                            yield read.read_id, read.sequence
                write_fasta(rescue_records(), rescue_fasta)
                if not unresolved:
                    rescued = []
                else:
                    pcr = run_in_silico_pcr_loci(rescue_fasta, loci, outdir / "in_silico_pcr",
                                                max_errors=max_primer_mismatches, threads=threads)
                    pcr_rows = read_pcr_results(pcr["stats"], pcr["products"])
                    # Failed rescue records have no sequence evidence to retain.
                    # Keep their compact assignments, releasing original bases.
                    rescued = assignments_from_pcr(unresolved, loci, pcr_rows, sample_id,
                                                   retain_unassigned_sequence=False)
                if primer_assignments is None:
                    primer_assignments = []
                primer_assignments.extend(rescued)
    if primer_assignments:
        from .locus_measurement import measure_locus_product

        by_locus = {locus.locus_id: locus for locus in loci}
        error_model = flank_error_model(inference_alignments, contexts) if technology in {"ont", "ont-hq"} else None
        refined = {locus.locus_id: {row.molecule_id for row in evidence
                   if row.locus_id == locus.locus_id and row.metadata.get("inference_tier") == 2}
                   for locus in loci}
        for assignment in primer_assignments:
            if not assignment.passes_assignment_qc or assignment.assigned_locus not in by_locus:
                continue
            measurement = measure_locus_product(
                assignment.oriented_sequence, by_locus[assignment.assigned_locus],
                assignment.oriented_quality, source="primer_recovery",
                round_tolerance=round_tolerance,
            )
            if measurement.called_allele is not None and measurement.status == "FULL_PRODUCT":
                molecule = normalize_read_id(assignment.read_id)[0]
                locus_id = assignment.assigned_locus
                metadata = long_read_measurement_metadata(measurement, by_locus[locus_id], error_model,
                    thresholds, molecule in refined[locus_id] or len(refined[locus_id]) < thresholds.maximum_refinement_reads)
                if metadata.get("inference_tier") == 2:
                    refined[locus_id].add(molecule)
                evidence.append(CandidateEvidence(
                    assignment.assigned_locus, measurement.called_allele,
                    molecule,
                    technology=technology, locus_confidence=1.0,
                    direct_product_measurement=True, full_repeat_span=True,
                    measured_repeat_count=measurement.called_allele,
                    evidence_tier="primer_recovery",
                    metadata={"product_sequence": measurement.product_sequence, **metadata},
                    raw_repeat_count=measurement.raw_repeat_count,
                    repeat_sigma=metadata.get("repeat_sigma"),
                    orientation=assignment.orientation,
                    quality_flags=tuple(flag for flag, present in (
                        ("boundary_ambiguous", metadata.get("boundary_ambiguous")),
                        ("refinement_limited", metadata.get("refinement_fallback")),
                    ) if present),
                ))
    for molecule, candidates in sorted(ambiguous_assignments.items()):
        for locus_id in candidates:
            evidence.append(CandidateEvidence(locus_id, None, molecule, technology=technology,
                candidate_loci=candidates, quality_flags=("ambiguous_locus",), evidence_tier="ambiguous_locus"))
    evidence_at = perf_counter()
    calls, molecule_calls = infer_alleles(
        evidence, loci, contexts, sample_id, technology,
        thresholds,
        require_measurement=bool(database_path),
        sample_mode=sample_mode, round_tolerance=round_tolerance,
    )
    fragment_paths = {}
    diagnostics = []
    products = []
    if technology == "illumina" and recruited_fragments is not None:
        from .short_read_fragments import resolve_fragment_calls

        calls, evidence, molecule_calls, diagnostics, products = resolve_fragment_calls(
            calls, evidence, molecule_calls, alignments, contexts, loci,
            recruited_fragments, sample_id, thresholds, round_tolerance,
            require_measurement=bool(database_path),
            sample_mode=sample_mode,
        )
        fragment_paths["short_read_diagnostics"] = outdir / "short_read_diagnostics.tsv"
        if keep_alignments:
            fragment_paths["short_read_microassembly"] = outdir / "short_read_microassembly.fasta"
            write_fasta(products, fragment_paths["short_read_microassembly"])
    query_products, reconstruction, evidence, changed = resolve_query_products(
        calls, evidence, loci, products, sample_id, thresholds, round_tolerance, max_primer_mismatches,
    )
    if changed:
        calls, molecule_calls = infer_alleles(
            evidence, loci, contexts, sample_id, technology, thresholds,
            require_measurement=bool(database_path), sample_mode=sample_mode,
            round_tolerance=round_tolerance,
        )
    estimate_unresolved_calls(calls, evidence, contexts, thresholds.temperature)
    apply_reconstructed_calls(calls, reconstruction, thresholds)
    finalize_repeat_count_calls(calls, evidence, loci, contexts, thresholds)
    inferred_at = perf_counter()
    if "short_read_diagnostics" in fragment_paths:
        for diagnostic, call in zip(diagnostics, calls):
            diagnostic["final_repeat_count"] = call["repeat_count"]
            if call.get("estimation_method"):
                diagnostic["call_method"] = call["estimation_method"]
        write_tsv(diagnostics, fragment_paths["short_read_diagnostics"], list(diagnostics[0]) if diagnostics else [])
    common_calls = outdir / "common_locus_calls.tsv"
    evidence_path = outdir / "molecule_candidate_evidence.tsv"
    from .evidence_reporting import write_evidence_reports
    supplemental_paths = write_evidence_reports(calls, evidence, loci, outdir, thresholds)
    if technology == "illumina":
        supplemental_paths["library_insert_distribution"] = outdir / "library_insert_distribution.tsv"
        library = dict(sample=sample_id, source="mapped_unique_flanks" if insert_calibration else "user_supplied" if insert_mean is not None else "unavailable",
            mean_bp=insert_mean if insert_mean is not None else "", sd_bp=insert_sd if insert_sd is not None else "",
            independent_fragments=insert_calibration[2] if insert_calibration else "",
            limitation="Mapping selection bias is not calibrated; paired-length genotypes remain estimates.")
        write_tsv([library], supplemental_paths["library_insert_distribution"], list(library))
    supplemental_paths["fastq_performance"] = outdir / "fastq_performance.tsv"
    performance = [
        {"metric": "competitive_mapping_seconds", "value": mapped_at - started},
        {"metric": "recruitment_and_extraction_seconds", "value": evidence_at - mapped_at},
        {"metric": "inference_and_reconstruction_seconds", "value": inferred_at - evidence_at},
        {"metric": "mapped_molecules", "value": len({row.molecule_id for row in alignments})},
        {"metric": "mapped_molecules_per_second", "value": len({row.molecule_id for row in alignments}) / max(1e-9, inferred_at - started)},
    ]
    for tier in (1, 2, 3):
        performance.append({"metric": f"tier_{tier}_locus_fraction", "value": sum(
            call.get("inference_tier") == tier for call in calls) / max(1, len(calls))})
    write_tsv(performance, supplemental_paths["fastq_performance"], ["metric", "value"])
    write_tsv(calls, common_calls, COMMON_LOCUS_CALL_FIELDS)
    write_tsv(
        (evidence_row(sample_id, row, molecule_calls.get((row.locus_id, row.molecule_id))) for row in evidence),
        evidence_path,
        EVIDENCE_FIELDS,
    )
    query_sequences = taxonomic_query_sequences(calls, query_products)
    query_path = outdir / "query_loci.fasta"
    write_fasta(query_sequences.items(), query_path)
    reconstruction_path = outdir / "query_locus_reconstruction.tsv"
    for row in reconstruction:
        row["exported"] = "yes" if row["locus_id"] in query_sequences else "no"
        if row["status"] == "RECONSTRUCTED" and row["exported"] == "no":
            row["status"] = "FINAL_CALL_DISAGREES_OR_UNRESOLVED"
    write_tsv(reconstruction, reconstruction_path, list(reconstruction[0]) if reconstruction else [])
    # Compatibility filename now contains sample sequence only. Rewrite even
    # when empty so an earlier reference proxy cannot survive a rerun.
    taxonomic_queries = outdir / "taxonomic_query_sequences.fasta"
    write_fasta(query_sequences.items(), taxonomic_queries)
    if not keep_alignments:
        bam.unlink(missing_ok=True)
    return calls, evidence, molecule_calls, {
        "query_loci": query_path,
        "query_locus_reconstruction": reconstruction_path,
        "common_locus_calls": common_calls,
        "molecule_evidence": evidence_path,
        "candidate_contexts": paths["fasta"],
        "candidate_metadata": paths["metadata"],
        "candidate_provenance": paths["provenance"],
        **({"taxonomic_query_sequences": taxonomic_queries} if database_path else {}),
        **({"candidate_alignments": bam} if keep_alignments else {}),
        **fragment_paths,
        **supplemental_paths,
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
        "estimated": "ESTIMATED",
    }
    output = []
    for row in calls:
        repeat = row["repeat_count"]
        output.append({
            "sample_id": row["sample"],
            "locus_id": row["locus"],
            "present": "no" if row["status"] == "not_found" else "yes",
            "repeat_count": repeat,
            "repeat_count_raw": row.get("repeat_count_raw", repeat),
            "product_size_bp": "",
            "read_depth": row.get("presence_support", row["molecule_support"]),
            "primary_read_depth": row["molecule_support"],
            "mean_coverage": "",
            "allele_confidence": row["best_probability"],
            "second_best_repeat_count": "",
            "second_best_probability": row["second_best_probability"],
            "inference_method": row.get("estimation_method") or row.get("inference_method", "shared_competitive_minimap2_inference"),
            **{field: row.get(field, "") for field in ESTIMATION_FIELDS},
            **{field: row.get(field, "") for field in CONFIDENCE_FIELDS + ["second_best_repeat_count"]},
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
                + (f"; low-confidence estimate: {row.get('estimation_reason', 'unresolved product')}; "
                   "model support is relative to tested candidates, not a calibrated probability of correctness"
                   if row["status"] == "estimated" else "")
            ),
        })
    return output
