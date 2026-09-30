"""Illumina orchestration using shared competitive alignment evidence.

The module deliberately separates locus detection from genotype resolution.  A
few flank-anchored molecules can establish presence without supplying enough
information to choose between adjacent repeat counts.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

from . import __version__
from .io import open_text, read_profiles, write_tsv
from .profile_matching import (
    PROFILE_MATCH_LOCUS_FIELDS,
    build_fingerprint,
    match_profiles,
    profile_match_locus_rows,
    sequence_reference_match_rows,
)


MAPPING_EVIDENCE_FIELDS = [
    "sample_id", "locus_id", "state", "repeat_count", "locus_length_bp",
    "confidence", "supporting_fragments", "proper_spanning_pairs",
    "junction_reads", "full_spanning_reads", "cigar_indel_reads",
    "mean_mapq", "median_mapq", "candidate_scores", "best_contexts",
    "context_taxa", "reason",
]


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
    reference_metadata_path: str | None,
    taxon_min_loci: int | None,
    taxon_identification: bool | None,
    missing_locus_min_depth: float = 3.0,
    missing_locus_min_fraction: float = 0.8,
    missing_locus_penalty: float = 8.0,
    show_progress: bool = True,
    classification_repeat_scale: float = 1.0,
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

    classification_paths: dict[str, Path] = {}
    if database_path:

        if show_progress:
            print(
                f"[{sample_id}] Classifying Illumina molecules against observed reference VNTRs"
            )
        from .mapping_classification import run_mapping_classification
        classification_paths.update(run_mapping_classification(
            database_path=database_path, loci=loci, outdir=output, sample_id=sample_id,
            reads1=filtered1, reads2=filtered2 if reads2_path else None, technology="illumina",
            sample_mode=sample_mode,
            threads=thread_count, minimap2_bin=minimap2_bin,
            locus_quality={str(row["locus"]): {"depth": row["molecule_support"], "status": row["status"]} for row in common_calls},
            query_repeat_counts={str(row["locus_id"]): row.get("repeat_count", "") for row in calls if row.get("status") in {"PASS", "LOW_DEPTH", "PRESENT"}},
            reference_metadata_path=reference_metadata_path,
            taxon_identification=taxon_identification, minimum_loci=taxon_min_loci or 2,
            repeat_scale=classification_repeat_scale,
            missing_locus_min_depth=missing_locus_min_depth,
            missing_locus_min_fraction=missing_locus_min_fraction,
            missing_locus_penalty=missing_locus_penalty,
        ))
    reference_rows = (
        read_profiles(classification_paths["mapping_reference_matches"])
        if classification_paths else []
    )
    closest_reference_bands = (
        read_profiles(classification_paths["closest_reference_bands"])
        if classification_paths else []
    )
    write_tsv(
        matches + sequence_reference_match_rows(reference_rows),
        output / "profile_matches.tsv",
        MATCH_FIELDS,
    )

    qc_values = dict(counters)
    qc_values.update({
        "insert_size_median": "", "insert_size_mean": "",
        "insert_size_standard_deviation": "", "insert_size_mad": "",
        "insert_size_pairs": 0,
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
                       "insert_size": {"median": None, "mean": None, "standard_deviation": None, "mad": None, "pairs_used": 0}, "parameters": {"minimum_mapq": short_min_mapping_quality,
                       "minimum_supporting_fragments": min_depth,
                       "minimum_spanning_pairs": short_min_spanning_pairs,
                       "confidence_threshold": short_confidence_threshold,
                       "maximum_candidate_repeat_count": short_max_candidate_repeat_count,
                       "secondary_alignments": short_consider_secondary}}, indent=2, sort_keys=True) + "\n")
    write_report(
        output, sample_id, allele_rows, loci, matches, profiles,
        reference_rows=reference_rows,
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
        **classification_paths,
        **unified_paths,
    }
