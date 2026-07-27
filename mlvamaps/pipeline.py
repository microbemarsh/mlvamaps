from __future__ import annotations

from pathlib import Path

from .bayesian_caller import call_loci
from .clustering import cluster_vntr_asvs
from .concurrency import DEFAULT_THREADS, resolve_threads
from .in_silico_pcr import read_pcr_results, run_in_silico_pcr_loci
from .io import read_fastq, read_profiles, write_fasta, write_fastq, write_tsv
from .locus_assignment import assignments_from_pcr
from .mapping import (
    MAPPING_SUMMARY_FIELDS,
    SNP_FIELDS,
    run_locus_mapping,
)
from .mixture import MIXTURE_FIELDS, estimate_variant_mixtures
from .ml_classifier import predict_read_alleles
from .profile_matching import build_fingerprint, match_profiles
from .phylogeny import dominant_read_query_sequences, run_phylogenetic_placement
from .progress import ProgressReporter
from .primers import read_loci_or_primers
from .qc import filter_reads
from .repeat_parser import extract_repeat_features
from .report import write_report
from .recruitment import (
    RECRUITMENT_READ_FIELDS,
    RECRUITMENT_SUMMARY_FIELDS,
    local_product_records,
    recruitment_fallback_evidence,
    recruitment_summary_rows,
    run_read_recruitment,
)


ASSIGNMENT_FIELDS = [
    "read_id",
    "sample_id",
    "assigned_locus",
    "assignment_score",
    "orientation",
    "primer_forward_detected",
    "primer_reverse_detected",
    "passes_assignment_qc",
    "product_size_bp",
]

FEATURE_FIELDS = [
    "read_id",
    "locus_id",
    "repeat_region_start",
    "repeat_region_end",
    "repeat_region_length_bp",
    "repeat_motif",
    "raw_repeat_count_estimate",
    "nearest_integer_repeat_count",
    "flank_quality_score",
    "repeat_pattern",
    "repeat_sequence",
    "mean_qscore",
    "mismatch_count_in_repeat_region",
    "motif_kmer_count",
    "left_primer_score",
    "right_primer_score",
    "left_flank_score",
    "right_flank_score",
    "product_size_bp",
    "repeat_measurement_method",
]

ASV_FIELDS = [
    "sample_id",
    "locus_id",
    "variant_id",
    "repeat_count",
    "support_reads",
    "unique_sequences",
    "frequency",
    "representative_read_id",
    "representative_pattern",
    "representative_sequence",
    "representative_length_bp",
    "reads_with_indels",
    "total_insertions",
    "total_deletions",
    "total_substitutions",
    "mean_edit_distance_to_representative",
    "max_edit_distance_to_representative",
]

ASV_MEMBERSHIP_FIELDS = [
    "sample_id",
    "read_id",
    "locus_id",
    "variant_id",
    "repeat_count",
    "repeat_sequence",
    "aligned_repeat_sequence",
    "aligned_representative_sequence",
    "insertions_vs_representative",
    "deletions_vs_representative",
    "substitutions_vs_representative",
    "edit_distance_to_representative",
]

PREDICTION_FIELDS = [
    "read_id",
    "locus_id",
    "predicted_repeat_count",
    "probability",
    "top_alt_repeat_count",
    "top_alt_probability",
    "variant_id",
    "insertions_vs_representative",
    "deletions_vs_representative",
    "substitutions_vs_representative",
    "evidence_weight",
    "raw_repeat_count_estimate",
    "measurement_sigma",
    "measurement_repeat_count_estimate",
]

ALLELE_FIELDS = [
    "sample_id",
    "locus_id",
    "called_repeat_count",
    "posterior_probability",
    "second_best_repeat_count",
    "second_best_posterior",
    "read_depth",
    "effective_read_depth",
    "primary_read_depth",
    "primary_effective_read_depth",
    "confidence_effective_depth",
    "num_vntr_asvs",
    "num_meaningful_variants",
    "num_candidate_variants",
    "num_confirmed_secondary_variants",
    "dominant_vntr_asv",
    "dominant_variant_fraction",
    "secondary_alleles",
    "allele_distribution",
    "call_status",
    "sample_mode",
    "calling_convention",
]

MATCH_FIELDS = [
    "sample_id",
    "best_profile_id",
    "strain_id",
    "distance",
    "matched_loci",
    "mismatched_loci",
    "confidence",
    "compared_loci",
    "mean_negative_log_likelihood",
    "profile_probability_score",
]


SIMPLE_CALL_FIELDS = [
    "sample_id",
    "locus_id",
    "present",
    "repeat_count",
    "repeat_count_raw",
    "product_size_bp",
    "read_depth",
    "primary_read_depth",
    "mean_coverage",
    "allele_confidence",
    "second_best_repeat_count",
    "second_best_probability",
    "inference_method",
    "dominant_variant_fraction",
    "num_candidate_variants",
    "num_confirmed_secondary_variants",
    "secondary_alleles",
    "allele_distribution",
    "status",
    "evidence",
]

REPEAT_COUNT_FIELDS = [
    "sample_id",
    "locus_id",
    "repeat_count",
    "repeat_count_raw",
    "read_depth",
    "primary_read_depth",
    "allele_confidence",
    "dominant_variant_fraction",
    "secondary_alleles",
    "status",
]

ALLELE_DISTRIBUTION_FIELDS = [
    "sample_id",
    "locus_id",
    "allele",
    "probability",
    "rank",
    "selected",
    "inference_method",
]


def allele_distribution_rows(sample_id: str, call_rows: list[dict]) -> list[dict]:
    """Expand compact posterior strings into an analysis-friendly long table."""
    rows = []
    for call in call_rows:
        selected = call.get("called_repeat_count", call.get("repeat_count", ""))
        method = call.get("inference_method", "read_distribution")
        for rank, entry in enumerate(
            filter(None, str(call.get("allele_distribution", "")).split(";")),
            start=1,
        ):
            allele, probability = entry.rsplit(":", 1)
            rows.append(
                {
                    "sample_id": sample_id,
                    "locus_id": call["locus_id"],
                    "allele": allele,
                    "probability": probability,
                    "rank": rank,
                    "selected": "yes" if str(selected) == allele else "no",
                    "inference_method": method,
                }
            )
    return rows


def simple_call_rows_from_alleles(sample_id: str, allele_rows: list[dict]) -> list[dict]:
    rows = []
    for row in allele_rows:
        read_depth = int(row.get("read_depth") or 0)
        present = read_depth > 0
        repeat_count = row.get("called_repeat_count", "") if present else ""
        dominant_fraction = float(row.get("dominant_variant_fraction") or 0)
        if present:
            dominant = row.get("dominant_vntr_asv", "")
            meaningful = int(row.get("num_meaningful_variants") or 0)
            evidence = (
                f"{read_depth} assigned reads; "
                f"{int(row.get('primary_read_depth') or 0)} primary reads; "
                f"{meaningful} meaningful variant(s); dominant {dominant} at "
                f"{dominant_fraction:.1%}"
            )
        else:
            evidence = "no assigned reads"
        rows.append(
            {
                "sample_id": sample_id,
                "locus_id": row["locus_id"],
                "present": "yes" if present else "no",
                "repeat_count": repeat_count,
                "repeat_count_raw": repeat_count,
                "product_size_bp": "",
                "read_depth": read_depth,
                "primary_read_depth": int(row.get("primary_read_depth") or 0),
                "mean_coverage": "",
                "allele_confidence": row.get("posterior_probability", 0.0),
                "second_best_repeat_count": row.get("second_best_repeat_count", ""),
                "second_best_probability": row.get("second_best_posterior", 0.0),
                "inference_method": (
                    "assembly_equivalent_read_distribution"
                    if row.get("calling_convention") == "assembly"
                    else "read_distribution"
                ),
                "dominant_variant_fraction": dominant_fraction,
                "num_candidate_variants": int(
                    row.get("num_candidate_variants") or 0
                ),
                "num_confirmed_secondary_variants": int(
                    row.get("num_confirmed_secondary_variants") or 0
                ),
                "secondary_alleles": row.get("secondary_alleles", ""),
                "allele_distribution": row.get("allele_distribution", ""),
                "status": row["call_status"],
                "evidence": evidence,
            }
        )
    return rows


def run_call(
    reads_path: str,
    loci_path: str | None,
    outdir: str,
    sample_id: str,
    primers_path: str | None = None,
    profiles_path: str | None = None,
    database_path: str | None = None,
    recruitment_database_path: str | None = None,
    min_read_length: int = 50,
    max_read_length: int = 100000,
    min_qscore: float = 17.0,
    max_primer_mismatches: int = 3,
    min_depth: int = 10,
    min_posterior: float = 0.75,
    min_cluster_size: int = 1,
    cluster_min_identity: float = 0.97,
    min_mixture_fraction: float = 0.01,
    min_secondary_reads: int = 2,
    vsearch_bin: str = "vsearch",
    amplirust_bin: str = "amplirust",
    minimap2_bin: str = "minimap2",
    mafft_bin: str = "mafft",
    raxml_ng_bin: str = "raxml-ng",
    epa_ng_bin: str = "epa-ng",
    raxml_model: str = "DNA",
    phylogeny_snp_weight: float = 1.0,
    phylogeny_repeat_weight: float = 1.0,
    reference_metadata_path: str | None = None,
    locus_mapping: bool = True,
    min_mapping_quality: int = 0,
    min_base_quality: int = 20,
    min_snp_depth: int = 3,
    min_snp_alternate_reads: int = 2,
    min_snp_frequency: float = 0.2,
    threads: int = DEFAULT_THREADS,
    show_progress: bool = False,
    sample_mode: str = "metagenome",
    assembly_equivalent_reads: bool = True,
    assembly_round_tolerance: float = 0.25,
    max_confidence_depth: float = 25.0,
    fastq_strategy: str = "recruit",
    recruitment_preset: str | None = None,
    recruitment_min_identity: float = 0.9,
    recruitment_min_aligned_bp: int = 100,
    recruitment_min_locus_margin: int = 10,
) -> dict[str, Path]:
    if fastq_strategy not in {"recruit", "primer"}:
        raise ValueError("fastq_strategy must be 'recruit' or 'primer'")
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    progress = ProgressReporter(enabled=show_progress)
    thread_count = resolve_threads(threads)
    progress.step(f"Starting FASTQ call for sample {sample_id!r} with {thread_count} worker(s)")

    progress.step("Loading panel")
    loci = read_loci_or_primers(loci_path, primers_path)
    profiles = read_profiles(profiles_path)
    progress.step(f"Loaded {len(loci):,} loci" + (f" and {len(profiles):,} reference profiles" if profiles else ""))

    progress.step(f"Reading reads from {reads_path}")
    reads = []
    for idx, read in enumerate(read_fastq(reads_path), start=1):
        reads.append(read)
        progress.count("Read FASTQ records", idx)
    progress.count("Read FASTQ records", len(reads), force=True)

    progress.step("Filtering reads")
    filtered_reads, qc_rows = filter_reads(reads, min_read_length, max_read_length, min_qscore, progress)
    progress.step(f"Kept {len(filtered_reads):,}/{len(reads):,} reads after QC")
    write_tsv(qc_rows, outdir_path / "qc_summary.tsv", ["metric", "value"])
    write_fastq(filtered_reads, outdir_path / "filtered_reads.fastq.gz")

    recruitment_paths = {
        "recruited_reads": outdir_path / "locus_recruited_reads.tsv",
        "locus_presence": outdir_path / "locus_presence.tsv",
        "local_products": outdir_path / "local_locus_products.fasta",
        "recruitment_references": (
            outdir_path / "recruitment" / "locus_recruitment_references.fasta"
        ),
        "recruitment_alignments": (
            outdir_path / "recruitment" / "read_recruitment.sam"
        ),
    }
    recruited_assignments = []
    recruited_rows: list[dict] = []
    if fastq_strategy == "recruit":
        progress.step("Recruiting reads competitively to locus product references")
        (
            recruited_rows,
            _presence_rows,
            recruited_assignments,
            recruitment_paths,
        ) = run_read_recruitment(
            filtered_reads,
            loci,
            outdir_path,
            sample_id,
            recruitment_database_path or database_path,
            thread_count,
            executable=minimap2_bin,
            preset=recruitment_preset,
            min_mapping_quality=min_mapping_quality,
            min_alignment_identity=recruitment_min_identity,
            min_aligned_bp=recruitment_min_aligned_bp,
            min_locus_score_margin=recruitment_min_locus_margin,
        )
    else:
        write_tsv(
            [],
            recruitment_paths["recruited_reads"],
            RECRUITMENT_READ_FIELDS,
        )
        write_tsv(
            [],
            recruitment_paths["locus_presence"],
            RECRUITMENT_SUMMARY_FIELDS,
        )
        recruitment_paths["local_products"].write_text("")

    if filtered_reads:
        progress.step("Assigning reads with MLVA_finder-compatible Sassy primer matching")
        assignment_fasta = outdir_path / "filtered_reads.fasta"
        write_fasta(((read.read_id, read.sequence) for read in filtered_reads), assignment_fasta)
        pcr_paths = run_in_silico_pcr_loci(
            assignment_fasta,
            loci,
            outdir_path / "in_silico_pcr",
            max_errors=max_primer_mismatches,
            threads=threads,
        )
        primer_assignments = assignments_from_pcr(
            filtered_reads,
            loci,
            read_pcr_results(pcr_paths["stats"], pcr_paths["products"]),
            sample_id,
            progress=progress,
        )
    else:
        primer_assignments = []
    recruited_by_read = {
        assignment.read_id: assignment for assignment in recruited_assignments
    }
    assignments = [
        (
            assignment
            if assignment.passes_assignment_qc
            else recruited_by_read.get(assignment.read_id, assignment)
        )
        for assignment in primer_assignments
    ]
    primer_read_ids = {assignment.read_id for assignment in primer_assignments}
    assignments.extend(
        assignment
        for assignment in recruited_assignments
        if assignment.read_id not in primer_read_ids
    )
    recruited_read_ids = {str(row["read_id"]) for row in recruited_rows}
    for assignment in assignments:
        if (
            not assignment.passes_assignment_qc
            or assignment.read_id in recruited_read_ids
        ):
            continue
        recruited_rows.append(
            {
                "read_id": assignment.read_id,
                "locus_id": assignment.assigned_locus,
                "reference_name": "",
                "reference_source": "primer_fallback",
                "candidate_allele": "",
                "mapping_quality": "",
                "alignment_identity": assignment.assignment_score,
                "locus_score_margin": "",
                "aligned_query_bp": len(assignment.oriented_sequence),
                "reference_coverage": "",
                "full_product": "yes",
                "genotype_informative": "yes",
                "evidence_class": "FULL_PRODUCT",
            }
        )
    presence_rows = recruitment_summary_rows(sample_id, loci, recruited_rows)
    write_tsv(
        recruited_rows,
        recruitment_paths["recruited_reads"],
        RECRUITMENT_READ_FIELDS,
    )
    write_tsv(
        presence_rows,
        recruitment_paths["locus_presence"],
        RECRUITMENT_SUMMARY_FIELDS,
    )
    write_fasta(
        local_product_records(
            [
                assignment
                for assignment in assignments
                if assignment.passes_assignment_qc
            ]
        ),
        recruitment_paths["local_products"],
    )
    assignment_rows = [{field: getattr(row, field) for field in ASSIGNMENT_FIELDS} for row in assignments]
    write_tsv(assignment_rows, outdir_path / "read_locus_assignments.tsv", ASSIGNMENT_FIELDS)
    assigned_count = sum(1 for row in assignments if row.passes_assignment_qc)
    progress.step(f"Assigned {assigned_count:,}/{len(assignments):,} reads to primer-supported loci")

    features = extract_repeat_features(assignments, loci, threads=threads, progress=progress)
    feature_rows = [{field: getattr(row, field) for field in FEATURE_FIELDS} for row in features]
    write_tsv(feature_rows, outdir_path / "read_repeat_features.tsv", FEATURE_FIELDS)
    progress.step(f"Extracted {len(features):,} repeat feature records")

    vsearch_dir = outdir_path / "vsearch"
    progress.step(
        f"Clustering VNTR reads per locus with VSEARCH using {thread_count} thread(s)"
    )
    asv_rows, fasta_records, asv_memberships = cluster_vntr_asvs(
        features,
        loci,
        vsearch_dir,
        threads=thread_count,
        min_cluster_size=min_cluster_size,
        min_identity=cluster_min_identity,
        executable=vsearch_bin,
    )
    fallback_asvs, fallback_predictions = recruitment_fallback_evidence(
        recruited_rows,
        loci,
        {feature.locus_id for feature in features},
        sample_id,
    )
    asv_rows.extend(fallback_asvs)
    fasta_records.extend(
        (
            str(row["variant_id"]),
            str(row["representative_sequence"]),
        )
        for row in fallback_asvs
    )
    for row in asv_rows:
        row["sample_id"] = sample_id
    for row in asv_memberships:
        row["sample_id"] = sample_id
    write_tsv(asv_rows, outdir_path / "vntr_asv_table.tsv", ASV_FIELDS)
    write_tsv(
        asv_memberships,
        outdir_path / "vntr_asv_memberships.tsv",
        ASV_MEMBERSHIP_FIELDS,
    )
    write_fasta(fasta_records, outdir_path / "vntr_asv_representatives.fasta")

    progress.step("Estimating meaningful variant fractions with count-based EM")
    mixture_rows = estimate_variant_mixtures(
        asv_rows,
        min_fraction=min_mixture_fraction,
        min_secondary_reads=min_secondary_reads,
    )
    write_tsv(
        mixture_rows,
        outdir_path / "vntr_mixture_abundance.tsv",
        MIXTURE_FIELDS,
    )

    mapping_rows: list[dict] = []
    snp_rows: list[dict] = []
    if locus_mapping:
        progress.step(
            "Mapping locus reads to dominant VSEARCH representatives with minimap2"
        )
        mapping_rows, snp_rows, _mapping_paths = run_locus_mapping(
            features,
            asv_rows,
            outdir_path,
            sample_id,
            thread_count,
            executable=minimap2_bin,
            min_mapping_quality=min_mapping_quality,
            min_base_quality=min_base_quality,
            min_snp_depth=min_snp_depth,
            min_snp_alternate_reads=min_snp_alternate_reads,
            min_snp_frequency=min_snp_frequency,
        )
        progress.step(
            f"Mapped reads at {len(mapping_rows):,} loci and retained {len(snp_rows):,} SNP call(s)"
        )
    write_tsv(
        mapping_rows,
        outdir_path / "locus_mapping_summary.tsv",
        MAPPING_SUMMARY_FIELDS,
    )
    write_tsv(snp_rows, outdir_path / "locus_snps.tsv", SNP_FIELDS)

    progress.step("Calling repeat counts")
    predictions = predict_read_alleles(
        features,
        loci,
        asv_memberships,
        assembly_equivalent=assembly_equivalent_reads,
        assembly_round_tolerance=assembly_round_tolerance,
    )
    predictions.extend(fallback_predictions)
    prediction_rows = [{field: getattr(row, field) for field in PREDICTION_FIELDS} for row in predictions]
    write_tsv(prediction_rows, outdir_path / "read_level_allele_predictions.tsv", PREDICTION_FIELDS)

    allele_rows = call_loci(
        predictions,
        loci,
        asv_rows,
        min_depth,
        min_posterior,
        mixture_rows=mixture_rows,
        sample_mode=sample_mode,
        calling_convention=(
            "assembly" if assembly_equivalent_reads else "probabilistic"
        ),
        max_confidence_depth=max_confidence_depth,
    )
    for row in allele_rows:
        row["sample_id"] = sample_id
    write_tsv(allele_rows, outdir_path / "allele_calls.tsv", ALLELE_FIELDS)
    allele_distribution_path = outdir_path / "allele_probability_distribution.tsv"
    write_tsv(
        allele_distribution_rows(sample_id, allele_rows),
        allele_distribution_path,
        ALLELE_DISTRIBUTION_FIELDS,
    )
    simple_call_rows = simple_call_rows_from_alleles(sample_id, allele_rows)
    write_tsv(simple_call_rows, outdir_path / "calls.tsv", SIMPLE_CALL_FIELDS)
    write_tsv(simple_call_rows, outdir_path / "locus_repeat_counts.tsv", REPEAT_COUNT_FIELDS)

    fingerprint_rows, probabilistic_rows = build_fingerprint(sample_id, allele_rows, loci)
    fingerprint_fields = ["sample_id"] + [locus.locus_id for locus in loci]
    write_tsv(fingerprint_rows, outdir_path / "mlva_fingerprint.tsv", fingerprint_fields)
    write_tsv(
        probabilistic_rows,
        outdir_path / "mlva_fingerprint_probabilistic.tsv",
        ["sample_id", "locus_id", "repeat_count", "posterior_probability"],
    )

    match_rows = match_profiles(
        sample_id,
        fingerprint_rows[0],
        profiles,
        allele_rows=allele_rows,
    )
    write_tsv(match_rows, outdir_path / "profile_matches.tsv", MATCH_FIELDS)
    phylogeny_paths: dict[str, Path] = {}
    phylogenetic_rows: list[dict] = []
    closest_reference_bands: list[dict] = []
    if database_path:
        progress.step(
            "Placing MAFFT-aligned queries with EPA-ng using reusable reference trees when available"
        )
        phylogeny_paths = run_phylogenetic_placement(
            dominant_read_query_sequences(features, asv_rows),
            database_path,
            outdir_path,
            sample_id,
            loci,
            thread_count,
            mafft_bin=mafft_bin,
            raxml_ng_bin=raxml_ng_bin,
            epa_ng_bin=epa_ng_bin,
            raxml_model=raxml_model,
            snp_weight=phylogeny_snp_weight,
            repeat_weight=phylogeny_repeat_weight,
            reference_metadata_path=reference_metadata_path,
            progress=progress,
        )
        phylogenetic_rows = read_profiles(phylogeny_paths["combined_marker_matches"])
        closest_reference_bands = read_profiles(
            phylogeny_paths["closest_reference_bands"]
        )
    progress.step("Writing HTML report")
    write_report(
        outdir_path,
        sample_id,
        allele_rows,
        loci,
        match_rows,
        profiles,
        asv_rows,
        mapping_rows,
        snp_rows,
        mixture_rows,
        phylogenetic_rows,
        closest_reference_bands,
        presence_rows,
    )
    progress.step(f"Done. Main calls: {outdir_path / 'calls.tsv'}")

    return {
        "outdir": outdir_path,
        "calls": outdir_path / "calls.tsv",
        "allele_calls": outdir_path / "allele_calls.tsv",
        "repeat_counts": outdir_path / "locus_repeat_counts.tsv",
        "allele_distribution": allele_distribution_path,
        "asv_table": outdir_path / "vntr_asv_table.tsv",
        "asv_memberships": outdir_path / "vntr_asv_memberships.tsv",
        "asv_representatives": outdir_path / "vntr_asv_representatives.fasta",
        "mixture_abundance": outdir_path / "vntr_mixture_abundance.tsv",
        "mapping_summary": outdir_path / "locus_mapping_summary.tsv",
        "mapping_snps": outdir_path / "locus_snps.tsv",
        "mapping_references": outdir_path / "locus_mapping_references.fasta",
        "mapping_alignments": outdir_path / "locus_read_alignments.sam",
        "minimap2": outdir_path / "minimap2",
        "vsearch": vsearch_dir,
        "in_silico_pcr": outdir_path / "in_silico_pcr",
        "fingerprint": outdir_path / "mlva_fingerprint.tsv",
        "report": outdir_path / "report.html",
        **recruitment_paths,
        **phylogeny_paths,
    }
