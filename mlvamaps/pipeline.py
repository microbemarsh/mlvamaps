from __future__ import annotations

from pathlib import Path

from .bayesian_caller import call_loci
from .clustering import cluster_vntr_asvs
from .concurrency import DEFAULT_THREADS, resolve_threads
from .in_silico_pcr import read_amplirust_results, run_amplirust_loci
from .io import read_fastq, read_profiles, write_fasta, write_fastq, write_tsv
from .locus_assignment import assignments_from_amplirust
from .mapping import (
    MAPPING_SUMMARY_FIELDS,
    SNP_FIELDS,
    run_locus_mapping,
)
from .ml_classifier import predict_read_alleles
from .novelty import score_novelty
from .profile_matching import build_fingerprint, match_profiles
from .progress import ProgressReporter
from .primers import read_loci_or_primers
from .qc import filter_reads
from .repeat_parser import extract_repeat_features
from .report import write_report


ASSIGNMENT_FIELDS = [
    "read_id",
    "sample_id",
    "assigned_locus",
    "assignment_score",
    "orientation",
    "primer_forward_detected",
    "primer_reverse_detected",
    "passes_assignment_qc",
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
    "num_vntr_asvs",
    "dominant_vntr_asv",
    "call_status",
]

MATCH_FIELDS = [
    "sample_id",
    "best_profile_id",
    "strain_id",
    "distance",
    "matched_loci",
    "mismatched_loci",
    "confidence",
]

NOVELTY_FIELDS = ["sample_id", "nearest_profile", "novelty_score", "interpretation"]

SIMPLE_CALL_FIELDS = [
    "sample_id",
    "locus_id",
    "present",
    "repeat_count",
    "repeat_count_raw",
    "product_size_bp",
    "read_depth",
    "mean_coverage",
    "status",
    "evidence",
]


def simple_call_rows_from_alleles(sample_id: str, allele_rows: list[dict]) -> list[dict]:
    rows = []
    for row in allele_rows:
        read_depth = int(row.get("read_depth") or 0)
        present = read_depth > 0
        repeat_count = row.get("called_repeat_count", "") if present else ""
        rows.append(
            {
                "sample_id": sample_id,
                "locus_id": row["locus_id"],
                "present": "yes" if present else "no",
                "repeat_count": repeat_count,
                "repeat_count_raw": repeat_count,
                "product_size_bp": "",
                "read_depth": read_depth,
                "mean_coverage": "",
                "status": row["call_status"],
                "evidence": f"{read_depth} assigned reads" if present else "no assigned reads",
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
    min_read_length: int = 50,
    max_read_length: int = 100000,
    min_qscore: float = 0.0,
    max_primer_mismatches: int = 3,
    min_depth: int = 10,
    min_posterior: float = 0.75,
    min_cluster_size: int = 2,
    cluster_min_identity: float = 0.97,
    vsearch_bin: str = "vsearch",
    amplirust_bin: str = "amplirust",
    minibwa_bin: str = "minibwa",
    locus_mapping: bool = True,
    min_mapping_quality: int = 0,
    min_base_quality: int = 20,
    min_snp_depth: int = 3,
    min_snp_alternate_reads: int = 2,
    min_snp_frequency: float = 0.2,
    threads: int = DEFAULT_THREADS,
    show_progress: bool = False,
) -> dict[str, Path]:
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

    if filtered_reads:
        progress.step("Assigning reads by degenerate primer pairs with Amplirust")
        assignment_fasta = outdir_path / "filtered_reads.fasta"
        write_fasta(((read.read_id, read.sequence) for read in filtered_reads), assignment_fasta)
        amplirust_paths = run_amplirust_loci(
            assignment_fasta,
            loci,
            outdir_path / "amplirust",
            max_errors=max_primer_mismatches,
            threads=threads,
            executable=amplirust_bin,
        )
        assignments = assignments_from_amplirust(
            filtered_reads,
            loci,
            read_amplirust_results(amplirust_paths["stats"], amplirust_paths["products"]),
            sample_id,
            progress=progress,
        )
    else:
        assignments = []
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
        alignment_work_dir=outdir_path / "minibwa" / "cluster_memberships",
        minibwa_executable=minibwa_bin,
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

    mapping_rows: list[dict] = []
    snp_rows: list[dict] = []
    if locus_mapping:
        progress.step(
            "Mapping locus reads to dominant VSEARCH representatives with minibwa"
        )
        mapping_rows, snp_rows, _mapping_paths = run_locus_mapping(
            features,
            asv_rows,
            outdir_path,
            sample_id,
            thread_count,
            executable=minibwa_bin,
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
    predictions = predict_read_alleles(features, loci, asv_memberships)
    prediction_rows = [{field: getattr(row, field) for field in PREDICTION_FIELDS} for row in predictions]
    write_tsv(prediction_rows, outdir_path / "read_level_allele_predictions.tsv", PREDICTION_FIELDS)

    allele_rows = call_loci(predictions, loci, asv_rows, min_depth, min_posterior)
    for row in allele_rows:
        row["sample_id"] = sample_id
    write_tsv(allele_rows, outdir_path / "allele_calls.tsv", ALLELE_FIELDS)
    write_tsv(simple_call_rows_from_alleles(sample_id, allele_rows), outdir_path / "calls.tsv", SIMPLE_CALL_FIELDS)

    fingerprint_rows, probabilistic_rows = build_fingerprint(sample_id, allele_rows, loci)
    fingerprint_fields = ["sample_id"] + [locus.locus_id for locus in loci]
    write_tsv(fingerprint_rows, outdir_path / "mlva_fingerprint.tsv", fingerprint_fields)
    write_tsv(
        probabilistic_rows,
        outdir_path / "mlva_fingerprint_probabilistic.tsv",
        ["sample_id", "locus_id", "repeat_count", "posterior_probability"],
    )

    match_rows = match_profiles(sample_id, fingerprint_rows[0], profiles)
    write_tsv(match_rows, outdir_path / "profile_matches.tsv", MATCH_FIELDS)
    novelty_rows = score_novelty(sample_id, allele_rows, match_rows)
    write_tsv(novelty_rows, outdir_path / "novelty_scores.tsv", NOVELTY_FIELDS)
    progress.step("Writing HTML report")
    write_report(
        outdir_path,
        sample_id,
        allele_rows,
        novelty_rows,
        loci,
        match_rows,
        profiles,
        asv_rows,
        mapping_rows,
        snp_rows,
    )
    progress.step(f"Done. Main calls: {outdir_path / 'calls.tsv'}")

    return {
        "outdir": outdir_path,
        "calls": outdir_path / "calls.tsv",
        "allele_calls": outdir_path / "allele_calls.tsv",
        "asv_table": outdir_path / "vntr_asv_table.tsv",
        "asv_memberships": outdir_path / "vntr_asv_memberships.tsv",
        "asv_representatives": outdir_path / "vntr_asv_representatives.fasta",
        "mapping_summary": outdir_path / "locus_mapping_summary.tsv",
        "mapping_snps": outdir_path / "locus_snps.tsv",
        "mapping_references": outdir_path / "locus_mapping_references.fasta",
        "mapping_alignments": outdir_path / "locus_read_alignments.sam",
        "minibwa": outdir_path / "minibwa",
        "vsearch": vsearch_dir,
        "amplirust": outdir_path / "amplirust",
        "fingerprint": outdir_path / "mlva_fingerprint.tsv",
        "report": outdir_path / "report.html",
    }
