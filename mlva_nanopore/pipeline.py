from __future__ import annotations

from pathlib import Path

from .bayesian_caller import call_loci
from .clustering import cluster_vntr_asvs
from .io import read_fastq, read_profiles, write_fasta, write_fastq, write_tsv
from .locus_assignment import assign_reads
from .ml_classifier import predict_read_alleles
from .novelty import score_novelty
from .profile_matching import build_fingerprint, match_profiles
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
    "mean_qscore",
    "indel_count_in_repeat_region",
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
    "frequency",
    "consensus_pattern",
]

PREDICTION_FIELDS = [
    "read_id",
    "locus_id",
    "predicted_repeat_count",
    "probability",
    "top_alt_repeat_count",
    "top_alt_probability",
]

ALLELE_FIELDS = [
    "sample_id",
    "locus_id",
    "called_repeat_count",
    "posterior_probability",
    "second_best_repeat_count",
    "second_best_posterior",
    "read_depth",
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
    threads: int = 0,
) -> dict[str, Path]:
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)

    loci = read_loci_or_primers(loci_path, primers_path)
    profiles = read_profiles(profiles_path)
    reads = list(read_fastq(reads_path))
    filtered_reads, qc_rows = filter_reads(reads, min_read_length, max_read_length, min_qscore)
    write_tsv(qc_rows, outdir_path / "qc_summary.tsv", ["metric", "value"])
    write_fastq(filtered_reads, outdir_path / "filtered_reads.fastq.gz")

    assignments = assign_reads(filtered_reads, loci, sample_id, max_primer_mismatches, threads=threads)
    assignment_rows = [{field: getattr(row, field) for field in ASSIGNMENT_FIELDS} for row in assignments]
    write_tsv(assignment_rows, outdir_path / "read_locus_assignments.tsv", ASSIGNMENT_FIELDS)

    features = extract_repeat_features(assignments, loci, threads=threads)
    feature_rows = [{field: getattr(row, field) for field in FEATURE_FIELDS} for row in features]
    write_tsv(feature_rows, outdir_path / "read_repeat_features.tsv", FEATURE_FIELDS)

    asv_rows, fasta_records = cluster_vntr_asvs(features)
    for row in asv_rows:
        row["sample_id"] = sample_id
    write_tsv(asv_rows, outdir_path / "vntr_asv_table.tsv", ASV_FIELDS)
    write_fasta(fasta_records, outdir_path / "vntr_asv_consensus.fasta")

    predictions = predict_read_alleles(features, loci)
    prediction_rows = [{field: getattr(row, field) for field in PREDICTION_FIELDS} for row in predictions]
    write_tsv(prediction_rows, outdir_path / "read_level_allele_predictions.tsv", PREDICTION_FIELDS)

    allele_rows = call_loci(predictions, loci, asv_rows, min_depth, min_posterior)
    for row in allele_rows:
        row["sample_id"] = sample_id
    write_tsv(allele_rows, outdir_path / "allele_calls.tsv", ALLELE_FIELDS)

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
    write_report(outdir_path, sample_id, allele_rows, novelty_rows, loci, match_rows, profiles, asv_rows)

    return {
        "outdir": outdir_path,
        "allele_calls": outdir_path / "allele_calls.tsv",
        "fingerprint": outdir_path / "mlva_fingerprint.tsv",
        "report": outdir_path / "report.html",
    }
