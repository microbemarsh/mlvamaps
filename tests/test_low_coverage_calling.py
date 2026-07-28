import pytest

from mlvamaps.bayesian_caller import call_loci
from mlvamaps.calling import (
    estimate_repeat_count_from_product_length,
    estimate_repeat_count_from_spanning_read,
)
from mlvamaps.ml_classifier import predict_read_alleles
from mlvamaps.mixture import estimate_variant_mixtures
from mlvamaps.models import Assignment, Locus, ReadPrediction, RepeatFeature
from mlvamaps.profile_matching import build_fingerprint, match_profiles
from mlvamaps.repeat_parser import extract_repeat_features


def _feature(raw_count: float) -> RepeatFeature:
    sequence = "ATGC" * 5
    return RepeatFeature(
        read_id="read1",
        locus_id="VNTR",
        repeat_region_start=0,
        repeat_region_end=len(sequence),
        repeat_region_length_bp=len(sequence),
        repeat_motif="ATGC",
        raw_repeat_count_estimate=raw_count,
        nearest_integer_repeat_count=round(raw_count),
        flank_quality_score=1.0,
        repeat_pattern="ATGC-ATGC-ATGC-ATGC-ATGC",
        repeat_sequence=sequence,
        mean_qscore=20.0,
        mismatch_count_in_repeat_region=0,
        motif_kmer_count=5,
        left_primer_score=1.0,
        right_primer_score=1.0,
        left_flank_score=1.0,
        right_flank_score=1.0,
        amplicon_sequence=sequence,
        amplicon_quality="5" * len(sequence),
        product_size_bp=len(sequence),
        repeat_measurement_method="assembly_product_length",
    )


def _membership() -> list[dict]:
    return [
        {
            "read_id": "read1",
            "locus_id": "VNTR",
            "variant_id": "VNTR_ASV1",
            "insertions_vs_representative": 0,
            "deletions_vs_representative": 0,
            "substitutions_vs_representative": 0,
            "aligned_repeat_sequence": "ATGC" * 5,
            "aligned_representative_sequence": "ATGC" * 5,
        }
    ]


def test_spanning_read_uses_same_product_calibration_as_assembly():
    locus = Locus(
        locus_id="VNTR",
        forward_primer="AAAA",
        reverse_primer="TTTT",
        repeat_unit_length_bp=4,
        expected_product_size_bp=40,
        nominal_repeat_units=5,
    )

    read_count, method = estimate_repeat_count_from_spanning_read(
        locus,
        product_size_bp=44,
        repeat_region_size_bp=20,
        flanks_resolved=True,
    )

    assert read_count == estimate_repeat_count_from_product_length(locus, 44) == 6
    assert method == "assembly_product_length"


def test_resolved_flanks_do_not_subtract_nonrepeat_sequence_twice():
    locus = Locus(locus_id="VNTR", repeat_unit_length_bp=4)

    count, method = estimate_repeat_count_from_spanning_read(
        locus,
        product_size_bp=60,
        repeat_region_size_bp=20,
        flanks_resolved=True,
    )

    assert count == 5
    assert method == "flank_bounded_repeat_length"


def test_read_measurement_uses_primer_match_product_size():
    locus = Locus(
        locus_id="VNTR",
        forward_primer="AAAA",
        reverse_primer="TTTT",
        repeat_motif="ATGC",
        repeat_unit_length_bp=4,
        expected_product_size_bp=40,
        nominal_repeat_units=5,
    )
    sequence = "AAAA" + ("ATGC" * 5) + "TTTT"
    assignment = Assignment(
        read_id="read1",
        sample_id="sample",
        assigned_locus="VNTR",
        assignment_score=1.0,
        orientation="forward",
        primer_forward_detected=True,
        primer_reverse_detected=True,
        passes_assignment_qc=True,
        oriented_sequence=sequence,
        oriented_quality="5" * len(sequence),
        forward_start=0,
        forward_end=4,
        reverse_start=len(sequence) - 4,
        reverse_end=len(sequence),
        forward_mismatches=0,
        reverse_mismatches=0,
        product_size_bp=44,
    )

    feature = extract_repeat_features([assignment], [locus], threads=1)[0]

    assert feature.product_size_bp == 44
    assert feature.raw_repeat_count_estimate == 6


def test_fastq_prediction_uses_assembly_rounding_but_retains_raw_value():
    locus = Locus(
        locus_id="VNTR",
        repeat_unit_length_bp=4,
        expected_min_repeats=3,
        expected_max_repeats=7,
    )

    prediction = predict_read_alleles(
        [_feature(5.25)],
        [locus],
        _membership(),
        assembly_equivalent=True,
    )[0]

    assert prediction.raw_repeat_count_estimate == 5.25
    assert prediction.measurement_repeat_count_estimate == 5.5
    assert prediction.predicted_repeat_count == 5.5


def test_single_spanning_read_remains_in_signature_with_low_depth_status():
    locus = Locus(
        locus_id="VNTR",
        repeat_unit_length_bp=4,
        expected_min_repeats=3,
        expected_max_repeats=7,
    )
    prediction = predict_read_alleles([_feature(5.0)], [locus], _membership())[0]

    row = call_loci([prediction], [locus], [], min_depth=3)[0]
    fingerprint, _probabilistic = build_fingerprint("sample", [row], [locus])

    assert row["called_repeat_count"] == 5
    assert row["call_status"] == "LOW_DEPTH"
    assert fingerprint[0]["VNTR"] == 5


def test_single_spanning_read_is_retained_by_default():
    locus = Locus(
        locus_id="VNTR",
        repeat_unit_length_bp=4,
        expected_min_repeats=3,
        expected_max_repeats=7,
    )
    prediction = predict_read_alleles([_feature(5.0)], [locus], _membership())[0]

    row = call_loci([prediction], [locus], [])[0]

    assert row["called_repeat_count"] == 5
    assert row["evidence_status"] == "SINGLE_MOLECULE_PROVISIONAL"
    assert row["call_status"] != "LOW_DEPTH"


def test_minimum_depth_cannot_be_below_one():
    with pytest.raises(ValueError, match="at least 1"):
        call_loci([], [], [], min_depth=0)


def test_metagenome_mode_flags_any_meaningful_secondary_allele():
    locus = Locus(locus_id="VNTR", expected_min_repeats=3, expected_max_repeats=7)
    predictions = [
        ReadPrediction(
            f"major{index}",
            "VNTR",
            5,
            1.0,
            5.5,
            0.0,
            "VNTR_ASV1",
            0,
            0,
            0,
            1.0,
            5.0,
            0.08,
            5.0,
        )
        for index in range(9)
    ] + [
        ReadPrediction(
            "minor",
            "VNTR",
            6,
            1.0,
            5.5,
            0.0,
            "VNTR_ASV2",
            0,
            0,
            0,
            1.0,
            6.0,
            0.08,
            6.0,
        )
    ]
    asvs = [
        {"locus_id": "VNTR", "variant_id": "VNTR_ASV1", "support_reads": 9, "frequency": 0.9},
        {"locus_id": "VNTR", "variant_id": "VNTR_ASV2", "support_reads": 1, "frequency": 0.1},
    ]
    mixtures = [
        {
            "locus_id": "VNTR",
            "variant_id": "VNTR_ASV1",
            "estimated_fraction": 0.9,
            "meaningful": "yes",
        },
        {
            "locus_id": "VNTR",
            "variant_id": "VNTR_ASV2",
            "estimated_fraction": 0.1,
            "meaningful": "yes",
        },
    ]

    isolate = call_loci(
        predictions,
        [locus],
        asvs,
        min_depth=1,
        mixture_rows=mixtures,
        sample_mode="isolate",
    )[0]
    metagenome = call_loci(
        predictions,
        [locus],
        asvs,
        min_depth=1,
        mixture_rows=mixtures,
    )[0]

    assert isolate["call_status"] == "PASS"
    assert metagenome["call_status"] == "MULTIPLE_VARIANTS"
    assert metagenome["called_repeat_count"] == 5


def test_repeated_primary_reads_increase_confidence_with_a_cap():
    locus = Locus(locus_id="VNTR", expected_min_repeats=3, expected_max_repeats=7)

    def predictions(count: int) -> list[ReadPrediction]:
        return [
            ReadPrediction(
                f"read{index}",
                "VNTR",
                5,
                0.6,
                5.5,
                0.4,
                "VNTR_ASV1",
                0,
                0,
                0,
                1.0,
                5.2,
                0.3,
                5.2,
            )
            for index in range(count)
        ]

    one = call_loci(predictions(1), [locus], [], min_depth=1)[0]
    five = call_loci(predictions(5), [locus], [], min_depth=1)[0]
    capped = call_loci(
        predictions(50),
        [locus],
        [],
        min_depth=1,
        max_confidence_depth=5,
    )[0]

    assert five["posterior_probability"] > one["posterior_probability"]
    assert capped["posterior_probability"] == five["posterior_probability"]
    assert capped["confidence_effective_depth"] == 5.0


def test_secondary_variant_does_not_muddy_primary_allele_posterior():
    locus = Locus(locus_id="VNTR", expected_min_repeats=3, expected_max_repeats=7)
    primary = [
        ReadPrediction(
            f"primary{index}",
            "VNTR",
            5,
            0.8,
            5.5,
            0.2,
            "VNTR_ASV1",
            0,
            0,
            0,
            1.0,
            5.1,
            0.25,
            5.1,
        )
        for index in range(8)
    ]
    secondary = [
        ReadPrediction(
            f"secondary{index}",
            "VNTR",
            6,
            0.8,
            5.5,
            0.2,
            "VNTR_ASV2",
            0,
            0,
            0,
            1.0,
            6.0,
            0.25,
            6.0,
        )
        for index in range(2)
    ]
    mixtures = [
        {
            "locus_id": "VNTR",
            "variant_id": "VNTR_ASV1",
            "estimated_fraction": 0.8,
            "evidence_class": "DOMINANT",
            "meaningful": "yes",
        },
        {
            "locus_id": "VNTR",
            "variant_id": "VNTR_ASV2",
            "estimated_fraction": 0.2,
            "evidence_class": "CONFIRMED_SECONDARY",
            "meaningful": "yes",
        },
    ]
    asvs = [
        {"locus_id": "VNTR", "variant_id": "VNTR_ASV1", "support_reads": 8, "frequency": 0.8},
        {"locus_id": "VNTR", "variant_id": "VNTR_ASV2", "support_reads": 2, "frequency": 0.2},
    ]

    primary_only = call_loci(primary, [locus], asvs[:1], min_depth=1)[0]
    mixed = call_loci(
        primary + secondary,
        [locus],
        asvs,
        min_depth=1,
        mixture_rows=mixtures,
    )[0]

    assert mixed["called_repeat_count"] == primary_only["called_repeat_count"] == 5
    assert mixed["posterior_probability"] == primary_only["posterior_probability"]
    assert mixed["primary_read_depth"] == 8
    assert mixed["read_depth"] == 10
    assert mixed["call_status"] == "MULTIPLE_VARIANTS"
    assert "VNTR_ASV2|6|0.200000|CONFIRMED_SECONDARY" in mixed["secondary_alleles"]


def test_singleton_secondary_is_candidate_not_signature_changing():
    asvs = [
        {
            "sample_id": "sample",
            "locus_id": "VNTR",
            "variant_id": "VNTR_ASV1",
            "repeat_count": 5,
            "support_reads": 9,
            "representative_sequence": "ATGC" * 5,
            "representative_length_bp": 20,
            "total_insertions": 0,
            "total_deletions": 0,
            "total_substitutions": 0,
        },
        {
            "sample_id": "sample",
            "locus_id": "VNTR",
            "variant_id": "VNTR_ASV2",
            "repeat_count": 6,
            "support_reads": 1,
            "representative_sequence": "TTTT" * 6,
            "representative_length_bp": 24,
            "total_insertions": 0,
            "total_deletions": 0,
            "total_substitutions": 0,
        },
    ]

    rows = estimate_variant_mixtures(
        asvs,
        min_fraction=0.01,
        min_secondary_reads=2,
    )
    by_variant = {row["variant_id"]: row for row in rows}

    assert by_variant["VNTR_ASV2"]["evidence_class"] == "CANDIDATE"
    assert by_variant["VNTR_ASV2"]["meaningful"] == "no"
    assert by_variant["VNTR_ASV2"]["abundance_class"] == "SECONDARY"


def test_partial_profile_ranking_uses_full_allele_probabilities():
    fingerprint = {"sample_id": "sample", "L1": 5, "L2": 4}
    profiles = [
        {"profile_id": "A", "strain_id": "A", "L1": "5", "L2": "4.5"},
        {"profile_id": "B", "strain_id": "B", "L1": "5.5", "L2": "4"},
    ]
    allele_rows = [
        {"locus_id": "L1", "allele_distribution": "5:0.51;5.5:0.49"},
        {"locus_id": "L2", "allele_distribution": "4:0.9;4.5:0.1"},
    ]

    matches = match_profiles(
        "sample", fingerprint, profiles, allele_rows=allele_rows
    )

    assert matches[0]["best_profile_id"] == "B"
    assert matches[0]["compared_loci"] == 2


def test_profile_distance_has_same_priority_with_or_without_read_probabilities():
    fingerprint = {"sample_id": "sample", "L1": 5}
    profiles = [
        {"profile_id": "ASSEMBLY_MATCH", "strain_id": "A", "L1": "5"},
        {"profile_id": "LIKELIHOOD_ONLY", "strain_id": "B", "L1": "5.5"},
    ]
    allele_rows = [
        {"locus_id": "L1", "allele_distribution": "5:0.1;5.5:0.9"}
    ]

    assembly_matches = match_profiles("sample", fingerprint, profiles)
    fastq_matches = match_profiles(
        "sample", fingerprint, profiles, allele_rows=allele_rows
    )

    assert assembly_matches[0]["best_profile_id"] == "ASSEMBLY_MATCH"
    assert fastq_matches[0]["best_profile_id"] == "ASSEMBLY_MATCH"
