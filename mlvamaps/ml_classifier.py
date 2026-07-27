from __future__ import annotations

import numpy as np

from .calling import (
    allele_grid,
    legacy_round_repeat_count,
    repeat_unit_length,
)
from .models import Locus, ReadPrediction, RepeatFeature


def predict_read_alleles(
    features: list[RepeatFeature],
    loci: list[Locus],
    cluster_memberships: list[dict],
    assembly_equivalent: bool = True,
    assembly_round_tolerance: float = 0.25,
) -> list[ReadPrediction]:
    by_locus = {locus.locus_id: locus for locus in loci}
    membership_by_read = {
        (row["locus_id"], row["read_id"]): row
        for row in cluster_memberships
    }
    predictions = []
    for feature in features:
        locus = by_locus[feature.locus_id]
        membership = membership_by_read.get((feature.locus_id, feature.read_id))
        if membership is None:
            # Clustering deliberately removes reads that do not support a retained variant.
            continue
        insertions = int(membership["insertions_vs_representative"])
        deletions = int(membership["deletions_vs_representative"])
        substitutions = int(membership["substitutions_vs_representative"])
        aligned_length = max(
            len(str(membership["aligned_repeat_sequence"])),
            len(str(membership["aligned_representative_sequence"])),
            1,
        )
        alignment_identity = max(0.0, 1 - ((insertions + deletions + substitutions) / aligned_length))
        quality_weight = 0.5 + max(0.0, min(feature.mean_qscore, 30)) / 60
        flank_weight = max(0.0, min(feature.flank_quality_score, 1.0))
        evidence_weight = flank_weight * quality_weight * alignment_identity
        # Q17 corresponds to roughly 98% per-base accuracy. At that quality,
        # primer-spanning reads should use a sharp assembly-like measurement;
        # lower qualities broaden the distribution without moving its center.
        error_probability = 10 ** (-max(feature.mean_qscore, 0.0) / 10)
        sigma = max(
            0.08,
            0.5 / max(repeat_unit_length(locus), 1),
            0.12 + (2.0 * error_probability),
        )
        counts = allele_grid(locus, step=0.5)
        measurement_center = (
            legacy_round_repeat_count(
                feature.raw_repeat_count_estimate,
                assembly_round_tolerance,
            )
            if assembly_equivalent
            else feature.raw_repeat_count_estimate
        )
        count_values = np.asarray(counts, dtype=np.float64)
        distances = float(measurement_center) - count_values
        squared_distances = distances * distances
        weights = np.exp(
            -(squared_distances - squared_distances.min()) / (2 * sigma * sigma)
        )
        probabilities = weights / weights.sum(dtype=np.float64)
        ranking = np.argsort(-probabilities, kind="stable")
        best_idx = int(ranking[0])
        if len(ranking) > 1:
            alt_idx = int(ranking[1])
            alt = (counts[alt_idx], float(probabilities[alt_idx]))
        else:
            alt = (None, 0.0)
        predictions.append(
            ReadPrediction(
                feature.read_id,
                feature.locus_id,
                counts[best_idx],
                round(float(probabilities[best_idx]), 6),
                alt[0],
                round(alt[1], 6),
                str(membership["variant_id"]),
                insertions,
                deletions,
                substitutions,
                round(evidence_weight, 6),
                feature.raw_repeat_count_estimate,
                round(sigma, 6),
                float(measurement_center),
            )
        )
    return predictions
