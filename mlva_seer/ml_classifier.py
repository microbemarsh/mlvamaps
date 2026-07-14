from __future__ import annotations

import numpy as np

from .models import Locus, ReadPrediction, RepeatFeature


def predict_read_alleles(
    features: list[RepeatFeature],
    loci: list[Locus],
    cluster_memberships: list[dict],
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
        sigma = max(0.18, 0.65 - min(feature.mean_qscore, 30) / 80)
        counts = np.arange(
            locus.expected_min_repeats,
            locus.expected_max_repeats + 1,
            dtype=np.int64,
        )
        distances = feature.raw_repeat_count_estimate - counts
        weights = np.exp(-(distances * distances) / (2 * sigma * sigma)) + 1e-9
        probabilities = weights / weights.sum(dtype=np.float64)
        ranking = np.argsort(-probabilities, kind="stable")
        best_idx = int(ranking[0])
        if len(ranking) > 1:
            alt_idx = int(ranking[1])
            alt = (int(counts[alt_idx]), float(probabilities[alt_idx]))
        else:
            alt = (None, 0.0)
        predictions.append(
            ReadPrediction(
                feature.read_id,
                feature.locus_id,
                int(counts[best_idx]),
                round(float(probabilities[best_idx]), 6),
                alt[0],
                round(alt[1], 6),
                str(membership["variant_id"]),
                insertions,
                deletions,
                substitutions,
                round(evidence_weight, 6),
            )
        )
    return predictions
