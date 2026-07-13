from __future__ import annotations

import math

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
        membership = membership_by_read[(feature.locus_id, feature.read_id)]
        insertions = int(membership["insertions_vs_consensus"])
        deletions = int(membership["deletions_vs_consensus"])
        substitutions = int(membership["substitutions_vs_consensus"])
        aligned_length = max(
            len(str(membership["aligned_repeat_sequence"])),
            len(str(membership["aligned_consensus_sequence"])),
            1,
        )
        alignment_identity = max(0.0, 1 - ((insertions + deletions + substitutions) / aligned_length))
        quality_weight = 0.5 + max(0.0, min(feature.mean_qscore, 30)) / 60
        flank_weight = max(0.0, min(feature.flank_quality_score, 1.0))
        evidence_weight = flank_weight * quality_weight * alignment_identity
        sigma = max(0.18, 0.65 - min(feature.mean_qscore, 30) / 80)
        weights = {}
        for count in range(locus.expected_min_repeats, locus.expected_max_repeats + 1):
            distance = feature.raw_repeat_count_estimate - count
            likelihood = math.exp(-(distance * distance) / (2 * sigma * sigma))
            weights[count] = likelihood + 1e-9
        total = sum(weights.values())
        ranked = sorted(((count, weight / total) for count, weight in weights.items()), key=lambda item: item[1], reverse=True)
        alt = ranked[1] if len(ranked) > 1 else (None, 0.0)
        predictions.append(
            ReadPrediction(
                feature.read_id,
                feature.locus_id,
                ranked[0][0],
                round(ranked[0][1], 6),
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
