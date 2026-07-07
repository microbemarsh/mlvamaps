from __future__ import annotations

import math

from .models import Locus, ReadPrediction, RepeatFeature


def predict_read_alleles(features: list[RepeatFeature], loci: list[Locus]) -> list[ReadPrediction]:
    by_locus = {locus.locus_id: locus for locus in loci}
    predictions = []
    for feature in features:
        locus = by_locus[feature.locus_id]
        sigma = max(0.18, 0.65 - min(feature.mean_qscore, 30) / 80)
        weights = {}
        for count in range(locus.expected_min_repeats, locus.expected_max_repeats + 1):
            distance = feature.raw_repeat_count_estimate - count
            likelihood = math.exp(-(distance * distance) / (2 * sigma * sigma))
            quality_weight = 0.5 + min(feature.flank_quality_score, 1.0) / 2
            indel_penalty = 1 / (1 + feature.indel_count_in_repeat_region)
            weights[count] = likelihood * quality_weight * indel_penalty + 1e-9
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
            )
        )
    return predictions
