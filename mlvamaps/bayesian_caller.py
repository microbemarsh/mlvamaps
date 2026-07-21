from __future__ import annotations

from collections import defaultdict

from .calling import allele_grid, gaussian_allele_probabilities
from .models import Locus, ReadPrediction


def call_loci(
    predictions: list[ReadPrediction],
    loci: list[Locus],
    asv_rows: list[dict],
    min_depth: int = 10,
    min_posterior: float = 0.75,
    mixture_rows: list[dict] | None = None,
) -> list[dict]:
    pred_by_locus: dict[str, list[ReadPrediction]] = defaultdict(list)
    for prediction in predictions:
        pred_by_locus[prediction.locus_id].append(prediction)

    asv_by_locus: dict[str, list[dict]] = defaultdict(list)
    for row in asv_rows:
        asv_by_locus[row["locus_id"]].append(row)

    mixture_by_locus: dict[str, list[dict]] = defaultdict(list)
    for row in mixture_rows or []:
        mixture_by_locus[row["locus_id"]].append(row)

    rows = []
    for locus in loci:
        preds = pred_by_locus.get(locus.locus_id, [])
        if not preds:
            rows.append(
                {
                    "sample_id": "",
                    "locus_id": locus.locus_id,
                    "called_repeat_count": "",
                    "posterior_probability": 0.0,
                    "second_best_repeat_count": "",
                    "second_best_posterior": 0.0,
                    "read_depth": 0,
                    "effective_read_depth": 0.0,
                    "num_vntr_asvs": 0,
                    "num_meaningful_variants": 0,
                    "dominant_vntr_asv": "",
                    "dominant_variant_fraction": 0.0,
                    "call_status": "LOCUS_DROPOUT",
                }
            )
            continue

        candidates = allele_grid(locus, step=0.5)
        weights = {count: 1e-6 for count in candidates}
        effective_read_depth = sum(pred.evidence_weight for pred in preds)
        measurement_groups: dict[tuple[float, float], float] = defaultdict(float)
        for pred in preds:
            if pred.raw_repeat_count_estimate is not None and pred.measurement_sigma is not None:
                measurement_groups[
                    (pred.raw_repeat_count_estimate, pred.measurement_sigma)
                ] += pred.evidence_weight
            else:
                weights[pred.predicted_repeat_count] = (
                    weights.get(pred.predicted_repeat_count, 1e-6)
                    + pred.probability * pred.evidence_weight
                )
                if pred.top_alt_repeat_count is not None:
                    weights[pred.top_alt_repeat_count] = (
                        weights.get(pred.top_alt_repeat_count, 1e-6)
                        + pred.top_alt_probability * pred.evidence_weight
                    )
        # Identical read lengths and quality-derived sigmas are common at high
        # depth. Collapse them before evaluating the grid to keep inference
        # proportional to distinct observations rather than total reads.
        for (raw_count, sigma), group_weight in measurement_groups.items():
            probabilities = gaussian_allele_probabilities(
                raw_count, candidates, sigma
            )
            for count, probability in zip(candidates, probabilities):
                weights[count] += probability * group_weight
        total = sum(weights.values())
        ranked = sorted(((count, weight / total) for count, weight in weights.items()), key=lambda item: item[1], reverse=True)
        best = ranked[0]
        second = ranked[1] if len(ranked) > 1 else ("", 0.0)
        locus_asvs = sorted(asv_by_locus.get(locus.locus_id, []), key=lambda row: row["support_reads"], reverse=True)
        locus_mixture = sorted(
            mixture_by_locus.get(locus.locus_id, []),
            key=lambda row: float(row.get("estimated_fraction") or 0),
            reverse=True,
        )
        meaningful_variants = [
            row for row in locus_mixture if str(row.get("meaningful", "")).lower() == "yes"
        ]
        if locus_mixture:
            dominant = str(locus_mixture[0]["variant_id"])
            dominant_freq = float(locus_mixture[0].get("estimated_fraction") or 0)
            meaningful_count = len(meaningful_variants)
        else:
            dominant = locus_asvs[0]["variant_id"] if locus_asvs else ""
            dominant_freq = float(locus_asvs[0]["frequency"]) if locus_asvs else 0.0
            meaningful_count = len(locus_asvs)
        status = "PASS"
        if len(preds) < min_depth:
            status = "LOW_DEPTH"
        elif best[1] < min_posterior or (best[1] - second[1]) < 0.2:
            status = "AMBIGUOUS"
        elif best[0] < locus.expected_min_repeats or best[0] > locus.expected_max_repeats:
            status = "OUT_OF_RANGE"
        elif meaningful_count > 1 and dominant_freq < 0.8:
            status = "MULTIPLE_VARIANTS"
        rows.append(
            {
                "sample_id": "",
                "locus_id": locus.locus_id,
                "called_repeat_count": best[0],
                "posterior_probability": round(best[1], 6),
                "second_best_repeat_count": second[0],
                "second_best_posterior": round(second[1], 6),
                "read_depth": len(preds),
                "effective_read_depth": round(effective_read_depth, 4),
                "num_vntr_asvs": len(locus_asvs),
                "num_meaningful_variants": meaningful_count,
                "dominant_vntr_asv": dominant,
                "dominant_variant_fraction": round(dominant_freq, 6),
                "allele_distribution": ";".join(
                    f"{count}:{probability:.6f}" for count, probability in ranked
                ),
                "call_status": status,
            }
        )
    return rows
