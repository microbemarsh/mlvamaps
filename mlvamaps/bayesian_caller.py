from __future__ import annotations

import math
from collections import defaultdict
from statistics import median

from .calling import allele_grid, gaussian_allele_probabilities
from .models import Locus, ReadPrediction


def _prediction_probabilities(
    prediction: ReadPrediction,
    candidates: list[int | float],
    measurement_override: float | None = None,
) -> list[float]:
    measurement = (
        measurement_override
        if measurement_override is not None
        else (
            prediction.measurement_repeat_count_estimate
            if prediction.measurement_repeat_count_estimate is not None
            else prediction.raw_repeat_count_estimate
        )
    )
    if measurement is not None and prediction.measurement_sigma is not None:
        return gaussian_allele_probabilities(
            measurement,
            candidates,
            prediction.measurement_sigma,
        )
    probabilities = {candidate: 1e-12 for candidate in candidates}
    probabilities[prediction.predicted_repeat_count] = max(
        prediction.probability, 1e-12
    )
    if prediction.top_alt_repeat_count is not None:
        probabilities[prediction.top_alt_repeat_count] = max(
            prediction.top_alt_probability, 1e-12
        )
    total = sum(probabilities.values())
    return [probabilities[candidate] / total for candidate in candidates]


def _combined_allele_posterior(
    predictions: list[ReadPrediction],
    candidates: list[int | float],
    max_confidence_depth: float,
    measurement_override: float | None = None,
) -> tuple[list[tuple[int | float, float]], float, float]:
    """Combine independent read evidence while capping correlated confidence."""
    raw_effective_depth = sum(
        max(0.0, prediction.evidence_weight) for prediction in predictions
    )
    capped_effective_depth = min(raw_effective_depth, max_confidence_depth)
    scale = (
        capped_effective_depth / raw_effective_depth
        if raw_effective_depth > 0
        else 0.0
    )
    log_likelihoods = [0.0 for _candidate in candidates]
    for prediction in predictions:
        probabilities = _prediction_probabilities(
            prediction,
            candidates,
            measurement_override=measurement_override,
        )
        evidence_weight = max(0.0, prediction.evidence_weight) * scale
        for index, probability in enumerate(probabilities):
            log_likelihoods[index] += evidence_weight * math.log(
                max(probability, 1e-300)
            )
    maximum = max(log_likelihoods, default=0.0)
    weights = [math.exp(value - maximum) for value in log_likelihoods]
    total = sum(weights)
    if total <= 0:
        normalized = [1.0 / len(candidates) for _candidate in candidates]
    else:
        normalized = [weight / total for weight in weights]
    ranked = sorted(
        zip(candidates, normalized),
        key=lambda item: (-item[1], float(item[0])),
    )
    return ranked, raw_effective_depth, capped_effective_depth


def _variant_allele(predictions: list[ReadPrediction]) -> int | float | str:
    support: dict[int | float, float] = defaultdict(float)
    for prediction in predictions:
        support[prediction.predicted_repeat_count] += max(
            0.0, prediction.evidence_weight
        )
    if not support:
        return ""
    return min(support, key=lambda allele: (-support[allele], float(allele)))


def call_loci(
    predictions: list[ReadPrediction],
    loci: list[Locus],
    asv_rows: list[dict],
    min_depth: int = 10,
    min_posterior: float = 0.75,
    mixture_rows: list[dict] | None = None,
    sample_mode: str = "metagenome",
    calling_convention: str = "probabilistic",
    max_confidence_depth: float = 25.0,
    primary_product_measurements: dict[str, dict] | None = None,
    read_evidence_rows: list[dict] | None = None,
) -> list[dict]:
    if sample_mode not in {"isolate", "metagenome"}:
        raise ValueError("sample_mode must be 'isolate' or 'metagenome'")
    if calling_convention not in {"assembly", "probabilistic"}:
        raise ValueError("calling_convention must be 'assembly' or 'probabilistic'")
    if max_confidence_depth <= 0:
        raise ValueError("max_confidence_depth must be positive")
    primary_product_measurements = primary_product_measurements or {}
    evidence_by_locus: dict[str, list[dict]] = defaultdict(list)
    for evidence in read_evidence_rows or []:
        evidence_by_locus[str(evidence.get("locus_id", ""))].append(evidence)
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
                    "primary_read_depth": 0,
                    "primary_effective_read_depth": 0.0,
                    "confidence_effective_depth": 0.0,
                    "num_vntr_asvs": 0,
                    "num_meaningful_variants": 0,
                    "num_candidate_variants": 0,
                    "num_confirmed_secondary_variants": 0,
                    "dominant_vntr_asv": "",
                    "dominant_variant_fraction": 0.0,
                    "secondary_alleles": "",
                    "call_status": "LOCUS_DROPOUT",
                    "sample_mode": sample_mode,
                    "calling_convention": calling_convention,
                    "primary_product_size_bp": "",
                    "primary_repeat_count_raw": "",
                    "primary_measurement_source": "",
                    "evidence_status": "NO_INFORMATIVE_READS",
                    "credible_alleles": "",
                    "full_product_reads": 0,
                    "repeat_informative_reads": 0,
                    "forward_strand_reads": 0,
                    "reverse_strand_reads": 0,
                    "median_mapping_quality": "",
                    "median_anchor_identity": "",
                    "consensus_read_agreement": "",
                }
            )
            continue

        locus_asvs = sorted(asv_by_locus.get(locus.locus_id, []), key=lambda row: row["support_reads"], reverse=True)
        locus_mixture = sorted(
            mixture_by_locus.get(locus.locus_id, []),
            key=lambda row: float(row.get("estimated_fraction") or 0),
            reverse=True,
        )
        meaningful_variants = [
            row for row in locus_mixture if str(row.get("meaningful", "")).lower() == "yes"
        ]
        candidate_variants = [
            row
            for row in locus_mixture
            if str(row.get("evidence_class", "")).upper() == "CANDIDATE"
        ]
        confirmed_secondary = [
            row
            for row in locus_mixture
            if str(row.get("evidence_class", "")).upper()
            == "CONFIRMED_SECONDARY"
        ]
        if locus_mixture:
            dominant = str(locus_mixture[0]["variant_id"])
            dominant_freq = float(locus_mixture[0].get("estimated_fraction") or 0)
            meaningful_count = len(meaningful_variants)
        else:
            dominant = locus_asvs[0]["variant_id"] if locus_asvs else ""
            dominant_freq = float(locus_asvs[0]["frequency"]) if locus_asvs else 0.0
            meaningful_count = 1 if dominant else 0
        primary_preds = [
            prediction for prediction in preds if prediction.variant_id == dominant
        ]
        if not primary_preds:
            primary_preds = preds
        candidates = allele_grid(locus, step=0.5)
        primary_measurement = primary_product_measurements.get(
            locus.locus_id, {}
        )
        measurement_override = primary_measurement.get("called_repeat_count")
        ranked, primary_effective_depth, confidence_effective_depth = (
            _combined_allele_posterior(
                primary_preds,
                candidates,
                max_confidence_depth,
                measurement_override=(
                    float(measurement_override)
                    if measurement_override not in ("", None)
                    else None
                ),
            )
        )
        best = ranked[0]
        second = ranked[1] if len(ranked) > 1 else ("", 0.0)
        effective_read_depth = sum(pred.evidence_weight for pred in preds)
        predictions_by_variant: dict[str, list[ReadPrediction]] = defaultdict(list)
        for prediction in preds:
            predictions_by_variant[prediction.variant_id].append(prediction)
        secondary_alleles = []
        for variant in locus_mixture[1:]:
            evidence_class = str(
                variant.get("evidence_class")
                or (
                    "CONFIRMED_SECONDARY"
                    if str(variant.get("meaningful", "")).lower() == "yes"
                    else "TRACE"
                )
            )
            if evidence_class == "TRACE":
                continue
            variant_id = str(variant.get("variant_id", ""))
            allele = _variant_allele(predictions_by_variant.get(variant_id, []))
            secondary_alleles.append(
                f"{variant_id}|{allele}|"
                f"{float(variant.get('estimated_fraction') or 0):.6f}|"
                f"{evidence_class}"
            )
        status = "PASS"
        locus_evidence = evidence_by_locus.get(locus.locus_id, [])
        full_product_reads = sum(row.get("measurement_status", row.get("evidence_class")) == "FULL_PRODUCT" for row in locus_evidence)
        repeat_informative_reads = sum(row.get("measurement_status", row.get("evidence_class")) == "REPEAT_INFORMATIVE" for row in locus_evidence)
        strand_counts = {
            "forward": sum(row.get("strand") == "forward" for row in locus_evidence),
            "reverse": sum(row.get("strand") == "reverse" for row in locus_evidence),
        }
        mapping_qualities = [float(row["mapping_quality"]) for row in locus_evidence if row.get("mapping_quality") not in ("", None)]
        anchor_identities = [
            (float(row["forward_anchor_identity"]) + float(row["reverse_anchor_identity"])) / 2
            for row in locus_evidence
            if row.get("forward_anchor_identity") not in ("", None)
            and row.get("reverse_anchor_identity") not in ("", None)
        ]
        credible, credible_mass = [], 0.0
        for allele, probability in ranked:
            credible.append(allele)
            credible_mass += probability
            if credible_mass >= 0.95:
                break
        if len(primary_preds) == 1:
            evidence_status = "SINGLE_MOLECULE_PROVISIONAL"
        elif len(primary_preds) < min_depth:
            evidence_status = "PROVISIONAL_LOW_DEPTH"
        elif best[1] < min_posterior or (best[1] - second[1]) < 0.2:
            evidence_status = "AMBIGUOUS"
        else:
            evidence_status = "CONFIDENT"
        if len(primary_preds) < min_depth:
            status = "LOW_DEPTH"
        elif best[1] < min_posterior or (best[1] - second[1]) < 0.2:
            status = "AMBIGUOUS"
        elif best[0] < locus.expected_min_repeats or best[0] > locus.expected_max_repeats:
            status = "OUT_OF_RANGE"
        elif meaningful_count > 1 and (
            sample_mode == "metagenome" or dominant_freq < 0.8
        ):
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
                "primary_read_depth": len(primary_preds),
                "primary_effective_read_depth": round(primary_effective_depth, 4),
                "confidence_effective_depth": round(confidence_effective_depth, 4),
                "num_vntr_asvs": len(locus_asvs),
                "num_meaningful_variants": meaningful_count,
                "num_candidate_variants": len(candidate_variants),
                "num_confirmed_secondary_variants": len(confirmed_secondary),
                "dominant_vntr_asv": dominant,
                "dominant_variant_fraction": round(dominant_freq, 6),
                "secondary_alleles": ";".join(secondary_alleles),
                "allele_distribution": ";".join(
                    f"{count}:{probability:.6f}" for count, probability in ranked
                ),
                "call_status": status,
                "sample_mode": sample_mode,
                "calling_convention": calling_convention,
                "primary_product_size_bp": primary_measurement.get(
                    "product_size_bp", ""
                ),
                "primary_repeat_count_raw": primary_measurement.get(
                    "raw_repeat_count", ""
                ),
                "primary_measurement_source": primary_measurement.get(
                    "source",
                    (
                        "recruitment_partial"
                        if dominant.endswith("_RECRUITED")
                        else "per_read_fallback"
                        if calling_convention == "assembly"
                        else "read_distribution"
                    ),
                ),
                "evidence_status": evidence_status,
                "credible_alleles": ",".join(str(allele) for allele in credible),
                "full_product_reads": full_product_reads,
                "repeat_informative_reads": repeat_informative_reads,
                "forward_strand_reads": strand_counts["forward"],
                "reverse_strand_reads": strand_counts["reverse"],
                "median_mapping_quality": round(median(mapping_qualities), 3) if mapping_qualities else "",
                "median_anchor_identity": round(median(anchor_identities), 6) if anchor_identities else "",
                "consensus_read_agreement": (
                    "" if primary_measurement.get("called_repeat_count", "") == ""
                    else "yes" if str(primary_measurement.get("called_repeat_count")) == str(best[0])
                    else "no"
                ),
            }
        )
    return rows
