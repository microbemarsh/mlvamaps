from __future__ import annotations

import math


RESERVED_PROFILE_COLUMNS = {"profile_id", "strain_id", "metadata"}


def build_fingerprint(sample_id: str, allele_rows: list[dict], loci: list) -> tuple[list[dict], list[dict]]:
    by_locus = {row["locus_id"]: row for row in allele_rows}
    fingerprint = {"sample_id": sample_id}
    probabilistic = []
    for locus in loci:
        call = by_locus.get(locus.locus_id, {})
        value = call.get("called_repeat_count", "")
        fingerprint[locus.locus_id] = value
        probabilistic.append(
            {
                "sample_id": sample_id,
                "locus_id": locus.locus_id,
                "repeat_count": value,
                "posterior_probability": call.get("posterior_probability", 0.0),
            }
        )
    return [fingerprint], probabilistic


def _allele_probability_map(row: dict) -> dict[float, float]:
    probabilities = {}
    for entry in filter(None, str(row.get("allele_distribution", "")).split(";")):
        try:
            allele, probability = entry.rsplit(":", 1)
            probabilities[float(allele)] = float(probability)
        except (TypeError, ValueError):
            continue
    return probabilities


def match_profiles(
    sample_id: str,
    fingerprint: dict,
    profiles: list[dict],
    allele_rows: list[dict] | None = None,
) -> list[dict]:
    if not profiles:
        return []
    probability_by_locus = {
        str(row["locus_id"]): _allele_probability_map(row)
        for row in allele_rows or []
    }
    locus_columns = [col for col in profiles[0] if col not in RESERVED_PROFILE_COLUMNS]
    rows = []
    for profile in profiles:
        matched = 0
        mismatched = []
        distance = 0.0
        compared = 0
        negative_log_likelihood = 0.0
        for locus in locus_columns:
            called = fingerprint.get(locus, "")
            expected = profile.get(locus, "")
            if called in ("", None) or expected in ("", None):
                continue
            compared += 1
            try:
                delta = abs(float(called) - float(expected))
            except ValueError:
                delta = 0.0 if str(called) == str(expected) else 1.0
            if delta == 0:
                matched += 1
            else:
                mismatched.append(locus)
                distance += delta
            distribution = probability_by_locus.get(locus, {})
            if distribution:
                try:
                    probability = distribution.get(float(expected), 1e-9)
                except ValueError:
                    probability = 1.0 if str(called) == str(expected) else 1e-9
                negative_log_likelihood -= math.log(max(probability, 1e-9))
        confidence = matched / compared if compared else 0.0
        mean_nll = negative_log_likelihood / compared if compared else float("inf")
        rows.append(
            {
                "sample_id": sample_id,
                "best_profile_id": profile.get("profile_id", ""),
                "strain_id": profile.get("strain_id", ""),
                "distance": round(distance, 4),
                "matched_loci": matched,
                "mismatched_loci": ",".join(mismatched),
                "confidence": round(confidence, 6),
                "compared_loci": compared,
                "mean_negative_log_likelihood": (
                    round(mean_nll, 6) if math.isfinite(mean_nll) else ""
                ),
                "profile_probability_score": (
                    round(math.exp(-mean_nll), 6) if math.isfinite(mean_nll) else 0.0
                ),
            }
        )
    return sorted(
        rows,
        key=lambda row: (
            row["distance"],
            -row["matched_loci"],
            float(row["mean_negative_log_likelihood"] or "inf"),
            str(row["best_profile_id"]),
        ),
    )[:20]
