from __future__ import annotations

import math


RESERVED_PROFILE_COLUMNS = {"profile_id", "strain_id", "metadata"}

PROFILE_MATCH_LOCUS_FIELDS = [
    "sample_id",
    "rank",
    "profile_id",
    "strain_id",
    "metadata",
    "locus_id",
    "called_repeat_count",
    "profile_repeat_count",
    "absolute_difference",
    "match_status",
    "profile_allele_probability",
]


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
        matched_locus_ids = []
        mismatched = []
        uncompared = []
        query_alleles = []
        profile_alleles = []
        locus_differences = []
        distance = 0.0
        compared = 0
        negative_log_likelihood = 0.0
        for locus in locus_columns:
            called = fingerprint.get(locus, "")
            expected = profile.get(locus, "")
            if called not in ("", None):
                query_alleles.append(f"{locus}={called}")
            if expected not in ("", None):
                profile_alleles.append(f"{locus}={expected}")
            if called in ("", None) or expected in ("", None):
                uncompared.append(locus)
                continue
            compared += 1
            try:
                delta = abs(float(called) - float(expected))
            except ValueError:
                delta = 0.0 if str(called) == str(expected) else 1.0
            if delta == 0:
                matched += 1
                matched_locus_ids.append(locus)
            else:
                mismatched.append(locus)
                distance += delta
            locus_differences.append(f"{locus}={round(delta, 6)}")
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
                "profile_id": profile.get("profile_id", ""),
                "best_profile_id": profile.get("profile_id", ""),
                "strain_id": profile.get("strain_id", ""),
                "metadata": profile.get("metadata", ""),
                "distance": round(distance, 4),
                "matched_loci": matched,
                "matched_locus_ids": ",".join(matched_locus_ids),
                "mismatched_loci": ",".join(mismatched),
                "uncompared_loci": ",".join(uncompared),
                "confidence": round(confidence, 6),
                "compared_loci": compared,
                "query_alleles": ";".join(query_alleles),
                "profile_alleles": ";".join(profile_alleles),
                "locus_differences": ";".join(locus_differences),
                "mean_negative_log_likelihood": (
                    round(mean_nll, 6) if math.isfinite(mean_nll) else ""
                ),
                "profile_probability_score": (
                    round(math.exp(-mean_nll), 6) if math.isfinite(mean_nll) else 0.0
                ),
            }
        )
    ranked = sorted(
        rows,
        key=lambda row: (
            row["distance"],
            -row["matched_loci"],
            float(row["mean_negative_log_likelihood"] or "inf"),
            str(row["best_profile_id"]),
        ),
    )
    for rank, row in enumerate(ranked, start=1):
        row["rank"] = rank
        row["is_best_match"] = "yes" if rank == 1 else "no"
    return ranked


def profile_match_locus_rows(
    sample_id: str,
    fingerprint: dict,
    profiles: list[dict],
    match_rows: list[dict],
    allele_rows: list[dict] | None = None,
) -> list[dict]:
    """Expand ranked profile comparisons to one machine-readable row per locus."""
    profile_by_id = {
        str(profile.get("profile_id", "")): profile for profile in profiles
    }
    probability_by_locus = {
        str(row["locus_id"]): _allele_probability_map(row)
        for row in allele_rows or []
    }
    rows = []
    for match in match_rows:
        profile_id = str(match.get("best_profile_id", ""))
        profile = profile_by_id.get(profile_id, {})
        for locus in (
            column
            for column in profile
            if column not in RESERVED_PROFILE_COLUMNS
        ):
            called = fingerprint.get(locus, "")
            expected = profile.get(locus, "")
            difference: float | str = ""
            status = "NOT_COMPARED"
            if called not in ("", None) and expected not in ("", None):
                try:
                    difference = abs(float(called) - float(expected))
                    status = "MATCH" if difference == 0 else "MISMATCH"
                except (TypeError, ValueError):
                    status = "MATCH" if str(called) == str(expected) else "MISMATCH"
                    difference = 0.0 if status == "MATCH" else 1.0
            expected_probability: float | str = ""
            if expected not in ("", None):
                try:
                    expected_probability = probability_by_locus.get(
                        locus, {}
                    ).get(float(expected), "")
                except (TypeError, ValueError):
                    expected_probability = ""
            rows.append(
                {
                    "sample_id": sample_id,
                    "rank": match.get("rank", ""),
                    "profile_id": profile_id,
                    "strain_id": profile.get("strain_id", ""),
                    "metadata": profile.get("metadata", ""),
                    "locus_id": locus,
                    "called_repeat_count": called,
                    "profile_repeat_count": expected,
                    "absolute_difference": difference,
                    "match_status": status,
                    "profile_allele_probability": expected_probability,
                }
            )
    return rows
