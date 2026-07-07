from __future__ import annotations


def score_novelty(sample_id: str, allele_rows: list[dict], profile_matches: list[dict]) -> list[dict]:
    uncertain = sum(1 for row in allele_rows if float(row.get("posterior_probability") or 0) < 0.75)
    non_pass = sum(1 for row in allele_rows if row.get("call_status") not in ("PASS", "LOW_DEPTH"))
    loci = max(len(allele_rows), 1)
    uncertainty_component = (uncertain + non_pass) / (2 * loci)
    if profile_matches:
        nearest = profile_matches[0]
        distance = float(nearest.get("distance", 0.0))
        distance_component = min(1.0, distance / max(loci, 1))
        nearest_profile = nearest.get("best_profile_id", "")
    else:
        distance_component = 0.5
        nearest_profile = ""
    novelty_score = min(1.0, 0.65 * distance_component + 0.35 * uncertainty_component)
    if novelty_score < 0.25:
        interpretation = "known-like"
    elif novelty_score < 0.6:
        interpretation = "uncertain"
    else:
        interpretation = "potentially novel profile"
    return [
        {
            "sample_id": sample_id,
            "nearest_profile": nearest_profile,
            "novelty_score": round(novelty_score, 6),
            "interpretation": interpretation,
        }
    ]
