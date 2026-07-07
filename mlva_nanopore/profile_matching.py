from __future__ import annotations


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


def match_profiles(sample_id: str, fingerprint: dict, profiles: list[dict]) -> list[dict]:
    if not profiles:
        return []
    locus_columns = [col for col in profiles[0] if col not in RESERVED_PROFILE_COLUMNS]
    rows = []
    for profile in profiles:
        matched = 0
        mismatched = []
        distance = 0.0
        compared = 0
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
        confidence = matched / compared if compared else 0.0
        rows.append(
            {
                "sample_id": sample_id,
                "best_profile_id": profile.get("profile_id", ""),
                "strain_id": profile.get("strain_id", ""),
                "distance": round(distance, 4),
                "matched_loci": matched,
                "mismatched_loci": ",".join(mismatched),
                "confidence": round(confidence, 6),
            }
        )
    return sorted(rows, key=lambda row: (row["distance"], -row["matched_loci"]))[:20]
