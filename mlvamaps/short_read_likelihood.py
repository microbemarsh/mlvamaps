"""Bounded conditional length inference from independently calibrated fragments."""

from collections import Counter, defaultdict
import math

from .calling import legacy_round_repeat_count, normalize_allele, repeat_unit_length


def _log_survival(value):
    if value < 10:
        return math.log(.5 * math.erfc(value / math.sqrt(2)))
    # Mills' expansion avoids erfc underflow in the calibrated library tails.
    inverse = 1 / value**2
    correction = 1 - inverse + 3 * inverse**2 - 15 * inverse**3 + 105 * inverse**4
    return -.5 * value**2 - math.log(value) - .5 * math.log(2 * math.pi) + math.log(correction)


def _log_normal_interval(low, high):
    if not low < high:
        return -math.inf
    if low < 0 < high:
        return math.log(.5 * (math.erf(high / math.sqrt(2)) - math.erf(low / math.sqrt(2))))
    first, second = ((_log_survival(low), _log_survival(high)) if low >= 0
                     else (_log_survival(-high), _log_survival(-low)))
    return first + math.log(-math.expm1(second - first)) if second < first else -math.inf


def _log_sum(values):
    maximum = max(values)
    if not math.isfinite(maximum):
        return maximum
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def _observation(metadata):
    kind = metadata.get("fragment_length_class")
    if kind is not None:
        fields = ("fragment_anchor_bp", "fragment_mate_length_bp", "fragment_boundary_min_bp",
                  "fragment_boundary_max_bp", "fragment_insert_mean_bp", "fragment_insert_sd_bp",
                  "fragment_allele_offset", "fragment_repeat_unit_bp")
        values = tuple(float(metadata[field]) for field in fields)
        anchor, length, low, high, mean, sigma, offset, unit = values
        opposite = float(metadata["fragment_opposite_anchor_bp"]) if kind == "boundary_pair" else 0.
        if (metadata.get("fragment_length_scheme") != "anchor_union_v1"
                or kind not in {"boundary_pair", "repeat_mate"} or not all(map(math.isfinite, values + (opposite,)))
                or anchor < 0 or length <= 0 or not 0 <= low < high <= length
                or mean <= 0 or sigma <= 0 or unit <= 0
                or not all(value.is_integer() for value in (anchor, length, low, high, opposite))
                or (kind == "boundary_pair" and not low <= opposite <= high)):
            raise ValueError("Invalid paired-fragment class geometry")
        return (kind, *values, opposite)
    value, sigma = float(metadata["pair_repeat_count"]), float(metadata["pair_repeat_sigma"])
    low, high = metadata.get("pair_selection_low_delta"), metadata.get("pair_selection_high_delta")
    if not math.isfinite(value) or not math.isfinite(sigma) or sigma <= 0:
        raise ValueError("Invalid paired-fragment count or uncertainty")
    if (low is None) != (high is None) or (low is not None and (
            not math.isfinite(low) or not math.isfinite(high) or not low <= 0 <= high or low >= high)):
        raise ValueError("Invalid paired-fragment selection interval")
    return ("gaussian", value, sigma) if low is None else ("selected_pair", value, sigma, low, high)


def _observation_score(state, observation, null_fraction=0.):
    kind, *values = observation
    if kind in {"gaussian", "selected_pair"}:
        value, sigma = values[:2]
        delta = (state - value) / sigma
        score = -.5 * delta**2 - math.log(sigma * math.sqrt(2 * math.pi))
        if kind == "selected_pair":
            mass = _log_normal_interval(delta + values[2] / sigma, delta + values[3] / sigma)
            return score - mass if math.isfinite(mass) else -math.inf
        return score
    anchor, length, low, high, mean, sigma, offset, unit, opposite = values
    repeat = unit * (state - offset)
    if repeat < 0:
        return -math.inf
    # Condition on the common retained set: a pure-repeat opposite mate OR a
    # mate crossing the opposite boundary. Other read classes are conditioned
    # out, rather than silently counted as missing observations.
    # Reads and insert sizes occupy integer base pairs. Integrating half-base
    # bins avoids shifting the inferred allele at inclusive selection edges.
    pure = _log_normal_interval((anchor + length - .5 - mean) / sigma,
                                (anchor + repeat + .5 - mean) / sigma)
    boundary_low = max(low, length - repeat)
    boundary = _log_normal_interval((anchor + repeat + boundary_low - .5 - mean) / sigma,
                                    (anchor + repeat + high + .5 - mean) / sigma)
    normalizer = _log_sum((pure, boundary))
    if kind == "repeat_mate":
        signal = pure - normalizer if math.isfinite(normalizer) else -math.inf
        null = math.log(.5)
    else:
        span = anchor + repeat + opposite
        density = _log_normal_interval((span - .5 - mean) / sigma, (span + .5 - mean) / sigma)
        signal = density - normalizer if math.isfinite(normalizer) and opposite >= boundary_low else -math.inf
        null = math.log(.5 / (high - low + 1))
    return signal if not null_fraction else _log_sum((math.log1p(-null_fraction) + signal,
                                                     math.log(null_fraction) + null))


def infer_short_read_length(rows, locus, thresholds, *, round_tolerance=.25, sample_mode="isolate"):
    """Return estimated call updates, with no measured or classification votes.

    Unique molecules contribute one likelihood. Reference/mate duplicates are
    collapsed, tied geometry alternatives share that molecule's likelihood,
    and partial observations constrain lengths without adding votes. The new
    two-class model requires independently calibrated inserts and deterministic
    strong-anchor selection. Legacy opposite-flank pairs can use their own
    target-selection correction. Pure-repeat mates alone remain unbounded.
    """
    molecules = defaultdict(list)
    for row in rows:
        molecules[row.molecule_id].append(row)
    usable = {name: values for name, values in molecules.items()
              if not any("ambiguous_locus" in row.quality_flags for row in values)}
    if any(row.measured_repeat_count is not None for values in usable.values() for row in values):
        return {}
    bounds, observations = {}, {}
    conditioned = False
    for name, values in usable.items():
        for row in values:
            for field in ("repeat_lower_bound", "pair_repeat_lower_bound"):
                bound = row.metadata.get(field)
                if bound is not None and math.isfinite(bound) and bound > 0:
                    bounds[name] = max(bound, bounds.get(name, 0))
        candidates = [row for row in values if "fragment_length_class" in row.metadata]
        candidates = candidates or [row for row in values if "pair_repeat_count" in row.metadata
                                    and "pair_repeat_sigma" in row.metadata]
        if not candidates:
            continue
        best = max(row.metadata.get("pair_alignment_score", 0) for row in candidates)
        alternatives = {_observation(row.metadata) for row in candidates
                        if row.metadata.get("pair_alignment_score", 0) == best}
        # Different observed classes for one molecule are not densities of the
        # same event. Keep its partial bounds without choosing a class by copies.
        if len({item[0] for item in alternatives}) != 1:
            continue
        observations[name] = tuple(sorted(alternatives))
        conditioned |= any(item[0] != "gaussian" for item in alternatives)
    lower = max(bounds.values(), default=0.)
    pairs = {name: values for name, values in observations.items() if values[0][0] != "repeat_mate"}
    rich = len(observations) - len(pairs)
    diagnostics = dict(length_model_status="lower_bound_only", finite_length_support=False,
                       upper_unbounded=True, paired_span_support=len(pairs), partial_bound_support=len(bounds),
                       repeat_rich_mate_support=rich)
    if sample_mode != "isolate" and conditioned:
        return {"length_model_status": "unsupported_sample_mode"}
    if not pairs:
        return {**diagnostics, "repeat_lower_bound": lower or ""} if bounds or rich else {}
    if not conditioned:
        return {}  # Public callers without selection geometry retain their existing fit.
    groups = Counter(observations.values())
    centers = []
    for values in pairs.values():
        modes = [(item[7] + (item[5] - item[1] - item[9]) / item[8], item[6] / item[8])
                 if item[0] == "boundary_pair" else (item[1], item[2]) for item in values]
        center = sum(value for value, _ in modes) / len(modes)
        variance = sum(sigma**2 + (value - center)**2 for value, sigma in modes) / len(modes)
        centers.append((center, math.sqrt(variance)))
    union = [item for values in groups for item in values if item[0] in {"boundary_pair", "repeat_mate"}]
    # ponytail: profile two contamination fractions instead of forcing a null
    # fraction onto clean data; calibrate this bounded grid on matched samples.
    null_fractions = (0., .01) if union else (0.,)
    cells = len(null_fractions) * sum(len(values) for values in groups)
    budget = min(thresholds.maximum_candidates, thresholds.maximum_likelihood_cells // max(1, cells))
    if budget < 3:
        precision = sum(1 / sigma**2 for _, sigma in centers)
        raw_best = max(lower, sum(value / sigma**2 for value, sigma in centers) / precision)
        best = legacy_round_repeat_count(raw_best, round_tolerance)
        return {**diagnostics, "length_model_status": "insufficient_budget",
                "repeat_count": best, "repeat_count_raw": raw_best, "best_candidate_repeat": best,
                "status": "estimated", "confidence": "low", "best_probability": 0.,
                "second_best_probability": 0., "margin": 0., "candidate_distribution": "",
                "repeat_interval_low": lower, "repeat_interval_high": "", "repeat_lower_bound": lower or "",
                "estimate_support": len(observations), "dominant_repeat": best, "dominant_fraction": 1.,
                "secondary_repeat": "", "secondary_fraction": "", "allele_fractions": f"{best}:1.00000000",
                "inference_method": "bounded_short_read_length_fallback",
                "probability_kind": "unscored_length_fallback", "estimation_method": "joint_short_read_length",
                "estimation_reason": "Insufficient budget for selection-aware likelihood; retain a lower-bound-constrained fragment mean with unbounded uncertainty and zero model support.",
                "computational_fallback": "short_read_length_grid_limit"}
    maximum_sigma = max(sigma for _, sigma in centers)
    high = max([lower + 8 * maximum_sigma] + [value + 8 * sigma for value, sigma in centers]
               + [item[7] + (max(item[5] - item[1], item[2]) + 8 * item[6]) / item[8] for item in union])
    if not math.isfinite(high) or high <= lower:
        return {**diagnostics, "length_model_status": "invalid_search_range"}
    step = 1 / max(1, repeat_unit_length(locus))
    requested = math.ceil((high - lower) / step) + 1
    count = min(budget, requested)
    limited = requested > budget
    step = (high - lower) / (count - 1) if limited else step
    states = [normalize_allele(lower + index * step) for index in range(count)]
    scores, nulls = [], []
    for state in states:
        fits = [(sum(number * (_log_sum([_observation_score(state, item, fraction) for item in alternatives])
                               - math.log(len(alternatives))) for alternatives, number in groups.items()), fraction)
                for fraction in null_fractions]
        score, fraction = max(fits, key=lambda item: (item[0], -item[1]))
        scores.append(score)
        nulls.append(fraction)
    finite = [(state, score) for state, score in zip(states, scores) if math.isfinite(score)]
    if not finite:
        return {**diagnostics, "length_model_status": "no_finite_likelihood"}
    raw_best, best_score = max(finite, key=lambda item: (item[1], -item[0]))
    supported = [state for state, score in finite if score >= best_score - 2]
    upper_unbounded = states[-1] in supported
    if union and all(item[0] in {"boundary_pair", "repeat_mate"} for values in groups for item in values):
        # A contamination-supported long-array plateau may remain plausible
        # beyond the finite grid. Its analytic limit must not become a false
        # finite upper bound, even under a coarse computational grid.
        tail = max(sum(number * (_log_sum([
            math.log1p(-fraction / 2) if item[0] == "repeat_mate" else
            math.log(fraction * .5 / (item[4] - item[3] + 1)) if fraction else -math.inf
            for item in alternatives]) - math.log(len(alternatives)))
            for alternatives, number in groups.items()) for fraction in null_fractions)
        upper_unbounded |= tail >= best_score - 2
    probabilities = defaultdict(float)
    weights = [math.exp(score - best_score) for _, score in finite]
    total = sum(weights)
    for (state, _), weight in zip(finite, weights):
        probabilities[legacy_round_repeat_count(state, round_tolerance)] += weight / total
    best = legacy_round_repeat_count(raw_best, round_tolerance)
    ranked = sorted(probabilities.items(), key=lambda item: (-item[1], item[0]))
    ranked.sort(key=lambda item: item[0] != best)
    probability = probabilities[best]
    alternative = max((value for state, value in ranked if state != best), default=0.)
    reason = ("Independent calibrated fragments scored conditionally on eligible boundary/repeat-mate classes "
              "or opposite-flank target selection, jointly constrained by anchored partial repeats. "
              "The isolate model compares 0% and 1% background assumptions; support and likelihood intervals "
              "are uncalibrated. Unassigned repeat-only reads do not supply locus abundance.")
    if upper_unbounded:
        reason += " Lengths at or beyond the upper search edge remain plausible; no finite upper length is established."
    return {**diagnostics, "length_model_status": "estimated",
            "finite_length_support": not upper_unbounded, "upper_unbounded": upper_unbounded,
            "length_null_fraction": nulls[states.index(raw_best)],
            "candidate_grid_low": states[0], "candidate_grid_high": states[-1],
            "repeat_count": best, "repeat_count_raw": raw_best, "best_candidate_repeat": best,
            "status": "estimated", "confidence": "low", "inference_tier": 3,
            "best_probability": round(probability, 8), "second_best_probability": round(alternative, 8),
            "margin": round(probability - alternative, 8),
            "candidate_distribution": ";".join(f"{state}:{value:.8f}" for state, value in ranked),
            "repeat_interval_low": min(supported), "repeat_interval_high": "" if upper_unbounded else max(supported),
            "repeat_lower_bound": lower or "", "estimate_support": len(observations),
            "dominant_repeat": best, "dominant_fraction": 1., "secondary_repeat": "", "secondary_fraction": "",
            "allele_fractions": f"{best}:1.00000000",
            "inference_method": "short_read_joint_length_likelihood", "probability_kind": "conditional_fragment_length_support",
            "estimation_method": "joint_short_read_length", "estimation_reason": reason,
            "computational_fallback": "short_read_length_grid_limit" if limited else ""}
