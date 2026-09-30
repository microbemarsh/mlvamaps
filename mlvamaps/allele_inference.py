"""Shared locus-level allele inference for short- and long-read evidence."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass

from .alignment_evidence import CandidateEvidence
from .candidate_contexts import CandidateContext
from .models import Locus
from .repeat_calibration import repeat_unit_length


def retain_detected_estimate(call: dict, locus: Locus, depth_estimate: dict | None = None) -> dict:
    """Retain a point estimate without claiming unsupported exact confidence.

    Depth ratios and panel priors have no validated probability of correctness.
    Their zero confidence means unassessed, not a probability that the count is
    wrong. They never generate sequence or become definitive calls.
    """
    if call['status'] == 'not_found' or call.get('repeat_count') not in ('', None):
        return call
    estimate = (depth_estimate or {}).get('repeat_count')
    reason = call.get('reason', '')
    lower = call.get('repeat_count_min')
    if (estimate is not None and 'repeat_length_lower_bound' in reason
            and lower not in ('', None) and estimate < float(lower)):
        estimate = None
        reason += '; depth_estimate_below_observed_lower_bound'
    if estimate is not None:
        interval = depth_estimate['repeat_interval']
        size = depth_estimate['amplicon_length']
        call.update(product_size_bp=size,
                    repeat_count_raw=depth_estimate.get('repeat_count_raw', ''),
                    inference_method='KMER_DEPTH', confidence_kind='unvalidated_depth',
                    repeat_count_interval_kind='coverage_sensitivity')
        reason = '; '.join(filter(None, (reason, 'depth_estimate_unvalidated_locus_ownership_and_coverage')))
    else:
        if not repeat_unit_length(locus) and not locus.nominal_repeat_units:
            call.update(confidence=0.0, best_probability=0.0, confidence_kind='missing_repeat_definition',
                        reason='; '.join(filter(None, (reason, 'missing_repeat_unit_or_nominal_count'))))
            return call
        nominal = bool(locus.expected_product_size_bp or locus.nominal_repeat_units)
        estimate = (locus.nominal_repeat_units if nominal else
                    (locus.expected_min_repeats + locus.expected_max_repeats) / 2)
        lower, upper = call.get('repeat_count_min'), call.get('repeat_count_max')
        # A censored boundary constrains a prior, but the search ceiling is not
        # an observed upper bound and must not truncate the nominal estimate.
        if lower not in ('', None):
            estimate = max(estimate, float(lower))
        if upper not in ('', None) and 'repeat_length_lower_bound' not in reason:
            estimate = min(estimate, float(upper))
        interval = (lower if lower not in ('', None) else min(estimate, locus.expected_min_repeats),
                    '' if 'repeat_length_lower_bound' in reason else
                    max(estimate, float(upper)) if upper not in ('', None) else max(estimate, locus.expected_max_repeats))
        call.update(repeat_count_raw='', inference_method='PANEL_PRIOR', confidence_kind='prior_only',
                    repeat_count_interval_kind='lower_bound' if interval[1] == '' else 'prior_range')
        reason = '; '.join(filter(None, (reason, 'nominal_count_prior' if nominal else 'panel_range_midpoint_prior')))
    call.update(repeat_count=estimate, status='estimated', confidence=0.0, best_probability=0.0,
                margin=0.0, best_candidate_repeat=estimate, dominant_repeat=estimate,
                second_best_repeat_count='', second_best_probability=0.0,
                candidate_distribution='', repeat_count_min=interval[0], repeat_count_max=interval[1],
                reason=reason)
    return call


@dataclass(frozen=True)
class InferenceThresholds:
    minimum_molecules: int = 0  # Compatibility only; no molecule-count cutoff.
    minimum_probability: float = 0.8
    minimum_margin: float = 0.2
    mixture_min_molecules: int = 0  # Compatibility only.
    mixture_min_fraction: float = 0.2
    temperature: float = 8.0


def _softmax(scores: dict[int | float, float]) -> dict[int | float, float]:
    if not scores:
        return {}
    maximum = max(scores.values())
    weights = {state: math.exp(max(-700.0, score - maximum)) for state, score in scores.items()}
    total = sum(weights.values())
    return {state: weight / total for state, weight in weights.items()}


def _tier_weight(evidence: CandidateEvidence) -> float:
    if evidence.direct_product_measurement:
        return 12.0
    if evidence.full_repeat_span:
        return 8.0
    if evidence.technology == "illumina" and not evidence.pair_geometry_support:
        # One repeat boundary establishes presence, not the unseen tract length.
        # Repeating that observation must not select a truncated candidate.
        return 0.0
    if evidence.left_boundary_span and evidence.right_boundary_span:
        return 5.0
    if evidence.left_boundary_span or evidence.right_boundary_span:
        return 2.5
    if evidence.repeat_indel_support:
        return 2.0
    if evidence.pair_geometry_support:
        return 1.5
    return 0.0  # Generic mapping establishes presence, never an allele.


def _molecule_distribution(
    rows: list[CandidateEvidence], states: list[int | float], temperature: float
) -> tuple[dict[int | float, float], int | float | None, float]:
    direct = [row.measured_repeat_count for row in rows if row.measured_repeat_count is not None]
    if direct:
        measured = Counter(direct).most_common(1)[0][0]
        scores = {state: -24.0 * abs(float(state) - float(measured)) for state in states}
        return _softmax(scores), measured, max(_tier_weight(row) for row in rows)
    # Max over duplicate reference contexts prevents reference-copy bias.
    by_state: dict[int | float, float] = {}
    weights: dict[int | float, float] = {}
    for row in rows:
        weight = _tier_weight(row)
        if weight <= 0:
            continue
        score = (
            row.alignment_score / max(temperature, 1e-6)
            + 2.0 * row.repeat_indel_support
            + row.pair_geometry_support
        )
        if score > by_state.get(row.repeat_count, -math.inf):
            by_state[row.repeat_count] = score
            weights[row.repeat_count] = weight
    if not by_state:
        return {}, None, 0.0
    probabilities = _softmax({state: by_state.get(state, min(by_state.values()) - 20.0) for state in states})
    winner = max(probabilities, key=lambda state: (probabilities[state], -float(state)))
    return probabilities, winner, weights.get(winner, 1.0)


def infer_alleles(
    evidence: list[CandidateEvidence],
    loci: list[Locus],
    contexts: list[CandidateContext],
    sample_id: str,
    technology: str,
    thresholds: InferenceThresholds | None = None,
) -> tuple[list[dict[str, object]], dict[tuple[str, str], int | float]]:
    """Infer repeat states while preserving presence, ambiguity, and low depth."""
    thresholds = thresholds or InferenceThresholds()
    states_by_locus: dict[str, list[int | float]] = defaultdict(list)
    for context in contexts:
        if context.repeat_count is not None and context.repeat_count not in states_by_locus[context.locus_id]:
            states_by_locus[context.locus_id].append(context.repeat_count)
    # Direct observations need not lie on the reference candidate grid.
    for row in evidence:
        if row.measured_repeat_count is not None and row.measured_repeat_count not in states_by_locus[row.locus_id]:
            states_by_locus[row.locus_id].append(row.measured_repeat_count)
    for states in states_by_locus.values():
        states.sort(key=float)
    by_locus_molecule: dict[tuple[str, str], list[CandidateEvidence]] = defaultdict(list)
    for row in evidence:
        by_locus_molecule[(row.locus_id, row.molecule_id)].append(row)

    calls: list[dict[str, object]] = []
    molecule_calls: dict[tuple[str, str], int | float] = {}
    for locus in loci:
        locus_items = {
            molecule: rows for (locus_id, molecule), rows in by_locus_molecule.items()
            if locus_id == locus.locus_id
        }
        states = states_by_locus[locus.locus_id]
        informative: list[tuple[dict[int | float, float], int | float, float, list[CandidateEvidence]]] = []
        for molecule, rows in locus_items.items():
            distribution, winner, weight = _molecule_distribution(rows, states, thresholds.temperature)
            if distribution and winner is not None:
                informative.append((distribution, winner, weight, rows))
                molecule_calls[(locus.locus_id, molecule)] = winner

        if not locus_items:
            ranked: list[tuple[int | float, float]] = []
            status = "not_found"
        elif not informative:
            ranked = []
            status = "detected_unresolved"
        else:
            log_scores = {state: 0.0 for state in states}
            for distribution, _winner, weight, _rows in informative:
                for state in states:
                    log_scores[state] += weight * math.log(max(distribution.get(state, 1e-12), 1e-12))
            posterior = _softmax(log_scores)
            ranked = sorted(posterior.items(), key=lambda item: (-item[1], float(item[0])))
            best_probability = ranked[0][1]
            second_probability = ranked[1][1] if len(ranked) > 1 else 0.0
            molecule_counts = Counter(winner for _distribution, winner, _weight, _rows in informative)
            mixture = molecule_counts.most_common()
            total = sum(molecule_counts.values())
            secondary_fraction = mixture[1][1] / total if len(mixture) > 1 else 0.0
            if (
                len(mixture) > 1
                and secondary_fraction >= thresholds.mixture_min_fraction
            ):
                status = "mixed"
            elif best_probability < thresholds.minimum_probability or best_probability - second_probability < thresholds.minimum_margin:
                status = "ambiguous"
            else:
                status = "called"

        best = ranked[0] if ranked else ("", 0.0)
        second = ranked[1] if len(ranked) > 1 else ("", 0.0)
        direct = {
            row.molecule_id for rows in locus_items.values() for row in rows
            if row.direct_product_measurement
        }
        full = {
            row.molecule_id for rows in locus_items.values() for row in rows
            if row.full_repeat_span
        }
        junction = {
            row.molecule_id for rows in locus_items.values() for row in rows
            if row.left_boundary_span or row.right_boundary_span
        }
        molecule_counts = Counter(
            winner for _distribution, winner, _weight, _rows in informative
        )
        mixture = molecule_counts.most_common()
        total_informative = sum(molecule_counts.values())
        calls.append(retain_detected_estimate({
            "sample": sample_id,
            "locus": locus.locus_id,
            "repeat_count": best[0],
            "status": status,
            "confidence": round(best[1], 8),
            "confidence_kind": "conditional_posterior",
            "best_probability": round(best[1], 8),
            "second_best_probability": round(second[1], 8),
            "margin": round(best[1] - second[1], 8),
            "molecule_support": len(locus_items),
            "direct_product_support": len(direct),
            "full_span_support": len(full),
            "junction_support": len(junction),
            "best_candidate_repeat": best[0],
            "candidate_distribution": ";".join(f"{state}:{probability:.8f}" for state, probability in ranked),
            "technology": technology,
            "dominant_repeat": mixture[0][0] if mixture else "",
            "secondary_repeat": mixture[1][0] if len(mixture) > 1 else "",
            "dominant_fraction": round(mixture[0][1] / total_informative, 8) if mixture else "",
            "secondary_fraction": round(mixture[1][1] / total_informative, 8) if len(mixture) > 1 else "",
        }, locus))
    return calls, molecule_calls


COMMON_LOCUS_CALL_FIELDS = [
    "sample", "locus", "repeat_count", "repeat_count_raw", "status", "confidence",
    "best_probability", "second_best_probability", "margin", "molecule_support",
    "direct_product_support", "full_span_support", "junction_support",
    "best_candidate_repeat", "candidate_distribution", "technology",
    "dominant_repeat", "secondary_repeat", "dominant_fraction", "secondary_fraction",
    "product_size_bp", "inference_method", "reason", "repeat_count_min", "repeat_count_max",
    "confidence_kind", "repeat_count_interval_kind",
]
