"""Shared locus-level allele inference for short- and long-read evidence."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass

from .alignment_evidence import CandidateEvidence
from .candidate_contexts import CandidateContext
from .models import Locus


@dataclass(frozen=True)
class InferenceThresholds:
    minimum_molecules: int = 3
    minimum_probability: float = 0.8
    minimum_margin: float = 0.2
    mixture_min_molecules: int = 2
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
                and mixture[1][1] >= thresholds.mixture_min_molecules
                and secondary_fraction >= thresholds.mixture_min_fraction
                and (technology != "illumina" or mixture[1][1] >= max(3, thresholds.mixture_min_molecules))
            ):
                status = "mixed"
            elif best_probability < thresholds.minimum_probability or best_probability - second_probability < thresholds.minimum_margin:
                status = "ambiguous"
            elif len(informative) < thresholds.minimum_molecules:
                status = "low_coverage"
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
        calls.append({
            "sample": sample_id,
            "locus": locus.locus_id,
            "repeat_count": best[0] if status in {"called", "low_coverage", "mixed"} else "",
            "status": status,
            "confidence": (
                "high" if status == "called" and best[1] >= 0.95
                else "moderate" if status == "called" else "low" if status == "low_coverage"
                else "unresolved"
            ),
            "best_probability": round(best[1], 8),
            "second_best_probability": round(second[1], 8),
            "margin": round(best[1] - second[1], 8),
            "molecule_support": len(informative),
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
        })
    return calls, molecule_calls


COMMON_LOCUS_CALL_FIELDS = [
    "sample", "locus", "repeat_count", "status", "confidence",
    "best_probability", "second_best_probability", "margin", "molecule_support",
    "direct_product_support", "full_span_support", "junction_support",
    "best_candidate_repeat", "candidate_distribution", "technology",
    "dominant_repeat", "secondary_repeat", "dominant_fraction", "secondary_fraction",
]