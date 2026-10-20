"""Shared locus-level allele inference for short- and long-read evidence."""

from __future__ import annotations

import math
from collections import Counter, defaultdict
from dataclasses import dataclass

import numpy as np

from .alignment_evidence import CandidateEvidence
from .calling import legacy_round_repeat_count, normalize_allele, repeat_unit_length
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
    rows: list[CandidateEvidence], states: list[int | float], temperature: float,
    require_measurement: bool = False,
) -> tuple[dict[int | float, float], int | float | None, float]:
    direct = [row.measured_repeat_count for row in rows if row.measured_repeat_count is not None]
    if require_measurement and len(set(direct)) != 1:
        # An observed reference is a recruitment target, not a measured allele.
        # Conflicting measurements of the same molecule cannot vote by copies.
        return {}, None, 0.0
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
    *, require_measurement: bool = False, sample_mode: str = "metagenome",
    round_tolerance: float = 0.25,
    consensus_measurements: dict | None = None,
) -> tuple[list[dict[str, object]], dict[tuple[str, str], int | float]]:
    """Infer repeat states while preserving presence, ambiguity, and low depth."""
    if sample_mode not in {"isolate", "metagenome"}:
        raise ValueError("sample_mode must be 'isolate' or 'metagenome'")
    if require_measurement:
        return _infer_measurements(
            evidence, loci, sample_id, technology, thresholds or InferenceThresholds(),
            sample_mode, round_tolerance, consensus_measurements or {},
        )
    thresholds = thresholds or InferenceThresholds()
    states_by_locus: dict[str, list[int | float]] = defaultdict(list)
    for context in contexts:
        if context.repeat_count is not None and context.repeat_count not in states_by_locus[context.locus_id]:
            states_by_locus[context.locus_id].append(context.repeat_count)
    by_locus_molecule: dict[tuple[str, str], list[CandidateEvidence]] = defaultdict(list)
    for row in evidence:
        by_locus_molecule[(row.locus_id, row.molecule_id)].append(row)
        # Mapping candidates recruit molecules; they must not censor measured
        # half alleles or products outside the reference bank's repeat range.
        measured = row.measured_repeat_count
        states = states_by_locus[row.locus_id]
        if measured is not None and measured not in states:
            states.append(measured)
    for states in states_by_locus.values():
        states.sort(key=float)

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
            distribution, winner, weight = _molecule_distribution(
                rows, states, thresholds.temperature, require_measurement,
            )
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
            "presence_support": len(locus_items),
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


def _infer_measurements(evidence, loci, sample_id, technology, thresholds,
                        sample_mode, round_tolerance, consensus_measurements):
    """Fit raw molecule lengths, then apply the assembly allele convention.

    Identical observations share a likelihood row and retain their molecule
    multiplicity. Reference copies, mates and derived consensus add no votes.
    """
    from .mapping_classification import fit_reference_mixture

    grouped = defaultdict(lambda: defaultdict(list))
    for row in evidence:
        grouped[row.locus_id][row.molecule_id].append(row)
    calls, molecule_calls = [], {}
    for locus in loci:
        molecules = grouped[locus.locus_id]
        observations = {}
        for molecule, rows in molecules.items():
            measured = [row for row in rows if row.measured_repeat_count is not None]
            if not measured or len({row.measured_repeat_count for row in measured}) != 1:
                continue
            # Complete products use the assembly scale; flank-only observations
            # cannot replace that scale on the same molecule.
            selected = [row for row in measured if row.direct_product_measurement] or measured
            raw_values = {row.raw_repeat_count for row in selected if row.raw_repeat_count is not None}
            if len(raw_values) > 1:
                continue
            raw = next(iter(raw_values)) if raw_values else selected[0].measured_repeat_count
            sigma = max((row.repeat_sigma if row.repeat_sigma is not None else 0.08
                         for row in selected), default=0.08)
            if not math.isfinite(raw) or raw < 0 or not math.isfinite(sigma) or sigma <= 0:
                raise ValueError("Molecule repeat measurements must be finite, nonnegative, with positive uncertainty")
            observations[molecule] = (float(raw), float(sigma))
            # Classification receives the original measurement, never a fitted
            # allele assignment or a consensus-derived second observation.
            molecule_calls[(locus.locus_id, molecule)] = selected[0].measured_repeat_count

        counts = Counter(observations.values())
        measured_support = len(observations)
        # Library geometry is an optional fallback, never a replacement for
        # directly measured products and never classification ground truth.
        if not observations:
            for molecule, rows in molecules.items():
                pairs = [row for row in rows if "pair_repeat_count" in row.metadata]
                if not pairs:
                    continue
                score = max(row.metadata["pair_alignment_score"] for row in pairs)
                values = {(row.metadata["pair_repeat_count"], row.metadata["pair_repeat_sigma"])
                          for row in pairs if row.metadata["pair_alignment_score"] == score}
                if len(values) == 1:
                    value, sigma = next(iter(values))
                    if not math.isfinite(value) or not math.isfinite(sigma) or sigma <= 0:
                        raise ValueError("Invalid paired-fragment length observation")
                    observations[molecule] = (value, sigma)
            counts = Counter(observations.values())
        pair_only = bool(observations) and not measured_support
        if pair_only:
            # Gaussian insert likelihoods have sufficient statistics: no dense
            # pair-by-allele matrix is needed even for a deeply covered library.
            # ponytail: target-selection bias is not modeled. Keep ESTIMATED;
            # use an ascertainment-aware likelihood before raising confidence.
            precision = sum(1 / sigma**2 for _, sigma in observations.values())
            mean = sum(value / sigma**2 for value, sigma in observations.values()) / precision
            counts = Counter({(mean, 1 / math.sqrt(precision)): 1})
        ranked, fractions = [], []
        raw_best, best_count = "", ""
        iterations, converged = 0, True
        consensus_status = "not_available"
        status = "not_found" if not molecules else "detected_unresolved"
        if counts:
            unit = max(repeat_unit_length(locus), 1)
            states = set()
            for raw, sigma in counts:
                # Local base-pair grids retain novel/fractional alleles without
                # expanding every allele between distant observed products.
                # Pair uncertainty can span many repeats; half-unit steps keep
                # that fallback bounded without pretending to base-pair precision.
                step = max(1 / unit, .5, 3 * sigma / 200) if pair_only else 1 / unit
                radius = max(1, math.ceil(3 * sigma / step))
                center = max(0, raw)
                states.update(normalize_allele(center + delta * step)
                              for delta in range(-radius, radius + 1)
                              if center + delta * step >= 0)
            consensus = consensus_measurements.get(locus.locus_id, {})
            candidate = (consensus.get("raw_repeat_count")
                         if consensus.get("supporting_reads", 0) >= thresholds.minimum_molecules else None)
            if candidate not in (None, "") and math.isfinite(float(candidate)) and float(candidate) >= 0:
                states.add(normalize_allele(candidate))
                consensus_status = "candidate_only"
            states = sorted(states)
            raw = np.asarray([item[0] for item in counts])[:, None]
            sigma = np.asarray([item[1] for item in counts])[:, None]
            multiplicities = np.asarray(list(counts.values()), dtype=float)
            delta = np.abs(raw - np.asarray(states)[None, :])
            # ponytail: a Gaussian core plus 1% broad error component handles
            # occasional indel outliers; technology-specific calibration needs
            # matched real data before these model supports are probabilities.
            core = -0.5 * (delta / sigma) ** 2 - np.log(sigma * math.sqrt(2 * math.pi))
            tail = -delta / np.maximum(1.0, sigma * 4) - np.log(2 * np.maximum(1.0, sigma * 4))
            log_likelihood = np.logaddexp(math.log(.99) + core, math.log(.01) + tail)
            joint = multiplicities @ log_likelihood
            posterior = _softmax(dict(zip(states, joint)))
            raw_best = max(states, key=lambda state: (posterior[state], -state))
            best_count = legacy_round_repeat_count(raw_best, round_tolerance)
            grouped_support = defaultdict(float)
            for state, probability in posterior.items():
                grouped_support[legacy_round_repeat_count(state, round_tolerance)] += probability
            ranked = sorted(grouped_support.items(), key=lambda item: (-item[1], item[0]))
            # Keep the selected raw mode and its rounded allele together.
            ranked.sort(key=lambda item: item[0] != best_count)
            mixture = False
            unresolved_mixture = False
            if sample_mode == "metagenome" and len(counts) > 1 and not pair_only:
                fitted, history = fit_reference_mixture(
                    log_likelihood, multiplicities, max_iterations=200, tolerance=1e-6,
                )
                iterations = len(history) - 1
                responsibilities = np.exp(log_likelihood - log_likelihood.max(axis=1, keepdims=True)) * fitted
                responsibilities /= responsibilities.sum(axis=1, keepdims=True)
                updated = multiplicities @ responsibilities / multiplicities.sum()
                converged = bool(np.max(np.abs(updated - fitted)) < 1e-6)
                grouped_fractions = defaultdict(float)
                for state, fraction in zip(states, fitted):
                    grouped_fractions[legacy_round_repeat_count(state, round_tolerance)] += float(fraction)
                fractions = sorted(grouped_fractions.items(), key=lambda item: (-item[1], item[0]))
                minimum = max(3 if technology == "illumina" else 2, thresholds.mixture_min_molecules)
                supported = [(allele, fraction) for allele, fraction in fractions
                             if fraction >= thresholds.mixture_min_fraction
                             and fraction * len(observations) >= minimum]
                # Require a penalized likelihood gain over the single-allele
                # model. Diffuse/identical mappings alone cannot establish a
                # mixture; many individually ambiguous molecules still can.
                parameters = max(1, int(np.count_nonzero(fitted * len(observations) >= 1)) - 1)
                gain = 2 * (history[-1] - float(joint.max()))
                mixture = (converged and len(supported) > 1
                           and gain > parameters * math.log(max(2, len(observations))))
                unresolved_mixture = not converged and len(supported) > 1
                if mixture:
                    best_count = supported[0][0]
                    raw_best = max((state for state in states
                                    if legacy_round_repeat_count(state, round_tolerance) == best_count),
                                   key=lambda state: fitted[states.index(state)])
            else:
                fractions = [(best_count, 1.0)]
            probability = dict(ranked).get(best_count, 0.0)
            alternative = max((value for state, value in ranked if state != best_count), default=0.0)
            if mixture:
                status = "mixed"
            elif unresolved_mixture or probability < thresholds.minimum_probability or probability - alternative < thresholds.minimum_margin:
                status = "ambiguous"
            elif len(observations) < thresholds.minimum_molecules:
                status = "low_coverage"
            else:
                status = "called"
            if consensus_status == "candidate_only":
                # Consensus proposes a sequence-supported state, but only the
                # original molecule likelihoods determine its support.
                consensus_status = ("agrees" if legacy_round_repeat_count(float(candidate), round_tolerance) == best_count
                                    else "disagrees")
                consensus_allele = legacy_round_repeat_count(float(candidate), round_tolerance)
                if (status == "ambiguous"
                        and grouped_support.get(consensus_allele, 0) >= max(grouped_support.values()) / 4):
                    best_count, raw_best = consensus_allele, float(candidate)
                    status, consensus_status = "estimated", "resolved_ambiguity"
                    if sample_mode == "isolate":
                        fractions = [(best_count, 1.0)]
            if pair_only:
                status = "estimated"

        ranked.sort(key=lambda item: item[0] != best_count)
        best_probability = dict(ranked).get(best_count, 0.0)
        second_probability = max((prob for state, prob in ranked if state != best_count), default=0.0)
        rows = [row for items in molecules.values() for row in items]
        calls.append({
            "sample": sample_id, "locus": locus.locus_id,
            "repeat_count": best_count if status in {"called", "low_coverage", "mixed", "estimated"} else "",
            "repeat_count_raw": raw_best,
            "status": status,
            "confidence": ("high" if status == "called" and best_probability >= .95
                           else "moderate" if status == "called"
                           else "low" if status in {"low_coverage", "estimated"} else "unresolved"),
            "best_probability": round(best_probability, 8),
            "second_best_probability": round(second_probability, 8),
            "margin": round(best_probability - second_probability, 8),
            "molecule_support": measured_support, "presence_support": len(molecules),
            "direct_product_support": len({row.molecule_id for row in rows if row.direct_product_measurement}),
            "full_span_support": len({row.molecule_id for row in rows if row.full_repeat_span}),
            "junction_support": len({row.molecule_id for row in rows if row.left_boundary_span or row.right_boundary_span}),
            "best_candidate_repeat": best_count,
            "candidate_distribution": ";".join(f"{state}:{prob:.8f}" for state, prob in ranked),
            "technology": technology,
            "dominant_repeat": fractions[0][0] if fractions else "",
            "secondary_repeat": fractions[1][0] if len(fractions) > 1 else "",
            "dominant_fraction": round(fractions[0][1], 8) if fractions else "",
            "secondary_fraction": round(fractions[1][1], 8) if len(fractions) > 1 else "",
            "allele_fractions": ";".join(f"{state}:{fraction:.8f}" for state, fraction in fractions),
            "inference_method": "molecule_likelihood_em" if sample_mode == "metagenome" else "molecule_likelihood",
            "probability_kind": "single_allele_model_support",
            "em_iterations": iterations, "em_converged": converged,
            "consensus_status": consensus_status,
            **({"estimation_method": "consensus_supported",
                "estimation_reason": "Supported consensus resolves competing molecule lengths; consensus adds no independent depth.",
                "estimate_support": len(observations)} if consensus_status == "resolved_ambiguity" else {}),
            **({"estimation_method": "paired_fragment_length",
                "estimation_reason": "Opposite-flank pairs scored using an independently calibrated insert-size distribution; repeat length is not directly measured.",
                "estimate_support": len(observations)} if pair_only else {}),
        })
    return calls, molecule_calls


def estimate_unresolved_calls(calls, evidence, contexts, temperature=8.0):
    """Fill unresolved calls with labelled estimates, without creating measurements.

    Prefer conflicting direct observations, then repeat-spanning CIGAR projections,
    then reference-conditioned alignment support. Average per-molecule support so
    duplicate mappings cannot multiply votes or manufacture depth confidence.
    """
    known = {context.candidate_id: context for context in contexts}
    by_locus = defaultdict(list)
    for row in evidence:
        by_locus[row.locus_id].append(row)
    for call in calls:
        if call["repeat_count"] != "" or call["status"] == "not_found":
            continue
        rows = by_locus[call["locus"]]
        support = call["molecule_support"]
        if call["best_candidate_repeat"] != "":
            # The strict caller already ranked measured but ambiguous evidence.
            method = "ambiguous_measurements"
            reason = "Competing molecule measurements; reporting the highest-support allele."
        else:
            measured = [row for row in rows if row.measured_repeat_count is not None]
            mapped = [row for row in rows if row.candidate_id in known
                      and row.repeat_count is not None and row.alignment_identity > 0
                      and math.isfinite(row.alignment_score)
                      and row.metadata.get("aligned_flank_bases", 0) >= max(1, min(21, max(
                          known[row.candidate_id].repeat_start,
                          len(known[row.candidate_id].sequence) - known[row.candidate_id].repeat_end,
                      )))]
            spanning = [row for row in mapped if row.metadata.get("spanning_repeat_indel_delta") is not None]
            selected = measured or spanning or mapped
            if not selected:
                continue
            method, reason = (
                ("conflicting_measurements", "Conflicting measurements within molecules; no unique direct count.")
                if measured else
                ("alignment_repeat_span", "Reference repeat count adjusted by repeat-region CIGAR indels; no direct product measurement.")
                if spanning else
                ("reference_conditioned", "Partial mappings support a reference allele; repeat length is not directly resolved and alternatives outside the catalog are untested.")
            )
            molecule_scores, weights = defaultdict(dict), {}
            for row in selected:
                count = row.measured_repeat_count if measured else row.repeat_count
                if spanning and not measured:
                    count += float(row.metadata["spanning_repeat_indel_delta"])
                if not math.isfinite(count) or count < 0:
                    continue
                count = int(count) if float(count).is_integer() else count
                score = 0.0 if measured else row.alignment_score / max(temperature, 1e-6)
                scores = molecule_scores[row.molecule_id]
                scores[count] = max(scores.get(count, -math.inf), score)
                weight = 1.0 if measured else max(0.01, min(1.0, row.alignment_identity * row.query_coverage))
                weights[row.molecule_id] = max(weights.get(row.molecule_id, 0), weight)
            if not molecule_scores:
                continue
            totals = defaultdict(float)
            for molecule, scores in molecule_scores.items():
                for count, probability in _softmax(scores).items():
                    totals[count] += weights[molecule] * probability
            total_weight = sum(weights.values())
            ranked = sorted(((count, value / total_weight) for count, value in totals.items()),
                            key=lambda item: (-item[1], item[0]))
            best = ranked[0]
            second = ranked[1][1] if len(ranked) > 1 else 0.0
            support = len(molecule_scores)
            if len(ranked) == 1:
                reason += " Only one allele is supported by the available evidence; model support does not establish accuracy."
            elif math.isclose(best[1], second):
                reason += " Top alleles are tied; the smaller count is shown deterministically."
            call.update(
                best_candidate_repeat=best[0], best_probability=round(best[1], 8),
                second_best_probability=round(second, 8), margin=round(best[1] - second, 8),
                candidate_distribution=";".join(f"{count}:{probability:.8f}" for count, probability in ranked),
            )
        call.update(repeat_count=call["best_candidate_repeat"], status="estimated", confidence="low",
                    estimation_method=method, estimation_reason=reason, estimate_support=support)


ESTIMATION_FIELDS = ["estimation_method", "estimation_reason", "estimate_support"]


COMMON_LOCUS_CALL_FIELDS = [
    "sample", "locus", "repeat_count", "status", "confidence",
    "best_probability", "second_best_probability", "margin", "molecule_support", "presence_support",
    "direct_product_support", "full_span_support", "junction_support",
    "best_candidate_repeat", "candidate_distribution", "technology",
    "dominant_repeat", "secondary_repeat", "dominant_fraction", "secondary_fraction",
] + ESTIMATION_FIELDS
COMMON_LOCUS_CALL_FIELDS += [
    "repeat_count_raw", "inference_method", "probability_kind", "allele_fractions",
    "em_iterations", "em_converged", "consensus_status",
]
