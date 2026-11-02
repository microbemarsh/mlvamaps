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
    maximum_candidates: int = 2048
    maximum_likelihood_cells: int = 262144
    maximum_iterations: int = 200
    maximum_refinement_reads: int = 100
    maximum_sequence_length: int = 10000
    maximum_alignment_cells: int = 2000000

    def __post_init__(self):
        for name in ("maximum_candidates", "maximum_likelihood_cells", "maximum_iterations",
                     "maximum_refinement_reads", "maximum_sequence_length", "maximum_alignment_cells"):
            if getattr(self, name) < 1:
                raise ValueError(f"{name} must be positive")


def _softmax(scores: dict[int | float, float]) -> dict[int | float, float]:
    if not scores:
        return {}
    maximum = max(scores.values())
    weights = {state: math.exp(max(-700.0, score - maximum)) for state, score in scores.items()}
    total = sum(weights.values())
    return {state: weight / total for state, weight in weights.items()}


def _tier_weight(evidence: CandidateEvidence) -> float:
    if "ambiguous_locus" in evidence.quality_flags:
        return 0.0
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
    if any("ambiguous_locus" in row.quality_flags for row in rows):
        return {}, None, 0.0
    direct = [row.measured_repeat_count for row in rows if row.measured_repeat_count is not None]
    if len(set(direct)) > 1 or (require_measurement and not direct):
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
        envelope = (min(states), max(states)) if states else ("", "")
        budget = max(1, min(thresholds.maximum_candidates,
                            thresholds.maximum_likelihood_cells // max(1, len(locus_items))))
        limited = len(states) > budget
        if limited:
            votes = Counter(value for rows in locus_items.values()
                            for value in {row.measured_repeat_count for row in rows if row.measured_repeat_count is not None})
            scores = defaultdict(float)
            for rows in locus_items.values():
                for row in rows:
                    scores[row.repeat_count] = max(scores[row.repeat_count], row.alignment_score)
            states = sorted(states, key=lambda state: (-votes[state], -scores[state], state))[:budget]
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
            if limited and status == "called":
                status = "estimated"

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
            "repeat_count": best[0],
            "computational_fallback": "candidate_grid_limit" if limited else "",
            "repeat_interval_low": envelope[0], "repeat_interval_high": envelope[1],
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
        derived = set()
        for molecule, rows in molecules.items():
            if any("ambiguous_locus" in row.quality_flags for row in rows):
                continue
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
            if all(row.evidence_tier == "microassembly" for row in selected):
                derived.add(molecule)
            else:
                molecule_calls[(locus.locus_id, molecule)] = selected[0].measured_repeat_count

        counts = Counter(value for molecule, value in observations.items() if molecule not in derived)
        for value in {observations[molecule] for molecule in derived}:
            counts.setdefault(value, 1)  # One reconstructed product is one likelihood, not N new reads.
        measured_support = len(observations)
        # Library geometry is an optional fallback, never a replacement for
        # directly measured products and never classification ground truth.
        if not observations and not any(row.measured_repeat_count is not None
                                        for rows in molecules.values() for row in rows):
            for molecule, rows in molecules.items():
                if any("ambiguous_locus" in row.quality_flags for row in rows):
                    continue
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
        limited, possible_secondary = "", False
        interval = ("", "")
        log_gap = ""
        if counts:
            interval = (max(0, min(raw - 3 * sigma for raw, sigma in counts)),
                        max(raw + 3 * sigma for raw, sigma in counts))
            unit = max(repeat_unit_length(locus), 1)
            states = set()
            for raw, sigma in sorted(counts):
                # Local base-pair grids retain novel/fractional alleles without
                # expanding every allele between distant observed products.
                # Pair uncertainty can span many repeats; half-unit steps keep
                # that fallback bounded without pretending to base-pair precision.
                step = max(1 / unit, .5, 3 * sigma / 200) if pair_only else 1 / unit
                width = 3 * sigma / step
                if not math.isfinite(width) or width > (thresholds.maximum_candidates - 1) // 2:
                    limited = "candidate_limit"
                    break
                radius = max(1, math.ceil(width))
                center = max(0, raw)
                states.update(normalize_allele(center + delta * step)
                              for delta in range(-radius, radius + 1)
                              if center + delta * step >= 0)
                if len(states) > thresholds.maximum_candidates:
                    limited = "candidate_limit"
                    break
            consensus = consensus_measurements.get(locus.locus_id, {})
            candidate = (consensus.get("raw_repeat_count")
                         if consensus.get("supporting_reads", 0) >= thresholds.minimum_molecules else None)
            if candidate not in (None, "") and math.isfinite(float(candidate)) and float(candidate) >= 0:
                states.add(normalize_allele(candidate))
                consensus_status = "candidate_only"
            states = sorted(states)
            if len(states) * len(counts) > thresholds.maximum_likelihood_cells:
                limited = "likelihood_cell_limit"
            if limited:
                states = []
                # Retain a bounded set of observed modes plus the Gaussian
                # sufficient-statistic estimate; score one column at a time.
                budget = max(1, min(thresholds.maximum_candidates,
                                    thresholds.maximum_likelihood_cells // len(counts)))
                centers = Counter()
                for (value, uncertainty), number in counts.items():
                    centers[max(0, value)] += number
                precision = sum(number / sigma**2 for (_, sigma), number in counts.items())
                mean = max(0, sum(number * value / sigma**2 for (value, sigma), number in counts.items()) / precision)
                candidates = [value for value, _ in sorted(centers.items(), key=lambda item: (-item[1], item[0]))]
                candidates.insert(min(1, len(candidates)), mean)
                states = list(dict.fromkeys(candidates))[:budget]
                scores = {}
                raw = np.asarray([item[0] for item in counts])
                sigma = np.asarray([item[1] for item in counts])
                multiplicities = np.asarray(list(counts.values()), dtype=float)
                for state in states:
                    delta = np.abs(raw - state)
                    core = -.5 * (delta / sigma)**2 - np.log(sigma * math.sqrt(2 * math.pi))
                    scale = np.maximum(1., sigma * 4)
                    tail = -delta / scale - np.log(2 * scale)
                    scores[state] = float(multiplicities @ np.logaddexp(math.log(.99) + core, math.log(.01) + tail))
                posterior = _softmax(scores)
                raw_best = max(posterior, key=lambda state: (posterior[state], -state))
                best_count = legacy_round_repeat_count(raw_best, round_tolerance)
                grouped_support = defaultdict(float)
                for state, probability in posterior.items():
                    grouped_support[legacy_round_repeat_count(state, round_tolerance)] += probability
                ranked = sorted(grouped_support.items(), key=lambda item: (-item[1], item[0]))
                status = "estimated"
        if counts and not limited:
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
            ordered_joint = sorted(joint, reverse=True)
            log_gap = float(ordered_joint[0] - ordered_joint[1]) if len(ordered_joint) > 1 else ""
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
                    log_likelihood, multiplicities, max_iterations=thresholds.maximum_iterations, tolerance=1e-6,
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
                if sample_mode == "isolate" and not pair_only:
                    modes = Counter()
                    uncertainties = defaultdict(list)
                    for (raw_value, uncertainty), number in counts.items():
                        allele = legacy_round_repeat_count(raw_value, round_tolerance)
                        modes[allele] += number
                        uncertainties[allele].append(uncertainty)
                    ordered_modes = sorted(modes, key=lambda allele: (-modes[allele], allele))
                    if len(ordered_modes) > 1:
                        first, second = ordered_modes[:2]
                        possible_secondary = (
                            modes[second] >= thresholds.mixture_min_molecules
                            and modes[second] / len(observations) >= thresholds.mixture_min_fraction
                            and abs(first - second) > max(.5, 3 * math.hypot(
                                max(uncertainties[first]), max(uncertainties[second])))
                        )
                        if possible_secondary:
                            fractions = [(allele, modes[allele] / len(observations)) for allele in (first, second)]
            probability = dict(ranked).get(best_count, 0.0)
            alternative = max((value for state, value in ranked if state != best_count), default=0.0)
            if mixture:
                status = "mixed"
            elif possible_secondary or unresolved_mixture or probability < thresholds.minimum_probability or probability - alternative < thresholds.minimum_margin:
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
            "repeat_count": best_count,
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
            "inference_method": ("bounded_measurement_likelihood" if limited else
                                 "molecule_likelihood_em" if sample_mode == "metagenome" else "molecule_likelihood"),
            "probability_kind": "single_allele_model_support",
            "em_iterations": iterations, "em_converged": converged,
            "consensus_status": consensus_status,
            "computational_fallback": limited,
            "possible_secondary": possible_secondary,
            "repeat_interval_low": interval[0], "repeat_interval_high": interval[1],
            "log_likelihood_gap": log_gap,
            "inference_tier": 1 if len(counts) == 1 and not pair_only else 3 if counts else 0,
            "derived_molecule_support": len(derived),
            "repeat_lower_bound": max((row.metadata.get("repeat_lower_bound") or 0 for row in rows), default=0) or "",
            **({"estimation_method": "consensus_supported",
                "estimation_reason": "Supported consensus resolves competing molecule lengths; consensus adds no independent depth.",
                "estimate_support": len(observations)} if consensus_status == "resolved_ambiguity" else {}),
            **({"estimation_method": "paired_fragment_length",
                "estimation_reason": "Opposite-flank pairs scored using a flank-calibrated or user-supplied insert-size distribution; repeat length is not directly measured and mapping selection bias is not calibrated.",
                "estimate_support": len(observations)} if pair_only else {}),
        })
    return calls, molecule_calls


def estimate_unresolved_calls(calls, evidence, contexts, temperature=8.0, thresholds=None):
    """Fill unresolved calls with labelled estimates, without creating measurements.

    Prefer conflicting direct observations, then repeat-spanning CIGAR projections,
    then reference-conditioned alignment support. Combine molecule likelihoods so
    duplicate mappings cannot multiply votes or manufacture depth confidence.
    """
    thresholds = thresholds or InferenceThresholds()
    known = {context.candidate_id: context for context in contexts}
    by_locus = defaultdict(list)
    for row in evidence:
        by_locus[row.locus_id].append(row)
    for call in calls:
        if (call["repeat_count"] != "" or call["status"] == "not_found" or call.get("computational_fallback")
                or call.get("possible_secondary")):
            continue
        ambiguous_molecules = {row.molecule_id for row in by_locus[call["locus"]] if "ambiguous_locus" in row.quality_flags}
        rows = [row for row in by_locus[call["locus"]] if row.molecule_id not in ambiguous_molecules]
        support = call["molecule_support"]
        if call["best_candidate_repeat"] != "":
            # The strict caller already ranked measured but ambiguous evidence.
            method = "ambiguous_measurements"
            call.update(estimation_method="ambiguous_measurements",
                        estimation_reason="Competing molecule measurements do not identify a unique count; retain the candidate distribution and interval.",
                        estimate_support=support)
            continue
        else:
            measured = [row for row in rows if row.measured_repeat_count is not None]
            paired = [row for row in rows if "pair_repeat_count" in row.metadata
                      and "pair_repeat_sigma" in row.metadata]
            best_pair_scores = {}
            for row in paired:
                best_pair_scores[row.molecule_id] = max(best_pair_scores.get(row.molecule_id, -math.inf),
                                                       row.metadata.get("pair_alignment_score", 0))
            paired = [row for row in paired if row.metadata.get("pair_alignment_score", 0) == best_pair_scores[row.molecule_id]]
            mapped = [row for row in rows if row.candidate_id in known
                      and row.repeat_count is not None and row.alignment_identity > 0
                      and math.isfinite(row.alignment_score)
                      and row.metadata.get("aligned_flank_bases", 0) >= max(1, min(21, max(
                          known[row.candidate_id].repeat_start,
                          len(known[row.candidate_id].sequence) - known[row.candidate_id].repeat_end,
                      )))]
            spanning = [row for row in mapped if row.metadata.get("spanning_repeat_indel_delta") is not None]
            if call.get("repeat_lower_bound") and not measured and not paired and not spanning:
                continue  # The finalizer conditions the locus prior on this bound.
            selected = measured or paired or spanning or mapped
            if not selected:
                continue
            method, reason = (
                ("conflicting_measurements", "Conflicting measurements within molecules; no unique direct count.")
                if measured else
                ("paired_fragment_length", "Competing calibrated paired-fragment lengths; alternatives within a molecule share one likelihood and do not establish a direct measurement.")
                if paired else
                ("alignment_repeat_span", "Reference repeat count adjusted by repeat-region CIGAR indels; no direct product measurement.")
                if spanning else
                ("reference_conditioned", "Partial mappings support a reference allele; repeat length is not directly resolved and alternatives outside the catalog are untested.")
            )
            molecule_scores, weights = defaultdict(dict), {}
            alternatives = defaultdict(set)
            for row in selected:
                count = row.measured_repeat_count if measured else max(0, row.metadata["pair_repeat_count"]) if paired else row.repeat_count
                if spanning and not measured and not paired:
                    count += float(row.metadata["spanning_repeat_indel_delta"])
                if not math.isfinite(count) or count < 0:
                    continue
                count = int(count) if float(count).is_integer() else count
                score = 0.0 if measured or paired else row.alignment_score / max(temperature, 1e-6)
                scores = molecule_scores[row.molecule_id]
                scores[count] = max(scores.get(count, -math.inf), score)
                if measured or paired:
                    raw = (row.raw_repeat_count if row.raw_repeat_count is not None else count) if measured else row.metadata["pair_repeat_count"]
                    sigma = (row.repeat_sigma if row.repeat_sigma is not None else .08) if measured else row.metadata["pair_repeat_sigma"]
                    if not math.isfinite(raw) or (measured and raw < 0) or not math.isfinite(sigma) or sigma <= 0:
                        raise ValueError("Invalid conflicting molecule measurement or uncertainty")
                    alternatives[row.molecule_id].add((raw, sigma))
                weight = 1.0 if measured or paired else max(0.01, min(1.0, row.alignment_identity * row.query_coverage))
                weights[row.molecule_id] = max(weights.get(row.molecule_id, 0), weight)
            if not molecule_scores:
                continue
            states = sorted({state for scores in molecule_scores.values() for state in scores})
            work = sum(len(values) for values in (alternatives or molecule_scores).values())
            budget = max(1, min(thresholds.maximum_candidates, thresholds.maximum_likelihood_cells // max(1, work)))
            if len(states) > budget:
                votes = Counter(state for scores in molecule_scores.values() for state in scores)
                states = sorted(states, key=lambda state: (-votes[state], state))[:budget]
                call["computational_fallback"] = "fallback_likelihood_grid_limit"
                reason += " Bounded candidate grid; omitted alleles remain in the evidence table."
            totals = dict.fromkeys(states, 0.)
            for molecule, scores in molecule_scores.items():
                for state in states:
                    if measured or paired:
                        values = sorted(alternatives[molecule])
                        likelihoods = [-.5 * ((state - raw) / sigma)**2 - math.log(sigma)
                                       for raw, sigma in values]
                        score = float(np.logaddexp.reduce(likelihoods)) - math.log(len(values))
                    else:
                        score = scores.get(state, min(scores.values()) - 20.) - max(scores.values())
                    totals[state] += weights[molecule] * score
            ranked = sorted(_softmax(totals).items(), key=lambda item: (-item[1], item[0]))
            best = ranked[0]
            second = ranked[1][1] if len(ranked) > 1 else 0.0
            support = len(molecule_scores)
            if alternatives:
                call["repeat_interval_low"] = max(0, min(raw - 3 * sigma for values in alternatives.values() for raw, sigma in values))
                call["repeat_interval_high"] = max(raw + 3 * sigma for values in alternatives.values() for raw, sigma in values)
            if len(ranked) == 1:
                reason += " Only one allele is supported by the available evidence; model support does not establish accuracy."
            elif math.isclose(best[1], second):
                reason += " Top alleles are tied; the smaller count is shown deterministically."
            call.update(
                best_candidate_repeat=best[0], best_probability=round(best[1], 8),
                second_best_probability=round(second, 8), margin=round(best[1] - second, 8),
                candidate_distribution=";".join(f"{count}:{probability:.8f}" for count, probability in ranked),
            )
        call.update(estimation_method=method, estimation_reason=reason, estimate_support=support)
        call.update(repeat_count=call["best_candidate_repeat"],
                    repeat_count_raw=call["best_candidate_repeat"],
                    status="estimated", confidence="low")


def locus_repeat_prior(locus, contexts, maximum_candidates=2048):
    """Equal mass on distinct observed alleles, then calibrated nominal/range.

    Database multiplicity is deliberately ignored. The default 0..100 search
    range is not a biological prior. Large ranges use a bounded uniform grid.
    """
    states = sorted({context.repeat_count for context in contexts
                     if context.locus_id == locus.locus_id
                     and context.observed_or_synthetic == "observed"
                     and context.repeat_count is not None
                     and math.isfinite(context.repeat_count) and context.repeat_count >= 0})
    source = "observed_reference_prior"
    if not states and locus.expected_product_size_bp > 0 and repeat_unit_length(locus) > 0:
        # Zero is also the legacy missing nominal default; require an observed
        # zero allele or an explicit zero range to make it an eligible prior.
        if math.isfinite(locus.nominal_repeat_units) and locus.nominal_repeat_units > 0:
            states = [locus.nominal_repeat_units]
            source = "calibrated_panel_prior"
    low, high = locus.expected_min_repeats, locus.expected_max_repeats
    limited = ""
    if not states and (low, high) != (0, 100) and 0 <= low <= high:
        size = high - low + 1
        count = min(size, maximum_candidates)
        states = sorted({low + round(i * (high - low) / max(1, count - 1)) for i in range(count)})
        source = "panel_range_prior"
        limited = "prior_grid_limit" if size > count else ""
    if not states:
        return {}, "", ""
    if len(states) > maximum_candidates:
        states = [states[round(i * (len(states) - 1) / max(1, maximum_candidates - 1))]
                  for i in range(maximum_candidates)]
        limited = "prior_grid_limit"
    return {state: 1 / len(states) for state in states}, source, limited


def finalize_repeat_count_calls(calls, evidence, loci, contexts, thresholds=None):
    """Complete FASTQ predictions after reconstruction, without creating reads."""
    thresholds = thresholds or InferenceThresholds()
    by_locus = defaultdict(list)
    contexts_by_locus = defaultdict(list)
    for row in evidence:
        by_locus[row.locus_id].append(row)
    for context in contexts:
        contexts_by_locus[context.locus_id].append(context)
    locus_by_id = {locus.locus_id: locus for locus in loci}
    if len(locus_by_id) != len(loci) or Counter(call["locus"] for call in calls) != Counter(locus_by_id.keys()):
        raise ValueError("FASTQ completion requires exactly one call per unique configured locus")
    for call in calls:
        locus = locus_by_id[call["locus"]]
        rows = by_locus[locus.locus_id]
        ambiguous = {row.molecule_id for row in rows if "ambiguous_locus" in row.quality_flags}
        usable = [row for row in rows if row.molecule_id not in ambiguous]
        refinements = sorted({str(row.metadata["refinement_fallback"]) for row in usable
                              if row.metadata.get("refinement_fallback")})
        if refinements and not call.get("computational_fallback"):
            call["computational_fallback"] = "motif_refinement:" + ";".join(refinements)
        source = call.get("estimation_method") or (
            "reconstructed_product" if call.get("inference_method") == "mapped_consensus_assembly_pcr" or call.get("derived_molecule_support")
            else "joint_measurement_likelihood" if any(row.measured_repeat_count is not None for row in usable)
            else "reference_conditioned")
        prior_only = bool(call.get("prior_only", False))
        if call["repeat_count"] in ("", None):
            prior, source, limited = locus_repeat_prior(locus, contexts_by_locus[locus.locus_id], thresholds.maximum_candidates)
            bounds = {}
            for row in usable:
                bound = row.metadata.get("repeat_lower_bound")
                if bound is not None and math.isfinite(bound) and bound > 0:
                    bounds[row.molecule_id] = max(bound, bounds.get(row.molecule_id, 0))
            if not prior and not bounds:
                raise ValueError(f"Locus {locus.locus_id}: no informative repeat evidence or eligible locus-specific prior; "
                                 "provide observed reference alleles, calibrated nominal units, or an explicit non-default repeat range")
            reason = "Prior-only prediction; no sample-derived repeat-length evidence. Presence is reported separately."
            if bounds:
                lower = max(bounds.values())
                scale = max(1., max(prior, default=lower) - min(prior, default=lower))
                # ponytail: an exponential continuation allows novel expansions
                # beyond a finite catalog. Its tail is a sensitivity assumption,
                # not a validated biological maximum or a calibrated interval.
                states = sorted({max(lower, state) for state in prior} |
                                {lower + scale * i for i in range(min(32, thresholds.maximum_candidates))})
                if len(states) > thresholds.maximum_candidates:
                    states = states[:thresholds.maximum_candidates]
                    limited = "partial_likelihood_grid_limit"
                prior = _softmax({state: -min((abs(state - center) for center in prior), default=state - lower) / scale
                                  for state in states})
                source = "partial_repeat_likelihood"
                reason = ("Joint anchored lower-bound constraints condition the locus prior with an exponential expansion tail. "
                          "No upper repeat boundary was measured; the reported interval is a bounded candidate envelope.")
                call.update(status="estimated", repeat_lower_bound=lower)
            best = min(prior, key=lambda state: (-prior[state], state))
            call.update(repeat_count=best, repeat_count_raw=best, best_candidate_repeat=best,
                        best_probability=prior[best], second_best_probability=sorted(prior.values(), reverse=True)[1] if len(prior) > 1 else 0.,
                        candidate_distribution=";".join(f"{state}:{prob:.8f}" for state, prob in prior.items()),
                        repeat_interval_low=min(prior), repeat_interval_high=max(prior),
                        computational_fallback=call.get("computational_fallback") or limited,
                        confidence="low" if bounds else "unresolved", estimation_method=source,
                        estimation_reason=reason, estimate_support=len(bounds))
            prior_only = not bounds
            if prior_only and limited:
                if source == "panel_range_prior":
                    call.update(repeat_interval_low=locus.expected_min_repeats, repeat_interval_high=locus.expected_max_repeats)
                elif source == "observed_reference_prior":
                    observed = [context.repeat_count for context in contexts_by_locus[locus.locus_id]
                                if context.observed_or_synthetic == "observed" and context.repeat_count is not None
                                and math.isfinite(context.repeat_count) and context.repeat_count >= 0]
                    call.update(repeat_interval_low=min(observed), repeat_interval_high=max(observed))
        value = call["repeat_count"]
        if not math.isfinite(float(value)) or float(value) < 0:
            raise ValueError(f"Locus {locus.locus_id}: invalid final repeat count {value!r}")
        if call.get("repeat_count_raw") in (None, ""):
            call["repeat_count_raw"] = value
        distribution = [(float(state), float(probability)) for state, probability in
                        (entry.split(":") for entry in call.get("candidate_distribution", "").split(";") if entry)]
        alternative = sorted(((state, probability) for state, probability in distribution if state != float(value)),
                             key=lambda item: (-item[1], item[0]))
        support = 0 if prior_only else int(call.get("estimate_support") or call["molecule_support"])
        model_support = float(call.get("best_probability", 0))
        # ponytail: an uncalibrated evidence-strength score, not accuracy.
        # Replace with empirical reliability calibration when truth data exist.
        weight = {"joint_measurement_likelihood": 1., "reconstructed_product": .5,
                  "conflicting_measurements": .25, "paired_fragment_length": .35,
                  "alignment_repeat_span": .4, "reference_conditioned": .1,
                  "partial_repeat_likelihood": .1}.get(source, .25)
        sigma = max((row.repeat_sigma or .08 for row in usable if row.measured_repeat_count is not None), default=.08)
        confidence = model_support * support / (support + max(1, thresholds.minimum_molecules)) * weight / (1. + sigma)
        if any(row.quality_flags for row in usable if row.measured_repeat_count is not None):
            confidence *= .5
        if call["status"] in {"ambiguous", "mixed"} or call.get("computational_fallback"):
            confidence *= .25
        if call.get("repeat_interval_low") in (None, ""):
            call["repeat_interval_low"] = min((state for state, _ in distribution), default=value)
        if call.get("repeat_interval_high") in (None, ""):
            call["repeat_interval_high"] = max((state for state, _ in distribution), default=value)
        call.update(estimate_source=source, prior_only=prior_only, independent_support=support,
                    model_support=model_support, data_evidence_confidence=round(confidence, 8),
                    second_best_repeat_count=alternative[0][0] if alternative else "",
                    second_best_support=alternative[0][1] if alternative else 0.,
                    confidence_reason=call.get("estimation_reason") or
                    "Sample molecule likelihood; model support and evidence confidence are uncalibrated.")


CONFIDENCE_FIELDS = ["estimate_source", "genotype_status", "model_support", "data_evidence_confidence",
                     "second_best_support", "repeat_interval_low", "repeat_interval_high",
                     "independent_support", "prior_only", "computational_fallback", "confidence_reason"]
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
    "computational_fallback", "possible_secondary", "repeat_interval_low", "repeat_interval_high",
    "log_likelihood_gap", "inference_tier",
    "derived_molecule_support", "repeat_lower_bound",
    "genotype_status",
]
COMMON_LOCUS_CALL_FIELDS += [field for field in CONFIDENCE_FIELDS + ["second_best_repeat_count"]
                            if field not in COMMON_LOCUS_CALL_FIELDS]
