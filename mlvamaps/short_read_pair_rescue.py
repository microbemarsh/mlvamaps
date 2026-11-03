"""Conditional repeat estimates from intact read pairs, without reconstruction.

Catalog flanks and periodic arrays are mapping hypotheses, never sample
sequences. Interrupted longer arrays remain possible: every rescue is unbounded.
"""
from collections import defaultdict
import math

from .calling import assembly_equivalent_product_allele
from .sequence import revcomp
from .short_read_fragments import _seed_matcher
from .short_read_likelihood import _log_normal_interval, _log_sum


class WeakPairCollector:
    """A separate bounded candidate pool; seed hits are not locus assignments."""

    def __init__(self, contexts, maximum_per_locus=2000, seed_length=12):
        if maximum_per_locus < 1 or seed_length != 12:
            raise ValueError("Pair recruitment needs a positive cap and exactly 12 bp seeds")
        self.maximum_per_locus, self.seed_length = maximum_per_locus, seed_length
        owners, repeats = defaultdict(set), set()
        self.backgrounds = defaultdict(set)
        self.unsupported_loci = set()
        for context in contexts:
            for flank in (context.sequence[:context.repeat_start], context.sequence[context.repeat_end:]):
                for start in range(len(flank) - seed_length + 1):
                    seed = flank[start:start + seed_length]
                    if set(seed) <= set("ACGT") and len(set(seed)) >= 3:
                        owners[seed].add(context.locus_id)
                        owners[revcomp(seed)].add(context.locus_id)
            repeats.add(context.sequence[context.repeat_start:context.repeat_end])
            motif = context.repeat_motif
            if motif and set(motif) <= set("ACGT"):
                repeats.add(motif * (2 + seed_length // len(motif)))
                self.backgrounds[context.locus_id].add((context.sequence[:context.repeat_start],
                    context.sequence[context.repeat_end:], motif))
            else:
                self.unsupported_loci.add(context.locus_id)
        matcher = _seed_matcher(owners)
        excluded = {match.group(1) for sequence in repeats for match in matcher.finditer(sequence)} if matcher else set()
        excluded |= {revcomp(seed) for seed in excluded}
        self.seeds = {seed: loci for seed, loci in owners.items() if seed not in excluded}
        self.matcher = _seed_matcher(self.seeds)
        self.fragments = defaultdict(dict)
        self.competitors = {}
        self.overflow_loci, self.conflicting_molecules = set(), set()
        self.diagnostics = dict(paired_molecules_seen=0, seed_matched_molecules=0,
                                overflow_loci=[], conflicting_molecules=0)

    def add(self, pair):
        if pair.read2 is None or self.matcher is None:
            return
        self.diagnostics["paired_molecules_seen"] += 1
        loci = {locus for read in (pair.read1, pair.read2)
                for match in self.matcher.finditer(read.sequence)
                for locus in self.seeds[match.group(1)]}
        if not loci or pair.molecule_id in self.conflicting_molecules:
            return
        self.diagnostics["seed_matched_molecules"] += 1
        previous = next((self.fragments[locus].get(pair.molecule_id)
                         for locus in self.competitors.get(pair.molecule_id, loci)
                         if pair.molecule_id in self.fragments[locus]), None)
        signature = lambda item: sorted((read.sequence, read.quality or "") for read in (item.read1, item.read2))
        if previous is not None and signature(previous) != signature(pair):
            self.conflicting_molecules.add(pair.molecule_id)
            self.diagnostics["conflicting_molecules"] = len(self.conflicting_molecules)
            return
        # A 12-mer alone is common in a large genome. Require its intact
        # read to cross an exact flank/periodic-array junction before retaining
        # the original pair; copy count and its opposite mate are not selected.
        eligible = {locus for locus in loci - self.overflow_loci
                    if _models(pair, self.backgrounds[locus], self, locus)}
        for locus in eligible:
            retained = self.fragments[locus]
            if pair.molecule_id not in retained and len(retained) >= self.maximum_per_locus:
                # ponytail: fail closed at the cap instead of subsampling a
                # competing hypothesis away. A native mapper is the next rung.
                self.overflow_loci.add(locus)
                retained.clear()
                self.diagnostics["overflow_loci"] = sorted(self.overflow_loci)
            else:
                retained.setdefault(pair.molecule_id, pair)
                self.competitors.setdefault(pair.molecule_id, set()).update(loci)


def _backgrounds(contexts, locus):
    result = set()
    for context in contexts:
        if context.locus_id != locus.locus_id:
            continue
        motif = context.repeat_motif or locus.repeat_motif
        left, right = context.sequence[:context.repeat_start], context.sequence[context.repeat_end:]
        if (motif and not set(motif + left + right) - set("ACGT")
                and left.startswith(locus.forward_primer) and right.endswith(revcomp(locus.reverse_primer))):
            result.add((left, right, motif))
    return tuple(sorted(result))


def _models(pair, backgrounds, collector, locus_id):
    """Exact primary anchors, with disjoint left- and right-primary classes."""
    result = set()
    for left, right, motif in backgrounds:
        for first, second in ((pair.read1, pair.read2), (pair.read2, pair.read1)):
            oriented = [("left", first.sequence, revcomp(second.sequence),
                         (second.quality or "5" * len(second.sequence))[::-1], left, right, motif, -1)]
            for phase in range(len(motif)):
                reverse_motif = revcomp(motif)
                shift = (len(motif) - phase) % len(motif)
                cycle = reverse_motif[shift:] + reverse_motif[:shift]
                oriented.append(("right", second.sequence, revcomp(first.sequence),
                    (first.quality or "5" * len(first.sequence))[::-1], revcomp(right), revcomp(left), cycle, phase))
            left_found = False
            for side, primary, opposite, quality, flank, other_flank, cycle, phase in oriented:
                if side == "right" and left_found:
                    continue
                if len(quality) != len(opposite) or set(primary + opposite) - set("ACGT"):
                    continue
                target = flank + (cycle * (1 + len(primary) // len(cycle)))[:len(primary)]
                start = max(0, len(flank) - len(primary) + 1)
                while True:
                    start = target.find(primary, start)
                    if start < 0 or start > len(flank) - collector.seed_length:
                        break
                    anchor = len(flank) - start
                    if any(locus_id in collector.seeds.get(primary[i:i + collector.seed_length], ())
                           for i in range(anchor - collector.seed_length + 1)):
                        result.add((side, anchor, len(primary) - anchor, opposite, quality,
                                    other_flank, cycle, phase, len(motif)))
                        left_found |= side == "left"
                    start += 1
    return tuple(sorted(result))


def _substitution_score(observed, expected, quality):
    if observed == expected:
        return 0.
    mismatch = -1
    for index, (one, two) in enumerate(zip(observed, expected)):
        if one != two:
            if mismatch >= 0:
                return -math.inf
            mismatch = index
    if mismatch < 0 or len(observed) != len(expected):
        return -math.inf
    # ponytail: one substitution and a Q30 ceiling avoid treating a weak
    # junction base as certain. A native pair-HMM is the upgrade for indels.
    error = 10 ** (-max(2, min(30, ord(quality[mismatch]) - 33)) / 10)
    return math.log(error / (3 * (1 - error)))


def _profile(model):
    """Cache sequence comparisons by periodic phase and boundary endpoint."""
    side, anchor, _partial, opposite, quality, flank, motif, _phase, _period = model
    size, period = len(opposite), len(motif)
    periodic = motif * (2 + size // period)
    pure, boundary = [], defaultdict(list)
    for phase in range(period):
        score = _substitution_score(opposite, periodic[phase:phase + size], quality)
        if math.isfinite(score):
            pure.append(((phase + anchor + size) % period, score))
        maximum = min(len(flank), size - 1, 11 if side == "right" else size - 1)
        for bases in range(1, maximum + 1):
            start = (phase + bases - size) % period
            continuation = periodic[start:start + size]
            expected = continuation[:size - bases] + flank[:bases]
            score = _substitution_score(opposite, expected, quality)
            if math.isfinite(score):
                boundary[phase].append((bases, score, score == 0 and opposite != continuation))
    return pure, boundary


def _phase_masses(fragment_weights, period):
    result = {}
    for phase in range(period):
        entries = [(length, math.exp(weight)) for length, weight in fragment_weights.items() if length % period == phase]
        lengths = [length for length, _ in entries]
        prefix, suffix = [0.], [0.]
        for _, weight in entries:
            prefix.append(prefix[-1] + weight)
        for _, weight in reversed(entries):
            suffix.append(suffix[-1] + weight)
        result[phase] = (lengths, prefix, suffix)
    return result


def _phase_mass(masses, phase, low, high):
    from bisect import bisect_left, bisect_right
    lengths, prefix, suffix = masses[phase]
    first, last = bisect_left(lengths, low), bisect_right(lengths, high)
    value = (suffix[len(lengths) - first] - suffix[len(lengths) - last]
             if first + last > len(lengths) else prefix[last] - prefix[first])
    return math.log(value) if value > 0 else -math.inf


def _model_score(model, length, mean, sigma, fragment_weights, profile=None, masses=None):
    side, anchor, partial, opposite, _quality, flank, motif, phase, period = model
    if length < partial or (phase >= 0 and length % period != phase):
        return -math.inf, False
    size = len(opposite)
    maximum = min(len(flank), size - 1, 11 if side == "right" else size - 1)
    low, high = max(anchor + size, anchor + partial), anchor + length + maximum
    if high < low:
        return -math.inf, False
    normalizer = _log_normal_interval((low - .5 - mean) / sigma, (high + .5 - mean) / sigma)
    if not math.isfinite(normalizer):
        return -math.inf, False
    pure, endpoints = profile if profile is not None else _profile(model)
    masses = masses if masses is not None else _phase_masses(fragment_weights, len(motif))
    scores = [_phase_mass(masses, residue, low, anchor + length) + score for residue, score in pure
              if low <= anchor + length]
    boundary = False
    for bases, score, informative in endpoints.get(length % len(motif), ()):
        fragment = anchor + length + bases
        if low <= fragment <= high and fragment in fragment_weights:
            scores.append(fragment_weights[fragment] + score)
            boundary |= informative and abs(fragment - mean) <= 4 * sigma
    return (_log_sum(scores) - normalizer, boundary) if scores else (-math.inf, False)


def infer_pair_repeat_rescue(loci, contexts, collector, calibration, thresholds, *,
                            sample_mode="isolate", round_tolerance=.25, target_loci=None):
    """Return call updates and audit rows, never measured/classification evidence."""
    updates, diagnostics = {}, []
    if collector is None:
        return updates, diagnostics
    try:
        mean, sigma = float(calibration["mean_bp"]), float(calibration["sd_bp"])
        valid = math.isfinite(mean) and math.isfinite(sigma) and mean > 0 and sigma > 0
    except (TypeError, KeyError, ValueError):
        valid = False
    if sample_mode != "isolate" or not valid:
        status = "unsupported_sample_mode" if sample_mode != "isolate" else "no_fragment_calibration"
        return {}, [dict(locus=locus.locus_id, status=status) for locus in loci]
    if not math.isfinite(mean + 8 * sigma) or not math.isfinite(mean - 8 * sigma):
        return {}, [dict(locus=locus.locus_id, status="invalid_fragment_calibration") for locus in loci]
    low_fragment, high_fragment = max(1, math.floor(mean - 8 * sigma)), math.ceil(mean + 8 * sigma)
    length_high = high_fragment
    if length_high + 1 > thresholds.maximum_candidates:
        return {}, [dict(locus=locus.locus_id, status="pair_likelihood_budget_exceeded") for locus in loci]
    weights = {length: _log_normal_interval((length - .5 - mean) / sigma, (length + .5 - mean) / sigma)
               for length in range(low_fragment, high_fragment + 1)}
    staged, compatible, blocked = {}, defaultdict(set), set(collector.conflicting_molecules)
    profiles, mass_tables = {}, {}
    # Competing loci remain even if their own call is already resolved.
    relevant_ids = {identifier for name, values in collector.fragments.items()
                    if target_loci is None or name in target_loci for identifier in values}
    for locus in loci:
        name = locus.locus_id
        diagnostic = dict(locus=name, status="no_eligible_pairs", recruited_pairs=len(collector.fragments.get(name, {})),
                          calibration_source=calibration.get("source", ""), upper_unbounded=True,
                          finite_length_support=False, whole_library=calibration.get("whole_library", False))
        diagnostics.append(diagnostic)
        if name in collector.overflow_loci:
            diagnostic["status"] = "candidate_pool_overflow"
            continue
        backgrounds = _backgrounds(contexts, locus)
        if not backgrounds or len({len(left) + len(right) for left, right, _ in backgrounds}) != 1:
            diagnostic["status"] = "unsupported_flank_calibration"
            blocked.update(collector.fragments.get(name, {}))
            continue
        models = {identifier: _models(pair, backgrounds, collector, name)
                  for identifier, pair in collector.fragments.get(name, {}).items()
                  if identifier in relevant_ids and identifier not in blocked
                  and not collector.competitors[identifier] & (collector.overflow_loci | collector.unsupported_loci)}
        models = {identifier: values for identifier, values in models.items() if values}
        distinct = {model for values in models.values() for model in values}
        profile_cells = sum(len(model[6]) * (1 + min(len(model[5]), len(model[3]) - 1,
                            11 if model[0] == "right" else len(model[3]) - 1)) * len(model[3])
                            for model in distinct if model not in profiles)
        # ponytail: bound Python sequence comparisons before creating profiles,
        # including patterns that eventually have no finite placements.
        if profile_cells > 16 * thresholds.maximum_likelihood_cells:
            diagnostic.update(status="pair_profile_budget_exceeded", requested_profile_cells=profile_cells)
            blocked.update(models)
            continue
        for model in distinct:
            if model not in profiles:
                profiles[model] = _profile(model)
            period = len(model[6])
            if period not in mass_tables:
                mass_tables[period] = _phase_masses(weights, period)
        cells = sum((length_high + 1) * max(1, len(profiles[model][0])
                    + max(map(len, profiles[model][1].values()), default=0))
                    for alternatives in set(models.values()) for model in alternatives)
        if cells > thresholds.maximum_likelihood_cells:
            diagnostic.update(status="pair_likelihood_budget_exceeded", requested_cells=cells)
            blocked.update(models)
            continue
        scores, boundaries, pattern_scores = {}, {}, {}
        for identifier, alternatives in models.items():
            if alternatives in pattern_scores:
                scores[identifier], boundaries[identifier] = pattern_scores[alternatives]
                compatible[identifier].add(name)
                continue
            values, sides = [], []
            for length in range(length_high + 1):
                scored = [_model_score(model, length, mean, sigma, weights, profiles[model], mass_tables[len(model[6])])
                          for model in alternatives]
                values.append(_log_sum([score for score, _ in scored]) - math.log(len(alternatives)))
                sides.append({model[0] for model, (_, boundary) in zip(alternatives, scored) if boundary})
            if any(math.isfinite(value) for value in values):
                scores[identifier], boundaries[identifier] = values, sides
                pattern_scores[alternatives] = (values, sides)
                compatible[identifier].add(name)
        staged[name] = (locus, backgrounds, models, scores, boundaries, diagnostic)
    for name, (locus, backgrounds, models, scores, boundaries, diagnostic) in staged.items():
        if target_loci is not None and name not in target_loci:
            diagnostic["status"] = "competitor_only"
            continue
        usable = sorted(identifier for identifier in scores if compatible[identifier] == {name} and identifier not in blocked)
        diagnostic.update(eligible_pairs=len(usable), competing_pairs=len(scores) - len(usable))
        if len(usable) < 4:
            diagnostic["status"] = "insufficient_independent_pairs"
            continue
        joint = [sum(scores[identifier][length] for identifier in usable) for length in range(length_high + 1)]
        best_length = max(range(len(joint)), key=lambda length: (joint[length], -length))
        best_score = joint[best_length]
        if not math.isfinite(best_score):
            diagnostic["status"] = "incompatible_pair_patterns"
            continue
        left = {identifier for identifier in usable if "left" in boundaries[identifier][best_length]}
        right = {identifier for identifier in usable if "right" in boundaries[identifier][best_length]}
        diagnostic.update(left_boundary_pairs=len(left), right_boundary_pairs=len(right), best_repeat_bases=best_length)
        if len(left) < 2 or len(right) < 2 or len(left | right) < 4:
            diagnostic["status"] = "insufficient_reciprocal_junctions"
            continue
        nonrepeat = len(backgrounds[0][0]) + len(backgrounds[0][1])
        raw, best = assembly_equivalent_product_allele(locus, nonrepeat + best_length, round_tolerance)
        if raw is None:
            diagnostic["status"] = "unsupported_allele_calibration"
            continue
        masses = defaultdict(float)
        for length, score in enumerate(joint):
            if math.isfinite(score):
                allele = assembly_equivalent_product_allele(locus, nonrepeat + length, round_tolerance)[1]
                masses[allele] += math.exp(score - best_score)
        total = sum(masses.values())
        probabilities = sorted(((allele, mass / total) for allele, mass in masses.items()), key=lambda item: (-item[1], item[0]))
        support = masses[best] / total
        alternative = max((probability for allele, probability in probabilities if allele != best), default=0.)
        partial = max(min(model[2] for model in models[identifier]) for identifier in usable)
        lower = assembly_equivalent_product_allele(locus, nonrepeat + partial, 0)[0]
        reason = ("Joint read-pair placement likelihood conditional on an exact flank anchor and eligible opposite-mate "
                  "placement, with discrete insert probabilities and one-substitution alternatives. Assumes a pure "
                  "periodic array and the calibrated fragment component; longer interrupted arrays or unobserved "
                  "fragment/background components remain possible. Support is conditional and uncalibrated; "
                  "no full repeat or finite upper length was measured.")
        diagnostic.update(status="estimated", repeat_count=best, conditional_support=support,
                          selected_pairs=len(usable), assumption="pure_motif_and_calibrated_fragment_component")
        updates[name] = dict(repeat_count=best, repeat_count_raw=raw, best_candidate_repeat=best,
            status="estimated", confidence="low", inference_tier=3, best_probability=round(support, 8),
            second_best_probability=round(alternative, 8), margin=round(support - alternative, 8),
            candidate_distribution=";".join(f"{allele}:{probability:.8f}" for allele, probability in probabilities),
            candidate_grid_low=assembly_equivalent_product_allele(locus, nonrepeat, 0)[0],
            candidate_grid_high=assembly_equivalent_product_allele(locus, nonrepeat + length_high, 0)[0],
            repeat_interval_low=lower, repeat_interval_high="", repeat_lower_bound=lower,
            upper_unbounded=True, finite_length_support=False, paired_span_support=len(left | right),
            partial_bound_support=len(usable), repeat_rich_mate_support=len(usable) - len(left | right),
            length_model_status="weak_pair_estimated", length_null_fraction="",
            estimate_support=len(usable), dominant_repeat=best, dominant_fraction=1., secondary_repeat="",
            secondary_fraction="", allele_fractions=f"{best}:1.00000000",
            inference_method="conditional_joint_pair_likelihood", probability_kind="conditional_fragment_component_support",
            estimation_method="joint_pair_read_likelihood", estimation_reason=reason, computational_fallback="",
            pair_rescue_support=len(usable), pair_rescue_boundary_support=len(left | right),
            pair_rescue_calibration_source=calibration.get("source", ""), pair_rescue_status="estimated",
            _molecule_ids=tuple(usable))
    return updates, diagnostics
