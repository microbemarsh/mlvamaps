"""Database-free competition over physical repeat lengths.

Locus assignment uses the combined native flank aligner. Once a read is
anchored, a one-sided repeat observation is censored length evidence, not an
alignment preference for a shorter reference. Paired fragments compete over
lengths using the measured library distribution. Neither operation needs an
allele database or a minimum molecule count.
"""
from collections import Counter
import math

import numpy as np


def candidate_likelihood(items, template, insert, maximum, minimum_probability):
    """Haploid fragment likelihood with phase-corrected, censored boundaries.

    Compute geometry once per molecule and score all lengths as NumPy vectors.
    A boundary-only read has equal likelihood for every compatible length;
    sequencing depth cannot identify its unobserved end. Identical fragment
    geometries share computation while retaining their molecule multiplicity.
    """
    from .targeted_reconstruction import LocusRecovery

    if maximum < 1 or maximum * template.unit > 1_000_000:
        raise ValueError('invalid candidate ceiling or candidate locus exceeds sequence safety limit')
    nonrepeat = len(template.left)+len(template.right)
    lower = 0.
    offsets = Counter()
    for item in items:
        if 'discordant' in item.classes:
            continue
        context = item.template or template
        phase = len(context.left)+len(context.right)-nonrepeat
        if item.lower_bound:
            lower = max(lower, item.lower_bound+phase/template.unit)
        if insert and 'FLANK_PAIR' in item.classes and item.fragment_offset is not None:
            offsets[item.fragment_offset-phase] += 1
    upper = min(maximum, max(10, template.locus.expected_max_repeats+2, math.ceil(lower)+3))
    if insert and offsets:
        upper = min(maximum, max(upper, math.ceil((insert.mean+4*insert.sd-min(offsets))/template.unit)))
    if lower and not offsets:
        upper = maximum  # A censored observation supplies no finite upper bound.
    step = .5 if template.unit > 1 else 1.
    while True:
        states = np.arange(0, upper+step/2, step)
        joint = np.zeros(len(states))
        for offset, count in offsets.items():
            # Student-t tails tolerate fragment outliers without an absolute
            # support threshold. Library uncertainty is retained at low depth.
            fragment = offset+states*template.unit
            joint -= count*2.5*np.log1p(((fragment-insert.mean)/insert.sd)**2/4)
        if lower:
            # This is a censored observation. Repeating the same boundary must
            # not concentrate probability at the shortest compatible allele.
            joint[states < lower-.25] = -np.inf
        finite = np.isfinite(joint)
        if not finite.any():
            return LocusRecovery(method='AMBIGUOUS', limit_reached=True,
                reason='candidate_limit_below_observed_lower_bound')
        posterior = np.exp(np.maximum(joint-joint[finite].max(), -700))
        posterior[~finite] = 0
        posterior /= posterior.sum()
        if upper >= maximum or posterior[-2:].sum() < .01:
            break
        upper = min(maximum, upper*2)
    ranking = np.argsort(-posterior, kind='stable')
    best, second = map(int, ranking[:2])
    selected = ranking[:np.searchsorted(np.cumsum(posterior[ranking]), .95)+1]
    limit = bool(lower and not offsets) or upper == maximum and posterior[-2:].sum() >= .01
    implied = Counter()
    for offset, count in offsets.items():
        implied[round((insert.mean-offset)/template.unit)] += count
    modes = implied.most_common()
    mixed = bool(len(modes)>1 and modes[1][1] >= .2*sum(implied.values())
                 and abs(modes[0][0]-modes[1][0])*template.unit > 4*insert.sd)
    confidence = float(posterior[best])
    informative = bool(lower or offsets)
    tied = np.count_nonzero(np.isclose(joint, joint[best], rtol=0, atol=1e-8)) > 1
    identifiable = bool(offsets) and confidence >= minimum_probability and not limit and not mixed
    return LocusRecovery(method='MIXED' if mixed else 'INFERRED' if identifiable else 'AMBIGUOUS',
        confidence=confidence if informative else 0, states=states, log_likelihoods=joint, posterior=posterior,
        best=float(states[best]) if informative else None, second=float(states[second]) if informative else None,
        interval=((float(states[finite].min()), float(maximum)) if lower and not offsets else
                  (float(states[selected].min()), float(states[selected].max()))) if informative else None,
        identifiable=identifiable, limit_reached=limit,
        reason='no_repeat_length_information' if not informative else
               'repeat_length_lower_bound' if not offsets else
               'minimum_likelihood_tie' if tied else
               'candidate_limit_reached' if limit else 'competitive_fragment_likelihood')
