#!/usr/bin/env python3
"""Compare exact repeat fitting with exhaustive versus collapsed native scoring.

Synthetic flanking reads exercise expansion to the repeat ceiling. The benchmark
excludes recruitment, QC and output writing, and checks every likelihood value,
posterior, component fraction and call field for exact equality.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import random
import sys
import time

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import parasail
from mlvamaps.models import Locus, ReadPair, ReadRecord
import mlvamaps.repeat_likelihood as model


def positive_int(value):
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError('must be positive')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--molecules', type=positive_int, default=500)
    parser.add_argument('--unit', type=positive_int, default=60)
    parser.add_argument('--maximum', type=positive_int, default=100)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rng = random.Random(271)
    def dna(n):
        return ''.join(rng.choices('ACGT', k=n))
    template = model.panel_template(Locus('L', forward_primer=dna(20), reverse_primer=dna(20),
        left_flank_sequence=dna(80), right_flank_sequence=dna(80), repeat_motif=dna(args.unit),
        repeat_unit_length_bp=args.unit, expected_max_repeats=min(80, args.maximum)))
    evidence = []
    for i in range(args.molecules):
        sequence = list(template.left + (template.motif * math.ceil(50 / args.unit))[:50])
        for position in rng.sample(range(10, 90), 2):
            sequence[position] = next(base for base in 'ACGT' if base != sequence[position])
        read = ''.join(sequence)
        item = model.classify_pair(ReadPair(str(i), ReadRecord(str(i), read, 'I'*len(read))), template)
        if item and 'F' in item.classes:
            evidence.append(item)
    if not evidence:
        raise RuntimeError('synthetic workload produced no flanking evidence')

    targets = {}
    def exhaustive(read, template, states, previous):
        scores = list(previous)
        profile = model._read_profile(read)
        for count in states[len(previous):]:
            if count not in targets:
                targets[count] = template.sequence(float(count))
            scores.append(parasail.sw_striped_profile_32(profile, targets[count], 5, 1).score)
        return np.asarray(scores, dtype=float)

    optimized = model._repeat_alignment_scores
    fits, results = [], []
    try:
        for name, score in [('exhaustive', exhaustive), ('optimized', optimized)]:
            model._repeat_alignment_scores = score
            model._read_profile.cache_clear()
            model._score_profile.cache_clear()
            started = time.perf_counter()
            fit = model.infer_repeat(evidence, template, maximum=args.maximum)
            elapsed = time.perf_counter() - started
            fits.append(fit)
            row = {'scoring': name, 'seconds': elapsed, 'molecules': len(evidence),
                   'candidate_states': len(fit.states), 'best_repeat': fit.best,
                   'identifiable': fit.identifiable,
                   'likelihood_sha256': hashlib.sha256(fit.molecule_likelihoods.tobytes()).hexdigest()}
            results.append(row)
            print(json.dumps(row), flush=True)
    finally:
        model._repeat_alignment_scores = optimized
    for name, expected in vars(fits[0]).items():
        observed = getattr(fits[1], name)
        if isinstance(expected, np.ndarray):
            np.testing.assert_array_equal(observed, expected)
        elif observed != expected:
            raise AssertionError(f'repeat inference differs in {name}')
    output = {'synthetic': True, 'unit': args.unit, 'exact_equivalence': True,
              'speedup': results[0]['seconds']/results[1]['seconds'], 'benchmarks': results}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + '\n')


if __name__ == '__main__':
    main()
