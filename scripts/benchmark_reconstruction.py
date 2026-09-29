#!/usr/bin/env python3
"""Synthetic end-to-end paired FASTQ benchmark; not a biological accuracy claim."""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, replace
import json
from pathlib import Path
import random
import resource
import sys
import time

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mlvamaps.models import Locus
from mlvamaps.sequence import revcomp
from mlvamaps.short_reads import run_short_read_call


def make_inputs(output, molecules, scenario):
    rng = random.Random(271)
    def dna(n):
        return ''.join(rng.choices('ACGT', k=n))
    output.mkdir(parents=True, exist_ok=True)
    loci, sequences = [], []
    for i in range(6):
        locus = Locus(f'L{i}', forward_primer=dna(20), reverse_primer=dna(20),
            left_flank_sequence=dna(160 if scenario in ('rescue', 'pileup') else 35),
            right_flank_sequence=dna(160 if scenario in ('rescue', 'pileup') else 35), repeat_motif='AATG',
            repeat_unit_length_bp=4, expected_max_repeats=15)
        repeat = (5+i) if scenario in ('direct', 'rescue', 'pileup') else (70+i)
        sequences.append(locus.forward_primer+locus.left_flank_sequence+'AATG'*repeat+
                         locus.right_flank_sequence+revcomp(locus.reverse_primer))
        if scenario in ('rescue', 'pileup'):
            locus = replace(locus, left_flank_sequence='', right_flank_sequence='', repeat_motif='NNNN',
                            expected_product_size_bp=len(sequences[-1]), nominal_repeat_units=repeat)
        loci.append(locus)
    panel = output/'panel.tsv'
    with panel.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(loci[0])), delimiter='\t')
        writer.writeheader()
        writer.writerows(asdict(l) for l in loci)
    contexts = [dna(80)+sequence+dna(80) for sequence in sequences] if scenario == 'pileup' else []
    first, second = output/'reads1.fq', output/'reads2.fq'
    with first.open('w') as a, second.open('w') as b:
        for i in range(molecules):
            # Include off-target background, with no per-locus input rescans.
            seq = dna(400) if i % 4 == 3 else sequences[i % len(sequences)]
            r1, r2 = seq[:150], revcomp(seq[-150:])
            if scenario == 'rescue' and i % 4 != 3:
                r1, mate_sequence = ((seq[:120], seq[70:210]),
                              (seq[-210:-70], seq[-120:]),
                              (seq[140:-140], seq[145:-135]))[(i//6) % 3]
                r2 = revcomp(mate_sequence)
            if scenario == 'pileup' and i % 4 != 3:
                context = contexts[i % len(contexts)]
                start = rng.randrange(len(context)-200+1)
                read = list(context[start:start+200])
                if i % 3:
                    position = rng.randrange(len(read))
                    read[position] = rng.choice([base for base in 'ACGT' if base != read[position]])
                r1 = ''.join(read)
                r2 = revcomp(r1)
            a.write(f'@r{i}/1\n{r1}\n+\n'+ 'I'*len(r1)+'\n')
            b.write(f'@r{i}/2\n{r2}\n+\n'+ 'I'*len(r2)+'\n')
    return first, second, panel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--molecules', type=int, default=10000)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--scenario', choices=['direct', 'unresolved', 'rescue', 'pileup'], default='direct')
    args = parser.parse_args()
    if args.molecules < 1 or args.threads < 1:
        parser.error('molecules and threads must be positive')
    first, second, panel = make_inputs(args.output, args.molecules, args.scenario)
    started = time.perf_counter()
    run_short_read_call(str(first), str(second), str(panel), str(args.output/'result'), 'bench',
                        threads=args.threads, show_progress=False)
    result = {'synthetic': True, 'scenario': args.scenario, 'molecules': args.molecules,
              'threads': args.threads, 'seconds': time.perf_counter()-started,
              'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024),
              'peak_child_rss_bytes': resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * (1 if sys.platform == 'darwin' else 1024),
              'reconstruction': json.loads((args.output/'result'/'reconstruction_metadata.json').read_text())['performance']}
    (args.output/'benchmark.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
