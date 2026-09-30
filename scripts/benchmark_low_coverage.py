#!/usr/bin/env python3
"""Synthetic sub-1x paired-read sweep against assembly-mode measurements.

This measures recovery and errors, not superiority over a whole-genome assembler.
The same reads are tested with a rich panel and a calibrated primer-only panel.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict, replace
import json
import math
from pathlib import Path
import random
import sys
import time

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mlvamaps.assembly_call import run_assembly_call
from mlvamaps.io import read_profiles, write_fasta, write_fastq, write_tsv
from mlvamaps.models import Locus, ReadRecord
from mlvamaps.sequence import revcomp
from mlvamaps.short_reads import run_short_read_call


def make_reference(output):
    rng = random.Random(418)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    loci, contigs = [], []
    for i, (unit, count) in enumerate(zip((4, 9, 12, 18, 24, 39), (5, 8.5, 10, 6, 12, 9))):
        left, right, motif = dna(40), dna(40), dna(unit)
        locus = Locus(f'L{i}', forward_primer=left[:20], reverse_primer=revcomp(right[-20:]),
            left_flank_sequence=left[20:], right_flank_sequence=right[:-20],
            repeat_motif=motif, repeat_unit_length_bp=unit, expected_min_repeats=2,
            expected_max_repeats=15, nominal_repeat_units=8, expected_product_size_bp=80+8*unit)
        product = left+(motif*20)[:round(count*unit)]+right
        loci.append(locus)
        contigs.append(dna(400)+product+dna(400))
    assembly = output/'truth.fasta'
    write_fasta(((f'contig{i}', sequence) for i, sequence in enumerate(contigs)), assembly)
    panels = {}
    for kind in ('rich', 'primers'):
        selected = loci if kind == 'rich' else [replace(l, left_flank_sequence='', right_flank_sequence='',
                                                       repeat_motif='N'*l.repeat_unit_length_bp) for l in loci]
        panels[kind] = output/f'{kind}.tsv'
        write_tsv([asdict(l) for l in selected], panels[kind], list(asdict(loci[0])))
    assemblies = {kind: output/f'assembly_{kind}' for kind in panels}
    truth = None
    for kind, panel in panels.items():
        paths = run_assembly_call(str(assembly), str(panel), str(assemblies[kind]), 'truth', threads=1)
        counts = {r['locus_id']: float(r['repeat_count']) for r in read_profiles(paths['calls'])}
        assert truth is None or counts == truth, 'Panels must preserve the same length calibration'
        truth = counts
    return contigs, panels, truth, assemblies


def make_reads(contigs, coverage, seed, error_rate):
    rng = random.Random(seed)
    total = sum(map(len, contigs))
    pairs = round(coverage*total/300)
    first, second = [], []
    quality = chr(33+max(2, min(40, round(-10*math.log10(error_rate))))) if error_rate else 'I'
    for index in range(pairs):
        fragment_size = max(150, min(450, round(rng.gauss(320, 30))))
        weights = [len(c)-fragment_size+1 for c in contigs]
        contig = rng.choices(contigs, weights=weights)[0]
        start = rng.randrange(len(contig)-fragment_size+1)
        fragment = contig[start:start+fragment_size]
        for mate, sequence, records in ((1, fragment[:150], first), (2, revcomp(fragment[-150:]), second)):
            sequence = ''.join(rng.choice([b for b in 'ACGT' if b != base])
                               if rng.random() < error_rate else base for base in sequence)
            records.append(ReadRecord(f'r{index}/{mate}', sequence, quality*len(sequence)))
    return first, second, pairs*300/total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--coverages', type=float, nargs='+', default=[.125, .25, .5, .75])
    parser.add_argument('--seeds', type=int, default=10)
    parser.add_argument('--error-rate', type=float, default=.001)
    args = parser.parse_args()
    if args.seeds < 1 or not 0 <= args.error_rate < 1 or any(not math.isfinite(c) or c <= 0 for c in args.coverages):
        parser.error('seeds and coverages must be positive; error rate must be in [0, 1)')
    args.output.mkdir(parents=True, exist_ok=True)
    contigs, panels, truth, assemblies = make_reference(args.output)
    results, summaries, manifest = [], defaultdict(lambda: defaultdict(float)), []
    for coverage in args.coverages:
        for seed in range(args.seeds):
            sample = f'coverage{coverage:g}_seed{seed}'
            directory = args.output/sample
            directory.mkdir(exist_ok=True)
            a, b, realized = make_reads(contigs, coverage, seed, args.error_rate)
            first, second = directory/'r1.fq', directory/'r2.fq'
            write_fastq(a, first)
            write_fastq(b, second)
            for kind, panel in panels.items():
                started = time.perf_counter()
                paths = run_short_read_call(str(first), str(second), str(panel), str(directory/kind), sample,
                                            threads=1, show_progress=False)
                elapsed = time.perf_counter()-started
                # Compare like-for-like SNP masking: assembly uses the same panel.
                comparison_sample = sample+'_'+kind
                manifest.append({'sample_id': comparison_sample, 'mode': 'assembly', 'outdir': str(assemblies[kind].resolve())})
                manifest.append({'sample_id': comparison_sample, 'mode': 'sr', 'outdir': str((directory/kind).resolve())})
                metrics = summaries[f'{coverage:g}:{kind}']
                metrics['runs'] += 1
                metrics['seconds'] += elapsed
                for row in read_profiles(paths['calls']):
                    observed = float(row['repeat_count']) if row['repeat_count'] else None
                    exact = observed == truth[row['locus_id']]
                    called = row['status'] == 'PASS'
                    metrics['loci'] += 1
                    metrics['called'] += called
                    metrics['correct_calls'] += called and exact
                    metrics['incorrect_calls'] += called and not exact
                    metrics['numeric_estimates'] += observed is not None
                    metrics['correct_estimates'] += exact
                    results.append({'sample_id': sample, 'panel': kind, 'coverage': coverage,
                        'realized_coverage': realized, 'locus_id': row['locus_id'], 'truth': truth[row['locus_id']],
                        'repeat_count': row['repeat_count'], 'status': row['status'], 'read_depth': row['read_depth']})
        print(f'Completed {coverage:g}x across {args.seeds} seeds and both panels', flush=True)
    for metrics in summaries.values():
        metrics['exact_call_recovery'] = metrics['correct_calls']/metrics['loci']
        metrics['error_rate_among_calls'] = metrics['incorrect_calls']/metrics['called'] if metrics['called'] else None
        metrics['missing_or_uncertain'] = metrics['loci']-metrics['called']
    write_tsv(results, args.output/'loci.tsv', list(results[0]))
    write_tsv(manifest, args.output/'runs.tsv', ['sample_id', 'mode', 'outdir'])
    report = {'synthetic': True, 'read_length': 150, 'fragment_mean': 320, 'fragment_sd': 30,
              'error_rate': args.error_rate, 'seeds': args.seeds, 'summaries': dict(summaries)}
    (args.output/'benchmark.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
