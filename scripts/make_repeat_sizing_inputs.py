#!/usr/bin/env python3
"""Create deterministic artificial DNA for assembly/FASTQ/MLVA_finder sizing checks.

Run: python scripts/make_repeat_sizing_inputs.py --output synthetic_sizing
No biological sequences, downloaded genomes, or known sample alleles are used.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import random
import sys

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mlvamaps.io import write_fastq
from mlvamaps.models import ReadRecord
from mlvamaps.sequence import revcomp


def make_inputs(root: Path):
    rng = random.Random(9302026)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    specifications = [(2, 0, 0, ''), (2, 1, 0, ''), (3, 9, 0, '')]
    specifications += [(4, 7, extra, '') for extra in range(4)]
    specifications += [(6, 13, 1, ''), (9, 3, 0, '')]
    specifications += [(12, 4, extra, '') for extra in (2, 3, 4, 8, 9, 10)]
    specifications += [(18, 9, 0, ''), (40, 3, 0, ''), (96, 3, 0, ''), (4, 80, 0, '')]
    specifications += [(12, 7, 0, change) for change in (
        'forward_insertion', 'forward_deletion', 'reverse_insertion', 'reverse_deletion',
        'flank_insertion', 'flank_deletion', 'reverse_orientation')]
    for sample in ('spanning', 'shotgun', 'unresolved'):
        output = root/sample
        assembly, database = output/'assemblies', output/'database'
        assembly.mkdir(parents=True, exist_ok=True)
        database.mkdir(exist_ok=True)
        primers, contigs, truth, single, first, second = [], [], [], [], [], []
        cases = specifications if sample == 'spanning' else [(u, n, extra, '') for u, n, extra in
            ((2, 0, 0), (2, 1, 0), (4, 7, 1), (4, 7, 3), (6, 9, 0), (9, 7, 0),
             (12, 4, 3), (12, 4, 9), (18, 4, 0), (40, 2, 0), (96, 1, 0))]
        if sample == 'unresolved':
            cases = [(4, 80, 0, '')]
        for i, (unit, count, extra, change) in enumerate(cases, 1):
            arm = 60 if sample == 'spanning' else 220
            left, right, motif = dna(arm), dna(arm), dna(unit)
            reference = left + motif*5 + right
            locus = f'S{i:02d}_{unit}bp_{len(reference)}bp_5U'
            primers.append(f'{locus}\t{left[:20]}\t{revcomp(right[-20:])}\n')
            (database/f'{locus}.fasta').write_text(f'>reference\n{reference}\n')
            sequence = left + motif*count + motif[:extra] + right
            if 'insertion' in change or 'deletion' in change:
                position = 8 if change.startswith('forward') else len(sequence)-8 if change.startswith('reverse') else 35
                sequence = (sequence[:position]+'A'+sequence[position:] if 'insertion' in change
                            else sequence[:position]+sequence[position+1:])
            calibrated_size = (2*arm + count*unit + extra +
                               (1 if change == 'flank_insertion' else -1 if change == 'flank_deletion' else 0))
            truth.append({'locus_id': locus, 'repeat_unit_bp': unit, 'repeat_tract_bp': count*unit+extra,
                          'product_size_bp': calibrated_size, 'sequence_size_bp': len(sequence),
                          'repeat_count_raw': 5+(calibrated_size-len(reference))/unit,
                          'change': change or 'none'})
            if change == 'reverse_orientation':
                sequence = revcomp(sequence)
            contigs.append(f'>{locus}\n{sequence}\n')
            if sample == 'spanning':
                for j in range(4):
                    single.append(ReadRecord(f'{locus}_{j}', sequence, 'I'*len(sequence)))
            elif sample == 'shotgun':
                # 150-bp mates from 250-bp fragments at staggered positions.
                # Neither reads nor fragments contain a complete PCR product.
                genome = dna(150)+sequence+dna(150)
                for j, start in enumerate(range(0, len(genome)-249, 10)):
                    fragment = genome[start:start+250]
                    first.append(ReadRecord(f'{locus}_{j}/1', fragment[:150], 'I'*150))
                    second.append(ReadRecord(f'{locus}_{j}/2', revcomp(fragment[-150:]), 'I'*150))
            else:
                # Opposite repeat boundaries, but no spanning sequence and no
                # independent fragment-length calibration: length is unknowable.
                for j in range(8):
                    first.append(ReadRecord(f'{locus}_{j}/1', sequence[arm-60:arm+90], 'I'*150))
                    second.append(ReadRecord(f'{locus}_{j}/2', revcomp(sequence[-arm-90:-arm+60]), 'I'*150))
        (assembly/f'{sample}.fasta').write_text(''.join(contigs))
        (output/'primers.tsv').write_text(''.join(primers))
        with (output/'truth.tsv').open('w') as handle:
            writer = csv.DictWriter(handle, delimiter='\t', fieldnames=list(truth[0]))
            writer.writeheader()
            writer.writerows(truth)
        if single:
            write_fastq(single, output/'single.fq')
        else:
            write_fastq(first, output/'r1.fq')
            write_fastq(second, output/'r2.fq')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    make_inputs(parser.parse_args().output)
