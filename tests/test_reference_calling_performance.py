"""Reference recruitment scores and calls must survive streaming/parallel fitting."""
from collections import defaultdict
import json
import os
import random
import shutil
import sys

import pytest

from mlvamaps.candidate_contexts import CandidateContext, generate_candidate_contexts, write_candidate_contexts
from mlvamaps.io import write_fastq
from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
from mlvamaps.minimap_mapping import build_minimap2_index, parse_candidate_alignments, run_minimap2_competitive_bam
from mlvamaps.models import Locus, ReadRecord
from mlvamaps.sequence import revcomp


@pytest.mark.parametrize('keep', [False, True])
def test_streamed_scores_match_full_alignment_reduction(tmp_path, keep):
    contexts = [CandidateContext(name, 'L', 'A'*100, 5, 4, 40, 60, 'A'*40, 'A'*40)
                for name in ('a', 'b', 'c')]
    groups = {'a': ('L', 'same'), 'b': ('L', 'same'), 'c': ('L', 'other')}
    records = [('pair/1', 65, 'a', 10), ('pair/1', 321, 'b', 20),
               ('pair/2', 129, 'a', 15), ('pair/2', 385, 'b', 12),
               ('pair/1', 321, 'c', 7), ('pair/1', 2048, 'c', 100),
               ('named/2', 0, 'c', -2), ('no_AS', 0, 'a', None)]
    sam = tmp_path/'input.sam'
    sam.write_text('@HD\tVN:1.6\n'+''.join(f'@SQ\tSN:{c.candidate_id}\tLN:100\n' for c in contexts)
        + ''.join(f'{name}\t{flag}\t{ref}\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII'
                  + (f'\tAS:i:{score}' if score is not None else '') + '\n'
                  for name, flag, ref, score in records)
        + 'unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\n')
    expected = {'pair': {('L', 'same'): {1: 20., 2: 15.}, ('L', 'other'): {1: 7.}},
                'named': {('L', 'other'): {2: -2.}}, 'no_AS': {('L', 'same'): {None: 10.}}}
    legacy = defaultdict(lambda: defaultdict(dict))
    for row in parse_candidate_alignments(sam, contexts):
        mates = legacy[row.molecule_id][groups[row.candidate_id]]
        mates[row.mate] = max(mates.get(row.mate, -float('inf')), row.alignment_score)
    assert legacy == expected
    bam, statistics, updates = tmp_path/'kept.bam', {}, []
    command = [sys.executable, '-c', 'import sys; sys.stdout.write(open(sys.argv[1]).read())', str(sam)]
    actual = run_minimap2_competitive_bam(command, bam if keep else None,
        score_groups=groups, statistics=statistics, on_progress=updates.append)
    assert actual == expected
    assert statistics == {'alignment_records': 9, 'mapped_molecules': 3}
    assert updates[-1] == 9
    assert bam.exists() == keep
    if keep:
        assert parse_candidate_alignments(bam, contexts) == parse_candidate_alignments(sam, contexts)


def test_reference_process_calls_match_serial_and_subset_rebuilds_index(tmp_path):
    minimap = shutil.which('minimap2')
    if minimap is None:
        pytest.skip('minimap2 is required for reference calling integration')
    rng = random.Random(88720)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    database = tmp_path/'database'
    database.mkdir()
    loci, first, second = [], [], []
    for i in range(4):
        left, right, motif = dna(180), dna(180), dna(9)
        loci.append(Locus(f'L{i}', forward_primer=left[:20], reverse_primer=revcomp(right[-20:]),
            repeat_motif=motif, repeat_unit_length_bp=9, expected_product_size_bp=432,
            nominal_repeat_units=8, expected_min_repeats=2, expected_max_repeats=35))
        reference = left+motif*8+right
        # Equally supported backgrounds differ outside the observed reads.
        alternate = reference[:38]+next(b for b in 'ACGT' if b != reference[38])+reference[39:]
        (database/f'L{i}.fasta').write_text(f'>reference\n{reference}\n>alternate\n{alternate}\n')
        sequence = left+motif*11+right
        for j in range(40):
            start = 150 if i % 2 == 0 else 110+j%50
            a = sequence[start:start+(180 if i % 2 == 0 else 150)]
            b = revcomp(sequence[start+70:start+220])
            if i % 2:
                a, b = left[-60:]+motif*7, revcomp(right[:80])
            first.append(ReadRecord(f'L{i}_{j}/1', a, 'I'*len(a)))
            second.append(ReadRecord(f'L{i}_{j}/2', b, 'I'*len(b)))
    contexts = generate_candidate_contexts(loci, database)
    assets = write_candidate_contexts(contexts, database/'competitive_mapping')
    build_minimap2_index(assets['fasta'], database/'competitive_mapping'/'short.mmi',
                        executable=minimap, kmer_size=21, window_size=11)
    write_fastq(first, tmp_path/'r1.fq')
    write_fastq(second, tmp_path/'r2.fq')
    def run(name, threads, panel):
        calls, _, _, paths = run_reconstructed_fastq_inference(reads1=tmp_path/'r1.fq', reads2=tmp_path/'r2.fq',
            loci=panel, database_path=database, outdir=tmp_path/name, sample_id='s', technology='illumina',
            minimum_molecules=0, minimum_probability=.8, threads=threads, minimap2_bin=minimap,
            show_progress=False)
        return calls, paths, json.loads(paths['reference_calling_performance'].read_text())
    serial, serial_paths, serial_stats = run('serial', 1, loci)
    parallel, parallel_paths, parallel_stats = run('parallel', 3, loci)
    assert serial == parallel
    assert [row['repeat_count'] for row in serial] == [11, '', 11, '']
    assert serial_stats['cached_index'] and parallel_stats['cached_index']
    assert 1 <= parallel_stats['workers_used'] <= 3
    assert all(row['worker_pid'] != os.getpid() for row in parallel_stats['loci'].values())
    assert serial_stats['workers_used'] == 1
    for key in ('reference_calling', 'reconstructed_locus_variants', 'mapped_read_memberships',
                'read_predictions'):
        assert serial_paths[key].read_bytes() == parallel_paths[key].read_bytes()
    assert not list((tmp_path/'parallel'/'reference_calling').glob('*.bam'))
    subset, _, subset_stats = run('subset', 1, loci[:1])
    assert not subset_stats['cached_index']  # Other panel loci must not enter this competition.
    assert subset == serial[:1]
