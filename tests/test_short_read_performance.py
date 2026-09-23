"""Equivalence and resource-bound checks for the optimized short-read path."""
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
import json
import random

import numpy as np
import parasail
import pytest

from mlvamaps.concurrency import bounded_ordered_map
from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.repeat_likelihood import (
    _MATRIX, _scores, classify_pair, flank_hit, infer_repeat, panel_template, repetitive,
)
from mlvamaps.short_read_recruitment import ShortReadRecruiter, recruit_short_reads
from mlvamaps.sequence import revcomp


def templates():
    rng = random.Random(19)
    def dna(n):
        return ''.join(rng.choices('ACGT', k=n))
    return {f'L{i}': panel_template(Locus(f'L{i}', forward_primer=dna(16), reverse_primer=dna(16),
            left_flank_sequence=dna(30), right_flank_sequence=dna(30), repeat_motif='AATG',
            repeat_unit_length_bp=4, expected_max_repeats=8)) for i in range(3)}


def reads(count=20):
    ts = list(templates().values())
    for i in range(count):
        t = ts[i % len(ts)]
        sequence = t.sequence(4 + i % 3)
        if i % 2:
            sequence = revcomp(sequence)
        yield ReadPair(str(i), ReadRecord(f'{i}/1', sequence, 'I' * len(sequence)),
                       ReadRecord(f'{i}/2', revcomp(sequence), 'I' * len(sequence)))


def test_vectorized_motif_test_matches_original():
    rng = random.Random(29)
    cases = [('', 'A'), ('A', ''), ('A'*9+'C', 'A'), ('A'*8+'CC', 'A'), ('AATG'*20, 'AATG'*8)]
    for length in (1, 10, 30, 150):
        for unit in (1, 4, 12, 60):
            motif = ''.join(rng.choices('ACGT', k=unit))
            for offset in (0, unit // 2):
                truth = ''.join(motif[(p+offset) % unit] for p in range(length))
                cases.append((truth, motif))
                changed = list(truth)
                for p in rng.sample(range(length), min(length, max(1, length//10))):
                    changed[p] = next(b for b in 'ACGT' if b != changed[p])
                cases.append((''.join(changed), motif))
                cases.append((''.join(rng.choices('ACGTN', k=length)), motif))
    # Exercise bounded phase blocks rather than the cached small matrix path.
    motif = ''.join(rng.choices('ACGT', k=1700))
    cases.extend([(motif[120:270], motif), ('N'*150, motif)])
    for sequence, motif in cases:
        expected = bool(sequence and motif) and max(
            sum(base == motif[(i+offset) % len(motif)] for i,base in enumerate(sequence))
            for offset in range(len(motif))) / len(sequence) >= .9
        assert repetitive(sequence, motif) == expected


def test_profile_alignments_preserve_original_traceback_and_likelihoods():
    t = templates()['L0']
    rng = random.Random(37)
    for sequence in [t.sequence(5), revcomp(t.sequence(6)), t.left[8:]+t.motif*4] + [
            ''.join(rng.choices('ACGT', k=90)) for _ in range(20)]:
        result = parasail.sw_trace_striped_32(sequence, t.left, 5, 1, _MATRIX)
        trace = result.traceback
        q, r = trace.query, trace.ref
        aligned = sum(a != '-' and b != '-' for a,b in zip(q,r))
        edits = sum(a != b for a,b in zip(q,r))
        observed = flank_hit(sequence, t.left)
        if aligned < 10 or edits/max(aligned,1) > .15:
            assert observed is None
        else:
            assert asdict(observed) == dict(query_start=result.end_query+1-len(q.replace('-','')),
                query_end=result.end_query+1, flank_start=result.end_ref+1-len(r.replace('-','')),
                flank_end=result.end_ref+1, score=result.score, edits=edits)
    seq = t.sequence(6)
    item = classify_pair(ReadPair('m', ReadRecord('m', seq)), t)
    states = np.arange(0, 10.5, .5)
    scalar = np.array([parasail.sw_striped_32(seq, t.sequence(r), 5, 1, _MATRIX).score for r in states])
    expected = -.5*((states-item.observed_repeat)/.20)**2 + (scalar-scalar.max())/8
    np.testing.assert_array_equal(_scores([item], t, states, None)[0], expected)


def test_range_expansion_reuses_existing_candidate_alignments(monkeypatch):
    t = templates()['L0']
    seq = t.left + t.motif*5
    item = classify_pair(ReadPair('m', ReadRecord('m', seq)), t)
    scores = []
    original = parasail.sw_striped_profile_32
    def tracked(profile, reference, *args):
        scores.append(reference)
        return original(profile, reference, *args)
    monkeypatch.setattr(parasail, 'sw_striped_profile_32', tracked)
    result = infer_repeat([item], t, maximum=40)
    assert len(result.states) == 81
    assert len(scores) == len(set(scores)) == 81


def test_bit_mask_seed_filter_preserves_exhaustive_matches():
    ts = templates()
    pairs = list(reads())
    chunk = ShortReadRecruiter(ts, 's').recruit(pairs)
    expected_evidence, expected_ambiguous = [], []
    for pair in pairs:
        matches = [item for _,t in sorted(ts.items()) if (item := classify_pair(pair,t))]
        if len(matches) == 1:
            expected_evidence.append(asdict(matches[0]))
        elif matches:
            expected_ambiguous.extend((item.molecule_id,item.locus_id) for item in matches)
    assert [asdict(e) for e in chunk.evidence] == expected_evidence
    assert [(row['molecule_id'],row['locus_id']) for row in chunk.ambiguous] == expected_ambiguous


def test_parallel_recruitment_is_ordered_and_matches_serial():
    pairs = list(reads())
    # Deliberately duplicate locus context to exercise ambiguous-molecule output.
    ts = templates()
    def collect(workers):
        stats = {}
        result = list(recruit_short_reads(iter(pairs), ts, 's', workers, chunk_size=4, statistics=stats))
        return [asdict(chunk) for chunk in result], stats
    serial, _ = collect(1)
    parallel, stats = collect(2)
    assert serial == parallel
    assert stats['workers'] == 2
    assert sum(c['examined'] for c in parallel) == len(pairs)


def test_bounded_map_does_not_eagerly_consume_input():
    consumed = []
    def source():
        for i in range(30):
            consumed.append(i)
            yield i
    with ThreadPoolExecutor(max_workers=2) as executor:
        mapped = bounded_ordered_map(executor, abs, source(), 4)
        assert next(mapped) == 0
        assert len(consumed) == 4
        assert list(mapped) == list(range(1,30))
    with ThreadPoolExecutor(max_workers=1) as executor:
        with pytest.raises(ValueError, match='max_pending'):
            list(bounded_ordered_map(executor, abs, [], 0))


def test_streamed_audit_matches_in_memory_files(tmp_path):
    from mlvamaps.io import write_fastq
    from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
    pairs = list(reads())
    first, second = tmp_path/'r1.fq', tmp_path/'r2.fq'
    write_fastq([p.read1 for p in pairs], first)
    write_fastq([p.read2 for p in pairs], second)
    results = []
    for stream in (False, True):
        result = run_reconstructed_fastq_inference(reads1=first, reads2=second,
            loci=[t.locus for t in templates().values()], database_path=None,
            outdir=tmp_path/str(stream), sample_id='s', technology='illumina',
            minimum_molecules=3, minimum_probability=.8, threads=1, stream_molecule_evidence=stream)
        results.append(result)
    assert results[0][0] == results[1][0]
    assert results[0][1] and not results[1][1]
    for key,path in results[0][3].items():
        if key != 'reconstruction_metadata':
            # Gzip header timestamps/filenames can differ; compare logical bytes.
            import gzip
            if path.suffix == '.gz':
                assert gzip.decompress(path.read_bytes()) == gzip.decompress(results[1][3][key].read_bytes())
            else:
                assert path.read_bytes() == results[1][3][key].read_bytes()
    stats = json.loads(results[1][3]['reconstruction_metadata'].read_text())['performance']
    assert stats['recruitment']['pairs_examined'] == len(pairs)
    assert stats['stage_seconds']['recruitment'] >= 0


def test_recruitment_processes_can_start_inside_batch_sample_threads():
    pairs = list(reads(12))
    ts = templates()
    def sample():
        return [asdict(chunk) for chunk in recruit_short_reads(pairs, ts, 's', threads=2, chunk_size=4)]
    with ThreadPoolExecutor(max_workers=2) as executor:
        first = executor.submit(sample)
        second = executor.submit(sample)
        assert first.result() == second.result()
