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


def test_long_motif_seed_filter_preserves_error_boundary_and_wrapping():
    from mlvamaps.repeat_likelihood import _seeded_repetitive
    rng = random.Random(813)
    for period in (128, 600, 1800):
        motif = ''.join(rng.choices('ACGT', k=period))
        for length in (9, 10, 19, 20, 149, 150, 151):
            phase = period - 4
            truth = ''.join(motif[(phase+i) % period] for i in range(length))
            for error_count in (length//10, length//10 + 1):
                sequence = list(truth)
                for index in rng.sample(range(length), error_count):
                    sequence[index] = next(base for base in 'ACGT' if base != sequence[index])
                sequence = ''.join(sequence)
                expected = max(sum(base == motif[(i+offset) % period] for i, base in enumerate(sequence))
                               for offset in range(period)) / length >= .9
                assert repetitive(sequence, motif) == expected
    motif = 'A'*1800 + 'C'
    sequence = 'C'*30 + 'A'*120
    assert _seeded_repetitive(sequence, motif) is None  # Bound dense seed-hit work.
    assert repetitive(sequence, motif) is False


def test_recovery_substage_progress_is_visible_in_stdout(tmp_path, capsys):
    from mlvamaps.io import write_fastq
    from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
    fastq = tmp_path/'reads.fq'
    write_fastq([p.read1 for p in reads(3)], fastq)
    _, _, _, paths = run_reconstructed_fastq_inference(reads1=fastq, reads2=None,
        loci=[t.locus for t in templates().values()], database_path=None,
        outdir=tmp_path/'out', sample_id='s', technology='illumina',
        minimum_molecules=3, minimum_probability=.8, threads=1, show_progress=True)
    output = capsys.readouterr().out
    assert output.index('Loading repeat templates') < output.index('Recruiting short-read molecules')
    assert 'Repeat fitting finished' in output
    assert 'Writing molecule memberships and allele predictions' in output
    assert 'Canonical genotype and evidence output finished' in output
    metadata = json.loads(paths['reconstruction_metadata'].read_text())
    assert metadata['performance']['stage_seconds']['template_loading'] >= 0
    assert set(metadata['performance']['output']) == {
        'genotyping_seconds', 'product_tables_seconds', 'canonical_alignments_seconds',
        'compatibility_tables_seconds'}


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
    original = parasail.sw_striped_profile_16
    def tracked(profile, reference, *args):
        scores.append(reference)
        return original(profile, reference, *args)
    monkeypatch.setattr(parasail, 'sw_striped_profile_16', tracked)
    result = infer_repeat([item], t, maximum=40)
    assert len(result.states) == 81
    assert 0 < len(scores) == len(set(scores)) < 81


def test_periodic_score_collapse_matches_full_native_alignments():
    from dataclasses import replace
    from mlvamaps.repeat_likelihood import _repeat_alignment_scores
    rng = random.Random(891)
    def dna(n):
        return ''.join(rng.choices('ACGT', k=n))
    base = templates()['L0']
    for unit, motif in ((1, 'A'), (3, 'ACG'), (4, 'AATG'), (7, dna(7)),
                        (60, dna(60)), (60, 'ACGT'*15), (5, 'ACGTACG')):
        template = replace(base, unit=unit, motif=motif)
        states = np.array([0, .5, 1, 2.5, 5, 8, 15, 30, 50, 100])
        target = template.sequence(8)
        cases = [target, target[:60], target[-60:], dna(75), 'N'*40,
                 target[:20] + 'NN' + target[23:80], target[:17] + 'A' + target[17:100]]
        # Deliberately mismatch declared unit and motif length, above, to test
        # right-boundary phase preservation as well as primitive motifs.
        for read in cases:
            expected = np.array([parasail.sw_striped_32(read, template.sequence(float(state)), 5, 1, _MATRIX).score
                                 for state in states], dtype=float)
            observed = _repeat_alignment_scores(read, template, states, np.empty(0))
            np.testing.assert_array_equal(observed, expected)
            prefix = _repeat_alignment_scores(read, template, states[:5], np.empty(0))
            np.testing.assert_array_equal(_repeat_alignment_scores(read, template, states, prefix), expected)


def test_repeat_fit_matches_exhaustive_scoring_for_mixtures_and_expansions(monkeypatch):
    from dataclasses import replace
    import mlvamaps.repeat_likelihood as module
    original = module._repeat_alignment_scores
    template = replace(templates()['L0'], unit=60, motif='AATG'*15)
    reads_at_locus = [template.sequence(2), template.sequence(3)] * 4
    reads_at_locus += [template.left + template.motif[:50], template.motif[:50] + template.right] * 4
    evidence = [classify_pair(ReadPair(str(i), ReadRecord(str(i), s, 'I'*len(s))), template)
                for i, s in enumerate(reads_at_locus)]
    evidence = [e for e in evidence if e and e.classes]
    def exhaustive(read, template, states, previous):
        return np.array([parasail.sw_striped_32(read, template.sequence(float(r)), 5, 1, _MATRIX).score
                         for r in states], dtype=float)
    for items in (evidence, [e for e in evidence if 'F' in e.classes]):
        monkeypatch.setattr(module, '_repeat_alignment_scores', exhaustive)
        before = infer_repeat(items, template, maximum=100)
        monkeypatch.setattr(module, '_repeat_alignment_scores', original)
        after = infer_repeat(items, template, maximum=100)
        for key, value in vars(before).items():
            if isinstance(value, np.ndarray):
                np.testing.assert_array_equal(getattr(after, key), value)
            else:
                assert getattr(after, key) == value


def test_projection_reuse_preserves_quality_filtering_and_molecule_support(monkeypatch):
    from dataclasses import replace
    import mlvamaps.locus_reconstruction as reconstruction
    from mlvamaps.repeat_likelihood import MoleculeEvidence, InsertDistribution
    template = templates()['L0']
    high = ('I'*len(template.left), 'I'*len(template.right))
    low = (high[0][:20] + '!' + high[0][21:], high[1])
    first = MoleculeEvidence('0', 'L0', ('S',), (template.left, template.right), ('+', '-'),
                             fragment_offset=len(template.left)+len(template.right), qualities=high)
    second = replace(first, molecule_id='1', qualities=low)
    cache = {}
    sequence = template.sequence(8)
    assert 20 in reconstruction._project_flanks(first, sequence, template, 8, cache)
    assert 20 not in reconstruction._project_flanks(second, sequence, template, 8, cache)
    items = [replace(first if i % 2 else second, molecule_id=str(i)) for i in range(50)]
    inference = infer_repeat(items, template, InsertDistribution(len(sequence), 1, 50, 'override'))
    calls = []
    original = reconstruction._project_flanks
    def tracked(*args):
        calls.append(True)
        return original(*args)
    monkeypatch.setattr(reconstruction, '_project_flanks', tracked)
    products = reconstruction.reconstruct_short_products(items, template, inference, 's')
    assert len(calls) == 2
    assert len(products) == 1
    assert products[0].sequence == sequence
    assert products[0].support_count == 50
    assert products[0].evidence['molecule_ids'] == [str(i) for i in range(50)]


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


@pytest.mark.parametrize('audit_mode', ['full', 'compact'])
def test_parallel_recruitment_is_ordered_and_matches_serial(audit_mode):
    pairs = list(reads())
    # Deliberately duplicate locus context to exercise ambiguous-molecule output.
    ts = templates()
    from dataclasses import replace
    ts['L1'] = replace(ts['L0'], locus=replace(ts['L0'].locus, locus_id='L1'))
    def collect(workers):
        stats = {}
        result = list(recruit_short_reads(iter(pairs), ts, 's', workers, chunk_size=4,
                                         statistics=stats, audit_mode=audit_mode))
        return [asdict(chunk) for chunk in result], stats
    serial, _ = collect(1)
    parallel, stats = collect(2)
    assert serial == parallel
    assert stats['workers'] == 2
    assert sum(c['examined'] for c in parallel) == len(pairs)


def test_anchor_existence_matches_complete_classification():
    from mlvamaps.repeat_likelihood import pair_has_anchor
    rng = random.Random(810)
    cases = list(reads())
    # Orphans, repeat-only reads, reverse reads, ambiguous bases and indels.
    t = templates()['L0']
    sequences = [t.motif*40, revcomp(t.motif*40), 'N'*150, '', t.left,
                 t.sequence(5)[:30] + 'N' + t.sequence(5)[32:]]
    sequences += [''.join(rng.choices('ACGTN', k=150)) for _ in range(100)]
    for i, sequence in enumerate(sequences):
        cases.append(ReadPair(str(i), ReadRecord(str(i), sequence)))
        cases.append(ReadPair(str(i), ReadRecord(str(i), 'AATG'*30), ReadRecord(str(i), sequence)))
    for pair in cases:
        for template in templates().values():
            assert pair_has_anchor(pair, template) == (classify_pair(pair, template) is not None)


def test_compact_recruitment_preserves_unique_evidence_and_insert_calibration():
    from dataclasses import replace
    ts = templates()
    ts['L1'] = replace(ts['L0'], locus=replace(ts['L0'].locus, locus_id='L1'))
    pairs = list(reads())
    flank = ts['L2'].left
    pairs.append(ReadPair('insert', ReadRecord('insert/1', flank[:25]),
                          ReadRecord('insert/2', revcomp(flank[-25:]))))
    full = ShortReadRecruiter(ts, 's', audit_mode='full').recruit(pairs)
    compact = ShortReadRecruiter(ts, 's', audit_mode='compact').recruit(pairs)
    assert full.evidence and full.ambiguous_pairs and full.unmatched_pairs and full.insert_lengths
    assert compact.evidence == full.evidence
    assert compact.insert_lengths == full.insert_lengths
    assert compact.ambiguous_pairs == full.ambiguous_pairs
    assert compact.unmatched_pairs == full.unmatched_pairs
    assert compact.examined == len(compact.evidence) + compact.ambiguous_pairs + compact.unmatched_pairs
    assert compact.locus_tests + compact.skipped_locus_tests == full.locus_tests
    assert compact.skipped_locus_tests == 0
    assert len(compact.ambiguity_witnesses) == compact.ambiguous_pairs
    full_hits = {}
    for row in full.ambiguous:
        full_hits.setdefault(row['molecule_id'], []).append(row['locus_id'])
    for molecule, first, second, complete in compact.ambiguity_witnesses:
        assert [first, second] == full_hits[molecule][:2]
        assert complete == 'yes'
    # Two candidates really do complete the search; absent templates do not count.
    compact = ShortReadRecruiter({'L0': ts['L0'], 'L1': ts['L1'], 'absent': None}, 's',
                                audit_mode='compact').recruit(pairs[:1])
    assert compact.ambiguity_witnesses == [('0', 'L0', 'L1', 'yes')]
    assert compact.skipped_locus_tests == 0


@pytest.mark.parametrize('stream', [False, True])
def test_compact_recovery_preserves_genotypes_and_writes_ambiguity_witnesses(tmp_path, stream):
    import csv
    import gzip
    from dataclasses import replace
    from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
    ts = templates()
    ts['L1'] = replace(ts['L0'], locus=replace(ts['L0'].locus, locus_id='L1'))
    pairs = list(reads(30))
    results = {}
    for mode in ('full', 'compact'):
        results[mode] = run_reconstructed_fastq_inference(reads1=None, reads2=None, pairs=iter(pairs),
            loci=[t.locus for t in ts.values()], database_path=None, outdir=tmp_path/mode,
            sample_id='s', technology='illumina', minimum_molecules=3, minimum_probability=.8,
            sr_recruitment_audit=mode, stream_molecule_evidence=stream)
    full, compact = results['full'], results['compact']
    assert full[0] == compact[0]
    for key, path in full[3].items():
        if key in {'reconstruction_metadata', 'molecule_evidence'}:
            continue
        def contents(p):
            return gzip.decompress(p.read_bytes()) if p.suffix == '.gz' else p.read_bytes()
        assert contents(path) == contents(compact[3][key]), key
    def rows(path):
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter='\t'))
    full_audit = rows(full[3]['molecule_evidence'])
    assert rows(compact[3]['molecule_evidence']) == [r for r in full_audit if r['classes'] != 'ambiguous_locus']
    witnesses = rows(compact[3]['ambiguous_molecules'])
    assert witnesses and all(r['sample_id'] == 's' for r in witnesses)
    assert {r['molecule_id'] for r in witnesses} == {r['molecule_id'] for r in full_audit if r['classes'] == 'ambiguous_locus'}
    stats = json.loads(compact[3]['reconstruction_metadata'].read_text())['performance']['recruitment']
    assert stats['ambiguous_pairs'] == len(witnesses)
    assert stats['audit_mode'] == 'compact'
    assert stats['pairs_per_second'] > 0
    assert stats['skipped_locus_tests'] == 0


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


def test_byte_lane_anchors_match_wide_alignments_with_indels_and_saturation(monkeypatch):
    from mlvamaps.repeat_likelihood import FlankHit
    rng = random.Random(481)
    def dna(n):
        return ''.join(rng.choices('ACGTN', k=n))
    fallbacks = []
    original = parasail.sw_trace_striped_profile_32
    def tracked(*args):
        fallbacks.append(True)
        return original(*args)
    monkeypatch.setattr(parasail, 'sw_trace_striped_profile_32', tracked)
    for i in range(300):
        flank = dna(rng.randint(4, 350))
        seq = list(flank)
        for _ in range(i % 17):
            position = rng.randrange(len(seq))
            if i % 2:
                seq.insert(position, rng.choice('ACGTN'))
            else:
                seq[position] = rng.choice('ACGTN')
        seq = dna(150) if i % 3 == 0 else ''.join(seq)
        minimum = (4, 10, 14, 20)[i % 4]
        result = parasail.sw_trace_striped_32(seq, flank, 5, 1, _MATRIX)
        trace = result.traceback
        q, r = trace.query, trace.ref
        aligned = sum(a != '-' and b != '-' for a, b in zip(q, r))
        edits = sum(a != b for a, b in zip(q, r))
        expected = None
        if aligned >= min(minimum, len(flank)) and edits / max(aligned, 1) <= .15:
            expected = FlankHit(result.end_query + 1 - len(q.replace('-', '')), result.end_query + 1,
                                result.end_ref + 1 - len(r.replace('-', '')), result.end_ref + 1,
                                result.score, edits)
        assert flank_hit(seq, flank, minimum) == expected
    assert fallbacks


def test_worker_formatted_audit_matches_original_tsv():
    import csv
    import io
    from mlvamaps.locus_reconstruction import MOLECULE_AUDIT_FIELDS
    ts = templates()
    pairs = list(reads(20))
    # Identical contexts guarantee ambiguous rows for every targeted molecule.
    from dataclasses import replace
    ts['duplicate'] = replace(ts['L0'], locus=replace(ts['L0'].locus, locus_id='duplicate'))
    expected = io.StringIO(newline='')
    writer = csv.DictWriter(expected, fieldnames=MOLECULE_AUDIT_FIELDS, delimiter='\t')
    for chunk in recruit_short_reads(pairs, ts, 's', 1, chunk_size=4):
        writer.writerows(chunk.ambiguous)
    assert expected.getvalue()
    actual = list(recruit_short_reads(pairs, ts, 's', 2, chunk_size=4, audit_fields=MOLECULE_AUDIT_FIELDS))
    assert ''.join(chunk.audit_tsv for chunk in actual) == expected.getvalue()
    assert not any(chunk.ambiguous for chunk in actual)


@pytest.mark.parametrize('materialize', [False, True])
def test_qc_stream_matches_disk_round_trip_with_orphans(tmp_path, materialize):
    from mlvamaps.io import read_fastq, read_fastq_pairs, write_fastq
    from mlvamaps.short_read_qc import filtered_pairs
    from mlvamaps.short_reads import qc_read_pairs
    from dataclasses import replace
    pairs = list(reads(12))
    pairs[1] = replace(pairs[1], read1=replace(pairs[1].read1, quality='!'*len(pairs[1].read1.sequence)))
    pairs[2] = replace(pairs[2], read2=replace(pairs[2].read2, quality='!'*len(pairs[2].read2.sequence)))
    pairs[3] = replace(pairs[3], read1=replace(pairs[3].read1, sequence='A', quality='I'),
                       read2=replace(pairs[3].read2, sequence='T', quality='I'))
    first, second = tmp_path/'r1.fq', tmp_path/'r2.fq'
    write_fastq([p.read1 for p in pairs], first)
    write_fastq([p.read2 for p in pairs], second)
    f1, f2, orphan = (tmp_path/name for name in ('f1.fq.gz', 'f2.fq.gz', 'o.fq.gz'))
    expected, metrics = qc_read_pairs(pairs, 40, 15, 0, .5)
    counters, statistics = {}, {}
    stream = filtered_pairs(first, second, f1, f2, orphan, materialize=materialize,
        counters=counters, statistics=statistics, min_length=40, min_mean_quality=15,
        trim_quality=0, min_pair_retention=.5, sample_id='s', chunk_size=2)
    first_pair = next(stream)
    assert counters['input_pairs'] == 2  # Inference can start before QC ends.
    actual = [first_pair, *stream]
    expected_order = [p for p in expected if p.read2] + [
        ReadPair(p.read1.read_id, p.read1) for p in expected if p.read2 is None]
    assert actual == expected_order
    assert counters == metrics
    assert statistics['seconds'] >= 0
    assert f1.exists() == f2.exists() == materialize
    if materialize:
        disk = list(read_fastq_pairs(f1, f2)) + [ReadPair(r.read_id, r) for r in read_fastq(orphan)]
        assert actual == disk


def test_qc_mean_quality_matches_original_numpy_reduction():
    from mlvamaps.sequence import mean_qscore
    rng = random.Random(102)
    assert mean_qscore(None) == mean_qscore('') == 0
    for length in (1, 40, 150, 151, 250, 10000):
        quality = ''.join(chr(rng.randint(33, 74)) for _ in range(length))
        expected = float(np.frombuffer(quality.encode('ascii'), dtype=np.uint8).mean(dtype=np.float64) - 33)
        assert mean_qscore(quality) == expected


def test_streaming_pipeline_preserves_database_inputs_and_results(tmp_path, monkeypatch):
    from mlvamaps.io import read_fastq_pairs, write_fastq, write_tsv
    from mlvamaps.short_reads import run_short_read_call
    pairs = list(reads(12))
    first, second = tmp_path/'r1.fq', tmp_path/'r2.fq'
    write_fastq([p.read1 for p in pairs], first)
    write_fastq([p.read2 for p in pairs], second)
    panel = tmp_path/'panel.tsv'
    loci = [asdict(t.locus) for t in templates().values()]
    write_tsv(loci, panel, list(loci[0]))
    classified = []
    def classify(**kwargs):
        classified.extend(read_fastq_pairs(kwargs['reads1'], kwargs['reads2']))
        return {}
    monkeypatch.setattr('mlvamaps.mapping_classification.run_mapping_classification', classify)
    monkeypatch.setattr('mlvamaps.minimap_mapping.minimap2_version', lambda _: 'test')
    for database in (None, 'database'):
        run_short_read_call(str(first), str(second), str(panel), str(tmp_path/str(database)), 's',
                            database_path=database, threads=1, show_progress=False)
    assert classified == pairs
    for name in ('calls.tsv', 'reconstructed_loci.fasta', 'reconstructed_locus_variants.tsv',
                 'locus_snps.tsv', 'short_read_qc_summary.tsv', 'molecule_candidate_evidence.tsv'):
        assert (tmp_path/'None'/name).read_bytes() == (tmp_path/'database'/name).read_bytes()
    metadata = json.loads((tmp_path/'database'/'short_read_run_metadata.json').read_text())
    assert metadata['performance']['stage_seconds']['total'] > 0
    assert metadata['performance']['qc']['materialized'] is True
