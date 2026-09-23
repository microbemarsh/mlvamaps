"""Regression coverage for incidental flank matches suppressing real loci."""
from dataclasses import replace
import random

import pytest

from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.repeat_likelihood import classify_pair, panel_template
from mlvamaps.sequence import revcomp
from mlvamaps.short_read_recruitment import ShortReadRecruiter, anchor_score, resolve_locus_matches


def panel():
    rng = random.Random(19)

    def dna(length):
        return ''.join(rng.choices('ACGT', k=length))

    return {f'L{i:02}': panel_template(Locus(
        f'L{i:02}', forward_primer=dna(20), reverse_primer=dna(20),
        left_flank_sequence=dna(80), right_flank_sequence=dna(80),
        repeat_motif='AATG', repeat_unit_length_bp=4, expected_max_repeats=8,
    )) for i in range(25)}


def paired_reads(template, name):
    sequence = template.sequence(6)
    return ReadPair(name, ReadRecord(name+'/1', sequence[:150], 'I'*150),
                    ReadRecord(name+'/2', revcomp(sequence[-150:]), 'I'*150))


@pytest.mark.parametrize('audit_mode', ['full', 'compact'])
def test_incidental_matches_do_not_hide_perfect_multilocus_calls(tmp_path, audit_mode):
    templates = panel()
    pairs = [paired_reads(t, f'{name}_{i}') for name, t in templates.items() for i in range(6)]
    # The original any-two-hits rule discards every one of these molecules.
    assert all(sum(classify_pair(pair, t) is not None for t in templates.values()) > 1
               for pair in pairs[::6])
    chunk = ShortReadRecruiter(templates, 's', audit_mode=audit_mode).recruit(pairs)
    assert chunk.ambiguous_pairs == chunk.unmatched_pairs == 0
    assert [(e.molecule_id, e.locus_id) for e in chunk.evidence] == [
        (pair.molecule_id, pair.molecule_id.split('_')[0]) for pair in pairs]
    calls, _, _, _ = run_reconstructed_fastq_inference(
        reads1=None, reads2=None, pairs=iter(pairs), loci=[t.locus for t in templates.values()],
        database_path=None, outdir=tmp_path, sample_id='s', technology='illumina',
        minimum_molecules=3, minimum_probability=.8, sr_recruitment_audit=audit_mode,
    )
    assert len(calls) == 25
    assert all(row['status'] == 'called' and row['repeat_count'] == 6 for row in calls)


@pytest.mark.parametrize('audit_mode', ['full', 'compact'])
@pytest.mark.parametrize('mismatches', [0, 1, 3])
def test_identical_and_near_tied_loci_remain_ambiguous(audit_mode, mismatches):
    original = panel()['L00']
    left = list(original.left)
    for position in (30, 40, 50)[:mismatches]:
        left[position] = next(base for base in 'ACGT' if base != left[position])
    related = replace(original, locus=replace(original.locus, locus_id='related'), left=''.join(left))
    chunk = ShortReadRecruiter({'original': original, 'related': related}, 's',
                              audit_mode=audit_mode).recruit([paired_reads(original, 'm')])
    assert not chunk.evidence
    assert chunk.ambiguous_pairs == 1


@pytest.mark.parametrize('audit_mode', ['full', 'compact'])
def test_uncontested_short_anchor_is_retained(audit_mode):
    template = panel()['L00']
    pair = ReadPair('short', ReadRecord('short', template.left[-10:] + 'AATG'*10))
    chunk = ShortReadRecruiter({'L00': template}, 's', audit_mode=audit_mode).recruit([pair])
    assert len(chunk.evidence) == 1
    assert 'F' in chunk.evidence[0].classes


@pytest.mark.parametrize('audit_mode', ['full', 'compact'])
def test_score_bounds_preserve_exhaustive_evidence_and_audits(audit_mode):
    from mlvamaps.repeat_likelihood import pair_anchor_score
    rng = random.Random(314)
    templates = panel()
    templates['L01'] = replace(templates['L00'], locus=templates['L01'].locus)
    pairs = [paired_reads(t, name) for name, t in templates.items()]
    # Include background, orphans, indels, low qualities and repeat-only mates.
    for i, pair in enumerate(list(pairs)):
        seq = pair.read1.sequence
        seq = seq[:40] + ('NN' if i % 2 else 'A') + seq[41:]
        pairs.append(replace(pair, molecule_id=f'error{i}',
                             read1=replace(pair.read1, sequence=seq, quality='!'*len(seq))))
        pairs.append(ReadPair(f'orphan{i}', ReadRecord(str(i), seq.lower())))
        pairs.append(ReadPair(f'background{i}', ReadRecord(str(i), ''.join(rng.choices('ACGT', k=150)))))
    t = templates['L00']
    pairs.extend([
        ReadPair('frr', ReadRecord('1', t.left), ReadRecord('2', revcomp(t.motif*40))),
        ReadPair('repeat_only', ReadRecord('1', t.motif*40)),
        ReadPair('unknown', ReadRecord('1', 'N'*150)),
    ])
    expected, ambiguous, witnesses = [], [], []
    unmatched = ambiguous_pairs = 0
    for pair in pairs:
        matches = []
        for _, template in sorted(templates.items()):
            item = classify_pair(pair, template)
            assert pair_anchor_score(pair, template) == (anchor_score(item) if item else 0)
            if item:
                matches.append(item)
        selected = resolve_locus_matches(matches)
        if selected:
            expected.append(selected)
        elif matches:
            ambiguous_pairs += 1
            ambiguous.extend((item.molecule_id, item.locus_id, item.alignment) for item in matches)
            ranked = sorted(matches, key=anchor_score, reverse=True)
            witnesses.append((pair.molecule_id, ranked[0].locus_id, ranked[1].locus_id, 'yes'))
        else:
            unmatched += 1
    actual = ShortReadRecruiter(templates, 's', audit_mode=audit_mode).recruit(pairs)
    assert actual.evidence == expected
    assert actual.ambiguous_pairs == ambiguous_pairs > 0
    assert actual.unmatched_pairs == unmatched > 0
    assert actual.skipped_locus_tests > 0
    # The all-N read has no flank seeds and never enters candidate scoring.
    assert actual.locus_tests + actual.skipped_locus_tests == (len(pairs)-1)*len(templates)
    if audit_mode == 'compact':
        assert actual.ambiguity_witnesses == witnesses
    else:
        import json
        assert [(r['molecule_id'], r['locus_id'], json.loads(r['alignment']))
                for r in actual.ambiguous] == ambiguous


def test_score_bound_matches_unclipped_alignment():
    import parasail
    from mlvamaps.repeat_likelihood import _MATRIX, flank_score_bound
    rng = random.Random(614)
    for length in (0, 9, 30, 150, 500):
        flank = ''.join(rng.choices('ACGTN', k=length))
        for seq in (flank, revcomp(flank), flank[:20] + 'AAA' + flank[22:]):
            expected = parasail.sw_striped_32(seq, flank, 5, 1, _MATRIX).score if seq and flank else 0
            assert flank_score_bound(seq, flank) == expected


def test_compact_path_only_classifies_retained_molecules(monkeypatch):
    import mlvamaps.short_read_recruitment as module
    templates = panel()
    # Two identical top contexts must remain ambiguous without building evidence.
    templates['L01'] = replace(templates['L00'], locus=templates['L01'].locus)
    pairs = [paired_reads(t, name) for name, t in templates.items()]
    classified = []
    original = module.classify_pair

    def tracked(pair, template):
        classified.append((pair.molecule_id, template.locus.locus_id))
        return original(pair, template)

    monkeypatch.setattr(module, 'classify_pair', tracked)
    chunk = ShortReadRecruiter(templates, 's', audit_mode='compact').recruit(pairs)
    assert classified == [(item.molecule_id, item.locus_id) for item in chunk.evidence]
    assert len(chunk.evidence) == 23
    assert chunk.ambiguous_pairs == 2
    assert chunk.locus_tests < len(pairs)*3


@pytest.mark.parametrize('audit_mode', ['full', 'compact'])
def test_pruned_recruitment_agrees_across_worker_counts(audit_mode):
    from dataclasses import asdict
    from mlvamaps.short_read_recruitment import recruit_short_reads
    templates = panel()
    templates['L01'] = replace(templates['L00'], locus=templates['L01'].locus)
    pairs = [paired_reads(t, name) for name, t in templates.items()]
    serial = list(recruit_short_reads(pairs, templates, 's', threads=1, chunk_size=4,
                                     audit_mode=audit_mode))
    parallel = list(recruit_short_reads(pairs, templates, 's', threads=2, chunk_size=4,
                                       audit_mode=audit_mode))
    assert [asdict(chunk) for chunk in serial] == [asdict(chunk) for chunk in parallel]
    assert sum(chunk.skipped_locus_tests for chunk in serial) > 0
