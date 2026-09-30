"""Sample-only locus rescue: recruitment, graph geometry and shared PCR."""
from dataclasses import replace
import csv
import json
import random

import pytest

from mlvamaps.locus_reconstruction import locus_templates, run_reconstructed_fastq_inference
from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.sequence import revcomp
from mlvamaps.short_read_recruitment import ShortReadRecruiter
from mlvamaps.sample_reconstruction import SampleRecruiter, learn_template, repeat_runs
from mlvamaps.targeted_reconstruction import recover_locus


def pair(sequence, mate=None, name='p'):
    return ReadPair(name, ReadRecord(name+'/1', sequence, 'I'*len(sequence)),
                    ReadRecord(name+'/2', revcomp(mate), 'I'*len(mate)) if mate else None)


def fixture(count=8):
    rng = random.Random(951)
    def dna(n):
        return ''.join(rng.choices('ACGT', k=n))
    left, motif, right = dna(180), 'ATGACCGTA', dna(180)
    sequence = left + motif*count + right
    locus = Locus('L_9bp_432bp_8U', forward_primer=left[:20], reverse_primer=revcomp(right[-20:]),
                  repeat_unit_length_bp=9, repeat_motif='N'*9,
                  expected_product_size_bp=432, nominal_repeat_units=8)
    template = locus_templates([locus])[locus.locus_id]
    pairs = [pair(sequence[:120], sequence[70:210], f'left{i}') for i in range(3)]
    pairs += [pair(sequence[-140:], sequence[-120:], f'right{i}') for i in range(3)]
    # Right primer anchors a mate that reaches back into the repeat.
    pairs += [pair(sequence[-210:-70], sequence[-120:], f'rmate{i}') for i in range(3)]
    return locus, template, sequence, pairs


@pytest.mark.parametrize('unit,nominal,product_bp,count,depth', [(12, 10, 314, 11, 2482), (9, 16, 229, 17, 14642)])
def test_high_depth_primer_only_pileup_reaches_assembly_caller(tmp_path, unit, nominal, product_bp, count, depth):
    from mlvamaps.assembly_call import run_assembly_call
    from mlvamaps.io import read_profiles, write_tsv
    from dataclasses import asdict
    rng = random.Random(123)
    def dna(n):
        return ''.join(rng.choices('ACGT', k=n))
    nonrepeat = product_bp-unit*nominal
    left, right, motif = dna(nonrepeat//2), dna(nonrepeat-nonrepeat//2), dna(unit)
    sequence = left+motif*count+right
    locus = Locus(f'L_{unit}bp_{product_bp}bp_{nominal}U', forward_primer=left[:20],
        reverse_primer=revcomp(right[-20:]), repeat_unit_length_bp=unit,
        repeat_motif='N'*unit, expected_product_size_bp=product_bp, nominal_repeat_units=nominal)
    world = dna(80)+sequence+dna(80)
    pairs = []
    for i in range(depth):
        start = rng.randrange(len(world)-200+1)
        read = list(world[start:start+200])
        if i % 3:
            position = rng.randrange(len(read))
            read[position] = rng.choice([base for base in 'ACGT' if base != read[position]])
        read = ''.join(read)
        pairs.append(pair(read, read, str(i)))
    template = locus_templates([locus])[locus.locus_id]
    recruited = ShortReadRecruiter({locus.locus_id: template}, 's').recruit(pairs).evidence
    assert len({seq for item in recruited for seq in item.sequences}) > 256
    assert not any(item.product_sequence for item in recruited)
    calls, _, _, paths = run_reconstructed_fastq_inference(reads1=None, reads2=None,
        pairs=iter(pairs), loci=[locus], database_path=None, outdir=tmp_path/'sr',
        sample_id='s', technology='illumina', minimum_molecules=2, minimum_probability=.8)
    assert calls[0]['repeat_count'] == count
    assert calls[0]['inference_method'] == 'RECONSTRUCTED'
    assert calls[0]['status'] == 'called'
    metadata = json.loads(paths['reconstruction_metadata'].read_text())
    assert metadata['performance']['recruitment']['pileup_backend'] == 'pysam_htslib'
    audit = read_profiles(paths['reconstruction_pcr'])[0]
    assert audit['measurement_source'] == 'sassy_assembly_pcr'
    assert int(audit['pcr_product_size_bp']) == len(sequence)
    assert int(audit['input_molecules']) >= depth//2
    panel, truth = tmp_path/'panel.tsv', tmp_path/'truth.fa'
    write_tsv([asdict(locus)], panel, list(asdict(locus)))
    truth.write_text('>truth\n'+sequence+'\n')
    assembly = run_assembly_call(str(truth), str(panel), str(tmp_path/'assembly'), 's')
    assert float(read_profiles(assembly['calls'])[0]['repeat_count']) == calls[0]['repeat_count']


def test_high_depth_pileup_cannot_bridge_an_unobserved_repeat():
    from mlvamaps.short_read_evidence import MoleculeEvidence
    from mlvamaps.targeted_reconstruction import microassemble
    locus, template, sequence, _ = fixture(40)
    rng = random.Random(14)
    items = []
    for i in range(1800):
        read = list(sequence[:220] if i % 2 else sequence[-220:])
        if i % 3:
            position = rng.randrange(len(read))
            read[position] = rng.choice([base for base in 'ACGT' if base != read[position]])
        items.append(MoleculeEvidence(str(i), locus.locus_id, ('UNINFORMATIVE',),
                                     (''.join(read),), ('+',)))
    assert len({item.sequences[0] for item in items}) > 256
    product, _, reason = microassemble(items, template)
    assert not product
    assert reason in {'no_complete_path', 'ambiguous_overlap_offsets'}


def test_internal_reads_rescue_a_primer_only_locus_and_reach_sassy(tmp_path, capsys):
    locus, template, sequence, pairs = fixture()
    recruited = ShortReadRecruiter({'L': template}, 's').recruit(pairs).evidence
    assert not recover_locus(recruited, template, 's').products
    bridges = [pair(sequence[140:292], name=f'bridge{i}') for i in range(3)]
    assert not ShortReadRecruiter({'L': template}, 's').recruit(bridges).evidence
    calls, _, _, paths = run_reconstructed_fastq_inference(reads1=None, reads2=None,
        pairs=iter(pairs+bridges), loci=[locus], database_path=None, outdir=tmp_path,
        sample_id='s', technology='illumina', minimum_molecules=2, minimum_probability=.8,
        show_progress=True)
    assert calls[0]['repeat_count'] == 8
    with paths['reconstruction_pcr'].open() as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    assert row['measurement_source'] == 'sassy_assembly_pcr'
    metadata = json.loads(paths['sample_repeat_graphs'].read_text())
    assert metadata['rounds'][0]['recruited'] == 3
    assert metadata['loci'][locus.locus_id]['graph_ready']
    progress = capsys.readouterr().out
    assert 'Sample recruitment round 1 finished' in progress
    assert f'Sample anchor learning {locus.locus_id} started: round 2; 12 molecules' in progress
    assert f'Final sample template {locus.locus_id} finished' in progress


def test_repeat_graph_learns_both_arms_without_reference_sequences():
    locus, template, sequence, pairs = fixture(30)
    items = ShortReadRecruiter({'L': template}, 's').recruit(pairs).evidence
    learned, info = learn_template(items, template)
    assert info['graph_ready']
    assert not learned.primer_only
    assert learned.left.startswith(locus.forward_primer)
    assert learned.right.endswith(revcomp(locus.reverse_primer))
    assert learned.sequence((len(sequence)-len(learned.left)-len(learned.right))/9) == sequence


@pytest.mark.parametrize('insert_known', [False, True])
def test_partial_graph_evidence_estimates_count_only_with_length_information(tmp_path, insert_known):
    locus, _, sequence, _ = fixture(30)
    pairs = [pair(sequence[:210], sequence[-210:], f'p{i}') for i in range(30)]
    calls, _, _, paths = run_reconstructed_fastq_inference(reads1=None, reads2=None,
        pairs=iter(pairs), loci=[locus], database_path=None, outdir=tmp_path,
        sample_id='s', technology='illumina', minimum_molecules=2, minimum_probability=.8,
        **({'insert_mean': len(sequence), 'insert_sd': 3} if insert_known else {}))
    if insert_known:
        assert calls[0]['repeat_count'] == 30
        assert calls[0]['inference_method'] == 'INFERRED'
        assert calls[0]['repeat_count_min'] <= 30 <= calls[0]['repeat_count_max']
    else:
        assert 0 < calls[0]['repeat_count'] < 30
        assert calls[0]['status'] == 'ambiguous'
        assert 'repeat_length_lower_bound' in calls[0]['reason']
        assert calls[0]['repeat_count_min'] <= 30 <= calls[0]['repeat_count_max']
    assert 'reconstruction_pcr' not in paths  # No fabricated contig enters PCR.


def test_sample_index_competes_across_loci_and_rejects_repeat_only_reads():
    locus, template, sequence, pairs = fixture()
    items = ShortReadRecruiter({'L': template}, 's').recruit(pairs).evidence
    competitor = replace(template, locus=replace(locus, locus_id='other'))
    other = [replace(item, locus_id='other', template=competitor) for item in items]
    recruiter = SampleRecruiter({locus.locus_id: items, 'other': other},
                                {locus.locus_id: template, 'other': competitor})
    item, ambiguous = recruiter.recruit_pair(pair(sequence[40:140]))
    assert item is None and set(ambiguous) == {locus.locus_id, 'other'}
    item, _ = recruiter.recruit_pair(pair('ATGACCGTA'*20))
    assert item is None


@pytest.mark.parametrize('unit', [4, 9, 21, 45])
def test_motif_discovery_handles_long_and_imperfect_units(unit):
    rng = random.Random(unit)
    motif = ''.join(rng.choices('ACGT', k=unit))
    repeat = list(motif*5)
    if unit >= 10:
        repeat[unit+3] = next(b for b in 'ACGT' if b != repeat[unit+3])
    runs = repeat_runs(''.join(repeat), unit)
    expected = min(motif[i:]+motif[:i] for i in range(unit))
    assert any(found == expected for _, _, found in runs)


def test_sample_recruitment_preserves_mates_and_filters_low_quality():
    locus, template, sequence, pairs = fixture()
    items = ShortReadRecruiter({'L': template}, 's').recruit(pairs).evidence
    recruiter = SampleRecruiter({locus.locus_id: items}, {locus.locus_id: template})
    candidate = pair(sequence[40:140], 'ATGACCGTA'*12)
    item, _ = recruiter.recruit_pair(candidate)
    assert item and len(item.sequences) == 2 and item.orientations == ('+', '-')
    poor = replace(candidate, read1=replace(candidate.read1, quality='!'*100))
    assert recruiter.recruit_pair(poor)[0] is None


def test_two_round_rescue_is_order_and_worker_independent(tmp_path):
    rng = random.Random(3491)
    def dna(n):
        return ''.join(rng.choices('ACGT', k=n))
    left, right, motif = dna(350), dna(350), 'ATGACCGTA'
    sequence = left+motif*8+right
    locus = Locus('L', forward_primer=left[:20], reverse_primer=revcomp(right[-20:]),
                  repeat_unit_length_bp=9, expected_product_size_bp=len(sequence), nominal_repeat_units=8)
    reads = []
    for i in range(3):
        reads += [pair(sequence[:120], sequence[70:210], f'left{i}'),
                  pair(sequence[-210:-70], sequence[-120:], f'right{i}'),
                  pair(sequence[170:350], name=f'extendleft{i}'),
                  pair(sequence[422:602], name=f'extendright{i}'),
                  pair(sequence[300:480], name=f'bridge{i}')]
    reads += [pair(dna(150), name=f'background{i}') for i in range(600)]
    results = []
    for threads, ordered in ((1, reads), (2, list(reversed(reads)))):
        calls, _, _, paths = run_reconstructed_fastq_inference(reads1=None, reads2=None,
            pairs=iter(ordered), loci=[locus], database_path=None, outdir=tmp_path/str(threads),
            sample_id='s', technology='illumina', minimum_molecules=2, minimum_probability=.8,
            threads=threads)
        assert calls[0]['repeat_count'] == 8
        metadata = json.loads(paths['sample_repeat_graphs'].read_text())
        assert [row['recruited'] for row in metadata['rounds']] == [6, 3]
        results.append(calls)
    assert results[0] == results[1]


def test_replay_qc_preserves_orphan_ids_and_does_not_double_count(tmp_path):
    from mlvamaps.io import write_fastq
    from mlvamaps.short_read_qc import filtered_pairs, replay_filtered_pairs
    first, second = tmp_path/'r1.fq', tmp_path/'r2.fq'
    write_fastq([ReadRecord('p/1', 'ACGT'*20, '!'*80)], first)
    write_fastq([ReadRecord('p/2', 'ACGT'*20, 'I'*80)], second)
    counters = {}
    options = dict(min_length=20, min_mean_quality=20, trim_quality=0, min_pair_retention=.5)
    original = list(filtered_pairs(first, second, tmp_path/'f1.gz', tmp_path/'f2.gz', tmp_path/'orphans.gz',
        materialize=False, counters=counters, statistics={}, sample_id='s', **options))
    replayed = list(replay_filtered_pairs(first, second, **options))
    assert original == replayed and original[0].molecule_id == 'p/2'
    assert counters['input_pairs'] == 1
