"""Scientific acceptance cases for reconstruction-first paired-end calling."""
from dataclasses import replace
import random

import pytest

from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.sequence import revcomp
from mlvamaps.short_read_evidence import classify_pair, panel_template, repeat_fraction, InsertDistribution
from mlvamaps.short_read_recruitment import ShortReadRecruiter
from mlvamaps.targeted_reconstruction import recover_locus, microassemble
from mlvamaps.locus_products import genotype_product, LocusProduct


def test_bounded_overlap_seed_search_matches_exhaustive_offsets():
    from mlvamaps.targeted_reconstruction import _overlaps
    from mlvamaps.short_read_evidence import repetitive
    rng = random.Random(2981)
    for size in (20, 49, 50, 99, 100, 149, 150, 200):
        overlap = ''.join(rng.choices('ACGT', k=size))
        for errors in (0, int(size*.02), int(size*.02)+1):
            changed = list(overlap)
            # Damage consecutive seed blocks to exercise the last required one.
            for offset in range(errors):
                position = offset*15
                changed[position] = next(b for b in 'ACGT' if b != changed[position])
            left, right = 'ACTG'*10+overlap, ''.join(changed)+'GTCA'*10
            expected = tuple(n for n in range(20, min(len(left), len(right))+1)
                if sum(a != b for a, b in zip(left[-n:], right[:n])) <= int(n*.02)
                and not repetitive(left[-n:], 'AATG') and not repetitive(right[:n], 'AATG'))
            assert _overlaps(left, right, 'AATG', unit=4, mismatch_fraction=.02) == expected


def test_native_overlap_end_masks_match_scalar_at_every_length():
    from mlvamaps.targeted_reconstruction import _nonrepeat_overlap_sizes
    from mlvamaps.short_read_evidence import repetitive
    rng = random.Random(843)
    cases = [('', ''), ('AAAA', ''), ('éééAééé', 'é'), ('N'*60, 'N')]
    for unit in (1, 4, 9, 39, 150, 1800):
        motif = ''.join(rng.choices('ACGT', k=unit))
        sequence = (motif*200)[3:153]
        cases += [(sequence, motif), ('ACGT'+sequence+'TGCA', motif)]
        mutated = list(sequence)
        for position in range(0, len(mutated), 10):
            mutated[position] = next(b for b in 'ACGT' if b != mutated[position])
        cases.append((''.join(mutated), motif))
    for sequence, motif in cases:
        for suffix in (True, False):
            expected = frozenset(n for n in range(1, len(sequence)+1)
                if not repetitive(sequence[-n:] if suffix else sequence[:n], motif))
            assert _nonrepeat_overlap_sizes(sequence, motif, suffix) == expected


@pytest.fixture
def template():
    rng = random.Random(761)
    def dna(n):
        return ''.join(rng.choices('ACGT', k=n))
    return panel_template(Locus('L', forward_primer=dna(20), reverse_primer=dna(20),
        left_flank_sequence=dna(80), right_flank_sequence=dna(80),
        repeat_motif='AATG', repeat_unit_length_bp=4, expected_max_repeats=12))


def pair(first, second=None, name='p', quality='I'):
    return ReadPair(name, ReadRecord(name+'/1', first, quality*len(first)),
                    ReadRecord(name+'/2', second, quality*len(second)) if second else None)


def evidence(template, sequences):
    return [classify_pair(pair(a, b, str(i)), template) for i,(a,b) in enumerate(sequences)]


@pytest.mark.parametrize('count', [2, 8, 40, 140])
def test_direct_novel_allele_bypasses_candidate_ceiling_and_inference(template, count, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('direct evidence must bypass assembly and inference')
    monkeypatch.setattr('mlvamaps.targeted_reconstruction.microassemble', forbidden)
    monkeypatch.setattr('mlvamaps.targeted_reconstruction.candidate_likelihood', forbidden)
    seq = template.sequence(count)
    result = recover_locus(evidence(template, [(seq, revcomp(seq))]*3), template, 's', maximum=100)
    assert result.method == 'DIRECT'
    assert genotype_product(result.products[0], template.locus).repeat_count == count
    assert result.counts['FULL_SPAN'] == 3


def test_unique_overlapping_pair_is_direct(template, monkeypatch):
    sequence = template.sequence(8)
    # Overlap includes both repeat boundaries, anchoring the offset uniquely.
    items = evidence(template, [(sequence[:160], revcomp(sequence[70:]))]*3)
    assert all(not e.product_sequence for e in items)
    result = recover_locus(items, template, 's')
    assert result.method == 'DIRECT'
    assert result.products[0].sequence == sequence
    assert result.counts['FULL_SPAN'] == 3


def test_targeted_three_read_stitch_preserves_repeat_and_flank_snps(template):
    sequence = template.sequence(8)
    sequence = sequence[:60] + ('A' if sequence[60] != 'A' else 'C') + sequence[61:]
    sequence = sequence[:111] + 'C' + sequence[112:]
    items = evidence(template, [(sequence[:100], None), (sequence[65:165], None), (sequence[140:], None)]*3)
    result = recover_locus(items, template, 's')
    assert result.method == 'RECONSTRUCTED'
    assert result.products[0].sequence == sequence
    expected = genotype_product(LocusProduct('s', 'L', 'assembly', 'a', sequence), template.locus)
    assert genotype_product(result.products[0], template.locus).combined_marker == expected.combined_marker


def test_reconstruction_does_not_merge_different_repeat_lengths(template):
    sequence = template.sequence(12)
    longer = template.sequence(13)
    items = evidence(template, [(sequence[:80], None), (sequence[50:190], None),
                                (longer[50:194], None), (sequence[160:], None)]*3)
    product, _, reason = microassemble(items, template)
    assert not product
    assert reason == 'multiple_reconstructions'


@pytest.mark.parametrize('count', [8, 12])
def test_disconnected_reconstructions_preserve_agreed_length(template, count, tmp_path):
    from mlvamaps.short_read_evidence import MoleculeEvidence
    sequence = template.sequence(8)
    alternative = list(template.sequence(count))
    # Two independently assembled paths with divergent interiors cannot form
    # overlap edges between them, but can still agree on primer-product length.
    for position in range(25, len(alternative)-25, 7):
        alternative[position] = next(b for b in 'ACGT' if b != alternative[position])
    alternative = ''.join(alternative)
    items = []
    for haplotype in (sequence, alternative):
        for read in (haplotype[:120], haplotype[60:190], haplotype[140:]):
            items.append(MoleculeEvidence(str(len(items)), 'L', ('UNINFORMATIVE',), (read,), ('+',)))
    product, members, reason = microassemble(items, template)
    if count == 8:
        assert product == ''.join(a if a == b else 'N' for a, b in zip(sequence, alternative))
        assert len(members) == 6 and not reason
        assert genotype_product(LocusProduct('s', 'L', 'short_read', 'v', product), template.locus).repeat_count == 8
        from dataclasses import asdict
        from mlvamaps.io import write_fastq, write_tsv, read_profiles
        from mlvamaps.short_reads import run_short_read_call
        locus = replace(template.locus, left_flank_sequence='', right_flank_sequence='',
                        repeat_motif='NNNN', expected_product_size_bp=len(sequence), nominal_repeat_units=8)
        panel = tmp_path/'panel.tsv'
        write_tsv([asdict(locus)], panel, list(asdict(locus)))
        pairs = [pair(read, revcomp(mate), str(i)) for i, (read, mate) in enumerate(
            (reads for haplotype in (sequence, alternative) for reads in
             ((haplotype[:120], haplotype[60:190]), (haplotype[140:], haplotype[140:]))))]
        first, second = tmp_path/'r1.fq', tmp_path/'r2.fq'
        write_fastq([p.read1 for p in pairs], first)
        write_fastq([p.read2 for p in pairs], second)
        output = run_short_read_call(str(first), str(second), str(panel), str(tmp_path/'out'),
                                     's', show_progress=False)
        call = read_profiles(output['calls'])[0]
        assert float(call['repeat_count']) == 8
        assert call['call_method'] == 'RECONSTRUCTED'
        assert '8 repeats' in output['report'].read_text()
    else:
        assert not product and reason == 'multiple_reconstructions'


@pytest.mark.parametrize('count', [30, 80, 160])
def test_repeat_only_overlap_never_selects_arbitrary_copy_number(template, count):
    sequence = template.sequence(count)
    items = evidence(template, [(sequence[:150], revcomp(sequence[-150:]))]*4)
    result = recover_locus(items, template, 's', maximum=200)
    assert result.method == 'AMBIGUOUS'
    assert not result.products
    assert result.counts['FLANK_PAIR'] == 4


def test_anchored_unmapped_mate_and_imperfect_cyclic_motif(template):
    motif = ('ATGA'*25)[:43]+'C'+('ATGA'*25)[43:]
    assert repeat_fraction(motif, template.motif) > .9
    assert repeat_fraction(revcomp(motif), template.motif) > .9
    item = classify_pair(pair(template.left, revcomp(motif)), template)
    assert {'ANCHORED_REPEAT', 'REPEAT_RICH'} <= set(item.classes)
    assert len(item.sequences) == 2
    result = recover_locus([item]*4, template, 's')
    assert not result.identifiable


@pytest.mark.parametrize('side', ['LEFT', 'RIGHT'])
def test_one_boundary_and_softclip_remain_ambiguous(template, side):
    seq = template.left + template.motif*20 if side == 'LEFT' else template.motif*20 + template.right
    item = classify_pair(pair(seq), template)
    assert side+'_BOUNDARY' in item.classes
    assert 'SOFTCLIP_'+side in item.classes
    result = recover_locus([item]*30, template, 's')
    assert result.method == 'AMBIGUOUS'
    assert result.confidence < .8


def test_repeat_fraction_is_configurable(template):
    seq = template.motif*20+'CGCTCTCGCCGCCGCGTCGC'
    low = classify_pair(pair(template.left, revcomp(seq)), template, .6)
    high = classify_pair(pair(template.left, revcomp(seq)), template, .99)
    assert 'REPEAT_RICH' in low.classes
    assert 'REPEAT_RICH' not in high.classes


def test_multiple_contexts_of_same_locus_do_not_compete_as_separate_loci(template):
    other = replace(template, left=template.left[:45]+'A'+template.left[46:])
    seq = other.sequence(8)
    chunk = ShortReadRecruiter({'a': template, 'b': other}, 's').recruit([pair(seq)])
    assert len(chunk.evidence) == 1
    assert chunk.ambiguous_pairs == 0
    assert chunk.evidence[0].template == other
    # A genuinely indistinguishable different locus stays ambiguous.
    identical = replace(other, locus=replace(other.locus, locus_id='other'))
    chunk = ShortReadRecruiter({'a': other, 'b': identical}, 's').recruit([pair(seq)])
    assert not chunk.evidence and chunk.ambiguous_pairs == 1


def test_low_quality_does_not_supply_direct_sequence(template):
    sequence = template.sequence(8)
    item = classify_pair(pair(sequence, quality='!'), template)
    assert not recover_locus([item], template, 's').products


def test_high_depth_direct_path_never_runs_assembly(template, monkeypatch):
    def forbidden(*args):
        pytest.fail('assembly on complete molecules')
    monkeypatch.setattr('mlvamaps.targeted_reconstruction.microassemble', forbidden)
    seq = template.sequence(8)
    result = recover_locus(evidence(template, [(seq, None)]*1000), template, 's')
    assert result.products[0].support_count == 1000
    assert result.confidence > .99


def test_bounded_graph_reports_limit(template):
    items = evidence(template, [(template.sequence(8), None)])
    assert microassemble(items, template, max_nodes=0)[2] == 'assembly_node_limit'


@pytest.mark.parametrize('primer_only', [False, True])
@pytest.mark.parametrize('reason', ['assembly_node_limit', 'assembly_path_limit'])
def test_assembly_limits_survive_likelihood_fallback(template, primer_only, reason):
    items = evidence(template, [(template.left, None)])
    result = recover_locus(items, replace(template, primer_only=primer_only), 's',
                           precomputed=([], ('', [], reason)))
    assert result.limit_reached
    assert reason in result.reason
    assert not result.products and result.best is None


def test_contained_flank_reads_do_not_vote_for_uncovered_bases(template):
    sequence = template.sequence(8)
    first = sequence[:90]+next(base for base in 'ACGT' if base != sequence[90])+sequence[91:120]
    items = evidence(template, [(first, None)] + [(sequence[60:190], None)]*3
                     + [(sequence[140:], None)]*5 + [(sequence[:60], None)]*1000)
    product, members, reason = microassemble(items, template)
    assert product == sequence
    assert len(members) == len(items)
    assert not reason


def test_htslib_pileup_counts_conflicting_mates_once_without_depth_truncation():
    from mlvamaps.targeted_reconstruction import _pileup_consensus
    # Exceed HTSlib's default 8,000-read cap. Conflicting observations consume
    # one N vote, not a second vote for either base, including at exactly 70%.
    rows = [('A'*100, 0, set(range(12000))), ('C'*100, 0, set(range(4000)))]
    assert _pileup_consensus(rows, 100) == 'N'*100
    rows[1] = ('C'*100, 0, set(range(3600)))
    assert _pileup_consensus(rows, 100) == 'A'*100


def test_htslib_pileup_preserves_uncovered_bases_and_overlapping_molecule_votes():
    from mlvamaps.targeted_reconstruction import _pileup_consensus
    rows = [('ACGT', 2, {'a', 'b'}), ('CGTC', 3, {'a'}), ('TTTT', 2, {'c'})]
    assert _pileup_consensus(rows, 9) == 'NNNNNTCNN'
    assert _pileup_consensus(reversed(rows), 9) == 'NNNNNTCNN'


def test_native_batch_preserves_threshold_and_unknown_votes(monkeypatch):
    from mlvamaps.targeted_reconstruction import _pileup_consensus, _pileup_consensuses
    import pysam
    calls = []
    original = pysam.samtools.consensus
    def consensus(*args, **kwargs):
        calls.append(args)
        return original(*args, **kwargs)
    monkeypatch.setattr(pysam.samtools, 'consensus', consensus)
    pools, expected = [], []
    for bases, base in [('N', 'N'), ('AAAAAAACCC', 'A'), ('AAAAAACCCC', 'N'),
                        ('AAAAAAANNN', 'A'), ('AAAAAANNNN', 'N')]:
        rows = [(b*4, 2, {i}) for i, b in enumerate(bases)]
        pools.append((rows, 8))
        expected.append('NN'+base*4+'NN')
    pools.append(([('AAAA', 2, set(range(8000))), ('CCCC', 2, set(range(8000, 12001)))], 8))
    expected.append('N'*8)  # Native consensus must include reads beyond depth 8,000.
    assert _pileup_consensuses(pools) == expected
    assert [_pileup_consensus(rows, length) for rows, length in pools] == expected
    assert len(calls) == 1


def test_pileup_preserves_trim_coordinates_support_and_ambiguous_seeds():
    from mlvamaps.targeted_reconstruction import _pileup_sequences
    rng = random.Random(918)
    world = ''.join(rng.choices('ACGT', k=400))
    sequences, expected = {}, {}
    for start in range(0, 100, 10):
        for length in (80, 100, 150, 200):
            truth = world[start:start+length]
            name = f'{start}:{length}'
            members = {name+':a', name+':b', name+':c'}
            sequences[truth] = members
            position = length//2
            error = truth[:position]+next(b for b in 'ACGT' if b != truth[position])+truth[position+1:]
            sequences[error] = {name+':error'}
            expected[truth] = members | sequences[error]
    # Repeated seeds cannot establish a unique placement, even at the same
    # length and coordinate. Preserve these distinct observations.
    for sequence in ('A'*100, 'A'*99+'C'):
        sequences[sequence] = {sequence}
        expected[sequence] = {sequence}
    assert _pileup_sequences(sequences) == expected


def test_trimmed_pileup_support_cannot_vote_outside_observed_bases():
    from mlvamaps.targeted_reconstruction import _pileup_consensus, _pileup_sequences
    rng = random.Random(189)
    truth = ''.join(rng.choices('ACGT', k=200))
    error = truth[:190]+next(b for b in 'ACGT' if b != truth[190])+truth[191:]
    sequences = {truth: {'a', 'b'}, error: {'c'}, truth[20:100]: set(range(1000))}
    rows = {}
    result = _pileup_sequences(sequences, rows)
    consensus = truth[:190]+'N'+truth[191:]
    assert list(result) == [consensus]
    assert _pileup_consensus(rows[consensus]+[(truth[100:], 100, {'extra'})], 200) == truth


def test_assembly_node_limit_applies_after_containment(template):
    from mlvamaps.short_read_evidence import MoleculeEvidence
    sequence = template.sequence(8)
    items = [MoleculeEvidence(str(i), 'L', ('UNINFORMATIVE',), (sequence[:end],), ('+',))
             for i, end in enumerate(range(30, len(sequence)+1))]
    assert microassemble(items, template, max_nodes=1)[0] == sequence
    items.append(MoleculeEvidence('separate', 'L', ('UNINFORMATIVE',), ('G'*100,), ('+',)))
    assert microassemble(items, template, max_nodes=1)[2] == 'assembly_node_limit'


def test_insert_likelihood_distinguishes_or_withholds_adjacent_lengths(template):
    items = evidence(template, [(template.left, revcomp(template.right))]*40)
    confident = recover_locus(items, template, 's', InsertDistribution(320, 3, 40, 'override'))
    uncertain = recover_locus(items, template, 's', InsertDistribution(320, 100, 40, 'override'))
    assert confident.method == 'INFERRED'
    assert genotype_product(confident.products[0], template.locus).repeat_count == 30
    assert uncertain.method == 'AMBIGUOUS'
    assert uncertain.confidence < confident.confidence


def test_no_evidence_is_no_call(template):
    result = recover_locus([], template, 's')
    assert result.method == 'NO_CALL' and not result.products and result.confidence == 0


def test_short_reads_have_no_engine_selector():
    from mlvamaps.cli import build_parser
    with pytest.raises(SystemExit):
        build_parser().parse_args(['call', '--sr-engine', 'competitive'])


def test_full_product_error_consensus_preserves_supported_snp(template):
    truth = template.sequence(8)
    truth = truth[:60]+'A'+truth[61:]
    reads = []
    for i in range(12):
        position = 65+i
        base = 'A' if truth[position] != 'A' else 'C'
        reads.append((truth[:position]+base+truth[position+1:], None))
    result = recover_locus(evidence(template, reads), template, 's')
    assert result.products[0].sequence == truth
    assert result.method == 'DIRECT'


def test_conflicting_singleton_lengths_retain_estimate_and_alternative(template):
    result = recover_locus(evidence(template, [(template.sequence(8), None),
                                               (template.sequence(12), None)]), template, 's')
    assert result.method == 'MIXED'
    assert len(result.products) == 2
    assert result.products[0].support_count == result.products[1].support_count == 1
    from mlvamaps.locus_reconstruction import _common_call
    call = _common_call(template.locus, result.products, 2, 's', 'illumina', 3, .8, result)
    assert call['status'] == 'mixed'
    assert call['repeat_count'] == 8
    assert call['second_best_repeat_count'] == 12
    assert (call['repeat_count_min'], call['repeat_count_max']) == (8, 12)


def test_single_flank_pair_can_infer_length_and_retain_observed_bases(template):
    items = evidence(template, [(template.left[-40:], revcomp(template.right[:40]))])
    assert items[0].fragment_offset is not None
    insert = InsertDistribution(items[0].fragment_offset + 8*template.unit, .1, 0, 'override')
    result = recover_locus(items, template, 's', insert=insert, minimum_spanning_pairs=100)
    assert result.method == 'INFERRED' and result.best == 8
    assert result.counts['FLANK_PAIR'] == 1
    assert template.left[-40:] in result.products[0].sequence
    assert template.right[:40] in result.products[0].sequence


@pytest.mark.parametrize('kind', ['insert', 'boundary', 'presence'])
def test_uncertain_length_estimates_survive_into_calls(template, kind):
    from mlvamaps.locus_reconstruction import _common_call
    from mlvamaps.unified_fastq import common_calls_to_compatibility
    from mlvamaps.profile_matching import build_fingerprint
    sequences = [(template.left, revcomp(template.right))] if kind == 'insert' else [
        (template.left + (template.motif*20 if kind == 'boundary' else ''), None)]
    result = recover_locus(evidence(template, sequences*4), template, 's',
                           InsertDistribution(320, 100, 40, 'override') if kind == 'insert' else None)
    assert not result.products
    common = _common_call(template.locus, result.products, 4, 's', 'illumina', 3, .8, result)
    call = common_calls_to_compatibility([common])[0]
    if kind == 'presence':
        assert call['repeat_count'] == ''
        assert call['status'] == 'PRESENT_COUNT_UNKNOWN'
        assert common['repeat_count_min'] == ''
        assert 'no_repeat_length_information' in call['evidence']
    else:
        assert call['status'] == 'AMBIGUOUS'
        assert call['repeat_count'] == (30 if kind == 'insert' else 20)
        assert common['repeat_count_min'] <= call['repeat_count'] <= common['repeat_count_max']
        assert call['allele_distribution']
        if kind == 'boundary':
            assert 'repeat_length_lower_bound' in call['evidence']
        fingerprint, _ = build_fingerprint('s', [{'locus_id': 'L', 'called_repeat_count': call['repeat_count']}], [template.locus])
        assert fingerprint[0]['L'] == call['repeat_count']


@pytest.mark.parametrize('primer_only', [False, True])
def test_reconstruction_never_loads_database_sequences(tmp_path, template, monkeypatch, primer_only):
    from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
    import json
    locus = template.locus
    if primer_only:
        locus = replace(locus, left_flank_sequence='', right_flank_sequence='', repeat_motif='NNNN',
                        expected_product_size_bp=len(template.sequence(8)), nominal_repeat_units=8)
    def forbidden(*args, **kwargs):
        pytest.fail('Illumina recovery accessed database sequences')
    monkeypatch.setattr('mlvamaps.candidate_contexts._base_contexts', forbidden)
    results = []
    for index, database in enumerate((None, tmp_path/'absent_database')):
        result = run_reconstructed_fastq_inference(reads1=None, reads2=None,
            pairs=iter([pair(template.sequence(12), name=str(i)) for i in range(3)]),
            loci=[locus], database_path=database, outdir=tmp_path/str(index), sample_id='s',
            technology='illumina', minimum_molecules=2, minimum_probability=.8)
        results.append(result[0])
        metadata = json.loads(result[3]['reconstruction_metadata'].read_text())
        assert metadata['performance']['recruitment']['template_source'] == 'locus_panel'
        assert metadata['performance']['recruitment']['template_contexts'] == 1
    assert results[0] == results[1]
    assert results[0][0]['repeat_count'] == 12


@pytest.mark.parametrize('scenario', ['span', 'overlap', 'reconstruct', 'repeat_only', 'unlinked'])
def test_legacy_primer_only_recovery(template, scenario, monkeypatch):
    from mlvamaps.locus_reconstruction import locus_templates
    from mlvamaps.primers import _locus_from_primer_row
    locus = _locus_from_primer_row('L_4bp_232bp_8U', template.locus.forward_primer,
                                   template.locus.reverse_primer)
    target = locus_templates([locus])[locus.locus_id]
    sequence = template.sequence(12)
    cases = {
        'span': [(sequence, None)],
        'overlap': [(sequence[:170], revcomp(sequence[70:]))],
        # Only the first mate in each pair has a primer. Rescued mates bridge
        # the gap using observed unique overlaps outside the repeat tract.
        'reconstruct': [(sequence[:80], revcomp(sequence[50:160])),
                        (revcomp(sequence[180:]), sequence[135:205])],
        'repeat_only': [(template.sequence(80)[:150], revcomp(template.sequence(80)[-150:]))],
        'unlinked': [(sequence[:70], revcomp(sequence[-70:]))],
    }
    def forbidden(*args, **kwargs):
        pytest.fail('Primer-only reads must not synthesize a candidate repeat interval')
    monkeypatch.setattr('mlvamaps.targeted_reconstruction.candidate_likelihood', forbidden)
    pairs = [pair(a, b, str(i)) for i, (a, b) in enumerate(cases[scenario]*3)]
    chunk = ShortReadRecruiter({locus.locus_id: target}, 's', audit_mode='compact').recruit(pairs)
    assert len(chunk.evidence) == len(pairs)
    insert = InsertDistribution(300, 10, 0, 'override') if scenario in ('repeat_only', 'unlinked') else None
    result = recover_locus(chunk.evidence, target, 's', insert=insert)
    if scenario in ('repeat_only', 'unlinked'):
        assert result.method == 'AMBIGUOUS'
        assert not result.products and not len(result.states)
        assert result.reason.startswith('primer_only_no_complete_product')
        from mlvamaps.locus_reconstruction import _common_call
        from mlvamaps.unified_fastq import common_calls_to_compatibility
        call = _common_call(locus, [], len(pairs), 's', 'illumina', 2, .8, result)
        assert result.reason in common_calls_to_compatibility([call])[0]['evidence']
    else:
        assert result.method == ('RECONSTRUCTED' if scenario == 'reconstruct' else 'DIRECT')
        assert result.products[0].sequence == sequence
        assert genotype_product(result.products[0], locus).repeat_count == 12


@pytest.mark.parametrize('duplicate_middle', [False, True])
def test_primer_only_reconstruction_uses_observed_reads_then_assembly_pcr(tmp_path, template, duplicate_middle):
    import csv
    from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
    from mlvamaps.assembly_call import legacy_assembly_call_rows, pcr_rows_to_products
    from mlvamaps.in_silico_pcr import read_pcr_results, run_in_silico_pcr_loci
    from mlvamaps.io import write_fasta
    locus = replace(template.locus, left_flank_sequence='', right_flank_sequence='',
        repeat_motif='NNNN', expected_product_size_bp=232, nominal_repeat_units=8)
    sequence = template.sequence(12)
    middle = sequence[60:190]
    # One sequencing substitution in a 50 bp overlap; no individual read or
    # read pair contains a complete primer-bounded product.
    middle = middle[:110] + next(b for b in 'ACGT' if b != middle[110]) + middle[111:]
    pairs = []
    for i in range(3):
        item = pair('TTTTT'+sequence[:120], revcomp(middle), str(i))
        pairs.append(replace(item, read1=replace(item.read1, quality='!!!!!'+'I'*120)))
    if duplicate_middle:
        # Same-coordinate reads with a substitution used to create reciprocal
        # full-length overlap edges and abort with cyclic_overlap_graph.
        altered = middle[:10] + next(b for b in 'ACGT' if b != middle[10]) + middle[11:]
        pairs += [pair(sequence[:120], revcomp(altered), f'error{i}') for i in range(3)]
    pairs += [pair(sequence[140:], name=f'right{i}') for i in range(5)]
    calls, _, _, paths = run_reconstructed_fastq_inference(reads1=None, reads2=None,
        pairs=iter(pairs), loci=[locus], database_path=None, outdir=tmp_path/'sr', sample_id='s',
        technology='illumina', minimum_molecules=2, minimum_probability=.8)
    assert calls[0]['repeat_count'] == 12
    assert calls[0]['inference_method'] == 'RECONSTRUCTED'
    with paths['reconstructed_locus_variants'].open() as handle:
        reconstructed = next(csv.DictReader(handle, delimiter='\t'))
    if duplicate_middle:
        assert reconstructed['sequence'][70] == 'N'
    assert len(reconstructed['sequence']) == len(sequence)
    with paths['reconstruction_pcr'].open() as handle:
        audit = next(csv.DictReader(handle, delimiter='\t'))
    assert audit['measurement_source'] == 'sassy_assembly_pcr'
    assert audit['pcr_product_size_bp'] == str(len(sequence))
    truth = tmp_path/'truth.fasta'
    write_fasta([('truth', sequence)], truth)
    pcr = run_in_silico_pcr_loci(truth, [locus], tmp_path/'assembly', threads=1)
    products = pcr_rows_to_products(read_pcr_results(pcr['stats'], pcr['products']), [locus], 's')
    assembly_call = legacy_assembly_call_rows([locus], products, 's')[0]
    assert calls[0]['repeat_count'] == assembly_call['repeat_count']


def test_verified_product_is_not_discarded_for_low_quality_outside_primers(template):
    sequence = template.sequence(8)
    item = pair('TTTT'+sequence+'AAAA')
    item = replace(item, read1=replace(item.read1, quality='!!!!'+'I'*len(sequence)+'!!!!'))
    observed = classify_pair(item, template)
    assert observed.product_sequence == sequence
    result = recover_locus([observed]*3, template, 's')
    assert result.method == 'DIRECT' and result.products[0].sequence == sequence


def test_short_read_primer_indel_uses_sassy_assembly_size_calibration(tmp_path, template):
    from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
    import csv
    locus = replace(template.locus, expected_product_size_bp=232, nominal_repeat_units=8)
    sequence = template.sequence(8)
    sequence = sequence[:10] + 'C' + sequence[10:]
    calls, _, _, paths = run_reconstructed_fastq_inference(reads1=None, reads2=None,
        pairs=iter([pair(sequence, name=str(i)) for i in range(3)]), loci=[locus],
        database_path=None, outdir=tmp_path, sample_id='s', technology='illumina',
        minimum_molecules=2, minimum_probability=.8)
    assert calls[0]['repeat_count'] == 8
    with paths['reconstruction_pcr'].open() as handle:
        audit = next(csv.DictReader(handle, delimiter='\t'))
    assert audit['poa_consensus_bp'] == '233'
    assert audit['pcr_product_size_bp'] == '232'



@pytest.mark.parametrize('unit', [0, 4])
def test_primer_only_without_length_calibration_preserves_sequence_but_not_count(tmp_path, template, unit):
    from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
    import csv
    locus = replace(template.locus, left_flank_sequence='', right_flank_sequence='',
                    repeat_motif='', repeat_unit_length_bp=unit)
    calls, _, _, paths = run_reconstructed_fastq_inference(reads1=None, reads2=None,
        pairs=iter([pair(template.sequence(8), name=str(i)) for i in range(3)]),
        loci=[locus], database_path=None, outdir=tmp_path, sample_id='s',
        technology='illumina', minimum_molecules=2, minimum_probability=.8)
    assert calls[0]['status'] == 'detected_unresolved'
    assert calls[0]['repeat_count'] == calls[0]['candidate_distribution'] == ''
    with paths['reconstructed_locus_variants'].open() as handle:
        product = next(csv.DictReader(handle, delimiter='\t'))
    assert product['sequence'] == template.sequence(8)
    assert product['repeat_count'] == ''
    with paths['reconstruction_pcr'].open() as handle:
        pcr = next(csv.DictReader(handle, delimiter='\t'))
    assert pcr['called_repeat_count'] == ''
    assert pcr['pcr_status'] == 'COUNT_UNCALIBRATED'


@pytest.mark.parametrize('primer,sequence', [('ACGTACGA', 'ACTTACGA'), ('ARYTACGTARYT', 'AGCTACGTAGCT')])
def test_short_and_degenerate_primer_matches_survive_seed_filter(primer, sequence, monkeypatch):
    from mlvamaps.locus_reconstruction import locus_templates
    def forbidden(*args, **kwargs):
        pytest.fail('Recruitment must not launch an external primer matcher')
    monkeypatch.setattr('mlvamaps.locus_measurement.find_anchor', forbidden)
    target = locus_templates([Locus('L', forward_primer=primer, reverse_primer='CCGTGGTCCGTT')])['L']
    reads = [pair(sequence, name=str(i)) for i in range(5)]
    expected = [classify_pair(p, target) for p in reads]
    assert all(expected)
    from mlvamaps.short_read_recruitment import recruit_short_reads
    for threads in (1, 2):
        chunks = list(recruit_short_reads(reads, {'L': target}, 's', threads=threads, chunk_size=2))
        assert [item for chunk in chunks for item in chunk.evidence] == expected


@pytest.mark.parametrize('change', ['mismatch', 'insertion', 'deletion', 'unknown'])
def test_primer_only_matching_accepts_indels_but_not_unknown_bases(template, change):
    from mlvamaps.short_read_evidence import _primer_hit
    primer = template.locus.forward_primer
    altered = {
        'mismatch': primer[:8] + next(b for b in 'ACGT' if b != primer[8]) + primer[9:],
        'insertion': primer[:8] + 'C' + primer[8:],
        'deletion': primer[:8] + primer[9:],
        'unknown': primer[:8] + 'N' + primer[9:],
    }[change]
    hit = _primer_hit('GGG' + altered + 'TTT', primer)
    if change == 'unknown':
        assert hit is None
    else:
        assert hit is not None and hit.edits == 1
        assert hit.query_start == 3 and hit.query_end == 3 + len(altered)


def test_primer_seed_prefilter_retains_matches_at_edit_limit():
    from mlvamaps.short_read_evidence import _primer_hit
    rng = random.Random(461)
    for length in (8, 12, 20, 31):
        for trial in range(30):
            primer = ''.join(rng.choices('ACGT', k=length))
            sequence = primer
            errors = min(2, max(1, length//5))
            for _ in range(errors):
                position = rng.randrange(len(sequence))
                if trial % 3 == 0:
                    sequence = sequence[:position] + sequence[position+1:]
                elif trial % 3 == 1:
                    sequence = sequence[:position] + rng.choice('ACGT') + sequence[position:]
                else:
                    sequence = sequence[:position] + rng.choice('ACGT') + sequence[position+1:]
            hit = _primer_hit(sequence, primer)
            assert hit is not None and hit.edits <= errors


def test_reconstruction_failure_cleans_temporary_fastqs(tmp_path, template, monkeypatch):
    import csv
    from dataclasses import asdict
    from mlvamaps.io import write_fastq
    from mlvamaps.short_reads import run_short_read_call
    panel = tmp_path/'panel.tsv'
    with panel.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(template.locus)), delimiter='\t')
        writer.writeheader()
        writer.writerow(asdict(template.locus))
    first, second = tmp_path/'r1.fq', tmp_path/'r2.fq'
    seq = template.sequence(8)
    write_fastq([pair(seq).read1], first)
    write_fastq([ReadRecord("p/2", revcomp(seq), "I"*len(seq))], second)
    def fail(**kwargs):
        next(kwargs['pairs'])
        raise RuntimeError('injected recovery failure')
    monkeypatch.setattr('mlvamaps.locus_reconstruction.run_reconstructed_fastq_inference', fail)
    out = tmp_path/'out'
    with pytest.raises(RuntimeError, match='injected'):
        run_short_read_call(str(first), str(second), str(panel), str(out), 's', show_progress=False)
    assert not list(out.glob('filtered*.fastq.gz'))


def test_benchmark_counts_uncalled_loci_and_joins_features(tmp_path):
    from scripts.benchmark_cross_mode import compare_runs
    for mode in ('assembly', 'sr'):
        folder = tmp_path/mode
        folder.mkdir()
        (folder/'calls.tsv').write_text('locus_id\trepeat_count\tstatus\nL\t\tNOT_FOUND\n')
    (tmp_path/'sr'/'short_read_repeat_evidence.tsv').write_text(
        'locus_id\tcall_method\tn_repeat_rich\tread_length\tvntr_length\nL\tNO_CALL\t3\t100\t200\n')
    metrics, rows = compare_runs([{'sample_id': 's', 'mode': mode, 'outdir': tmp_path/mode}
                                 for mode in ('assembly', 'sr')])
    assert metrics['assembly:sr']['per_locus_call_rate']['L'] == {'assembly': 0, 'sr': 0}
    assert metrics['assembly:sr']['ambiguous_no_call_rate_by_mode']['sr'] == 1
    assert rows[0]['n_repeat_rich_b'] == '3'
    assert rows[0]['vntr_read_ratio_b'] == 2
    assert rows[0]['comparison'] == 'missing_in_both'


def test_uncertain_fastq_count_is_exported_and_compared_with_assembly(tmp_path, template):
    from dataclasses import asdict
    from mlvamaps.assembly_call import run_assembly_call
    from mlvamaps.io import write_tsv, write_fastq, read_profiles
    from mlvamaps.short_reads import run_short_read_call
    from scripts.benchmark_cross_mode import compare_runs
    # A wide independently measured insert distribution supports a best length
    # but cannot identify one allele with high confidence.
    panel = tmp_path/'panel.tsv'
    write_tsv([asdict(template.locus)], panel, list(asdict(template.locus)))
    for mate, sequence in [(1, template.left), (2, revcomp(template.right))]:
        write_fastq([ReadRecord(str(i), sequence, 'I'*len(sequence)) for i in range(4)], tmp_path/f'r{mate}.fq')
    output = run_short_read_call(str(tmp_path/'r1.fq'), str(tmp_path/'r2.fq'), str(panel),
        str(tmp_path/'sr'), 's', insert_mean=320, insert_sd=100, show_progress=False)
    call = read_profiles(output['calls'])[0]
    assert call['status'] == 'AMBIGUOUS'
    assert float(call['repeat_count']) == 30
    assert float(read_profiles(output['fingerprint'])[0]['L']) == 30
    assert '30 repeats' in output['report'].read_text()
    import csv
    with output['myoga_loci'].open() as handle:
        assert next(csv.DictReader(handle))['repeat_count_min'] == call['repeat_count_min']
    truth = tmp_path/'truth.fa'
    truth.write_text('>truth\n'+template.sequence(30)+'\n')
    run_assembly_call(str(truth), str(panel), str(tmp_path/'assembly'), 's')
    metrics, rows = compare_runs([{'sample_id': 's', 'mode': mode, 'outdir': tmp_path/mode}
                                 for mode in ('assembly', 'sr')])
    comparison = metrics['assembly:sr']
    assert comparison['shared_callable_loci'] == 0
    assert comparison['best_estimates_including_ambiguous']['exact_repeat_concordance'] == 1
    assert rows[0]['best_estimate_absolute_difference'] == 0
    assert rows[0]['status_b'] == 'AMBIGUOUS'


def test_common_interpreter_reports_repeat_snp_without_changing_masked_marker(template):
    sequence = template.sequence(8)
    sequence = sequence[:111]+'C'+sequence[112:]
    genotypes = [genotype_product(LocusProduct('s', 'L', source, 'v', sequence), template.locus)
                 for source in ('assembly', 'long_read', 'short_read')]
    assert genotypes[0] == replace(genotypes[1], product=genotypes[0].product)
    assert len(genotypes[2].repeat_motif_edits) == 1
    assert genotypes[2].repeat_motif_edits[0]['query_position'] == 12
    assert genotypes[2].repeat_motif_edits[0]['alternate'] == 'C'
    inferred = genotype_product(LocusProduct('s', 'L', 'short_read', 'i',
        template.left+'N'*32+template.right), template.locus)
    assert inferred.repeat_motif_edits == ()


def test_direct_flank_indel_uses_observed_primer_product_length(template):
    sequence = template.sequence(8)
    sequence = sequence[:30]+'CGTC'+sequence[30:]
    result = recover_locus(evidence(template, [(sequence, None)]*4), template, 's')
    assert result.products[0].sequence == sequence
    expected = genotype_product(LocusProduct('s', 'L', 'assembly', 'v', sequence), template.locus)
    actual = genotype_product(result.products[0], template.locus)
    assert actual.repeat_count == expected.repeat_count


def test_uninformative_flank_molecule_cannot_bias_candidate_likelihood(template):
    import numpy as np
    from mlvamaps.targeted_reconstruction import candidate_likelihood
    flank = classify_pair(pair(template.left), template)
    assert flank.classes == ('UNINFORMATIVE',)
    without = candidate_likelihood([], template, None, 20, .8)
    with_flank = candidate_likelihood([flank]*100, template, None, 20, .8)
    np.testing.assert_array_equal(without.posterior, with_flank.posterior)
    assert not with_flank.identifiable



def test_censored_boundaries_do_not_favor_short_candidates_at_high_depth(template, monkeypatch):
    import numpy as np
    from mlvamaps.short_read_evidence import MoleculeEvidence
    from mlvamaps.targeted_reconstruction import candidate_likelihood
    def forbidden(*args):
        pytest.fail('a censored boundary does not require per-candidate alignment')
    monkeypatch.setattr('parasail.sw_striped_profile_16', forbidden)
    boundary = MoleculeEvidence('one', 'L', ('LEFT_BOUNDARY',), ('AATG'*25,), (), lower_bound=25)
    one = candidate_likelihood([boundary], template, None, 100, .8)
    many = candidate_likelihood([replace(boundary, molecule_id=str(i)) for i in range(1000)], template, None, 100, .8)
    np.testing.assert_array_equal(one.posterior, many.posterior)
    assert not one.identifiable and one.reason == 'repeat_length_lower_bound'
    assert np.ptp(one.posterior[one.states >= 25]) == 0
    assert all(one.posterior[one.states < 25] == 0)


def test_competitive_fragment_likelihood_is_invariant_to_learned_flank_phase(template):
    import numpy as np
    from mlvamaps.short_read_evidence import MoleculeEvidence
    from mlvamaps.targeted_reconstruction import candidate_likelihood
    insert = InsertDistribution(320, 5, 1, 'override')
    original = MoleculeEvidence('one', 'L', ('FLANK_PAIR',), (), (),
                                lower_bound=10, fragment_offset=200, template=template)
    shifted = replace(original, template=replace(template, left=template.left+'AA'),
                      lower_bound=9.5, fragment_offset=202)
    first = candidate_likelihood([original], template, insert, 100, .8)
    second = candidate_likelihood([shifted], template, insert, 100, .8)
    np.testing.assert_array_equal(first.states, second.states)
    np.testing.assert_allclose(first.posterior, second.posterior)
    assert first.best == 30


def test_candidate_limit_below_observed_boundary_is_never_a_call(template):
    from mlvamaps.short_read_evidence import MoleculeEvidence
    from mlvamaps.targeted_reconstruction import candidate_likelihood
    item = MoleculeEvidence('one', 'L', ('LEFT_BOUNDARY',), (), (), lower_bound=120)
    result = candidate_likelihood([item], template, None, 100, .8)
    assert result.best is None and result.limit_reached and not result.identifiable


def test_single_calibration_fragment_does_not_imply_known_library_variance(template):
    from mlvamaps.short_read_evidence import estimate_insert_distribution, MoleculeEvidence
    from mlvamaps.targeted_reconstruction import candidate_likelihood
    learned = estimate_insert_distribution([320])
    assert learned is not None and learned.pairs_used == 1 and learned.sd >= 32
    pair = MoleculeEvidence('one', 'L', ('FLANK_PAIR',), (), (), fragment_offset=200)
    result = candidate_likelihood([pair], template, learned, 100, .8)
    assert result.best == 30 and not result.identifiable
    assert result.interval[1]-result.interval[0] > 10
    measured = estimate_insert_distribution([], 320, 1)
    assert measured.sd == 1


def test_boundary_uncertainty_remains_open_to_the_candidate_ceiling(template):
    from mlvamaps.short_read_evidence import MoleculeEvidence
    from mlvamaps.targeted_reconstruction import candidate_likelihood
    item = MoleculeEvidence('one', 'L', ('LEFT_BOUNDARY',), (), (), lower_bound=25)
    result = candidate_likelihood([item], template, None, 500, .8)
    assert result.interval == (25, 500)
    assert result.states[-1] == 500 and result.limit_reached
    assert not result.identifiable
