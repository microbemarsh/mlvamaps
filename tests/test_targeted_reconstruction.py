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


def test_removed_engine_flag_is_rejected():
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


def test_conflicting_singleton_lengths_withhold_precise_call(template):
    result = recover_locus(evidence(template, [(template.sequence(8), None),
                                               (template.sequence(12), None)]), template, 's')
    assert result.method == 'AMBIGUOUS'
    assert not result.products
    assert result.confidence < .8


def test_database_natural_motif_and_flanks_are_retained(tmp_path, template):
    from mlvamaps.locus_reconstruction import recruitment_templates
    natural = template.sequence(8)
    natural = natural[:110]+'C'+natural[111:]
    (tmp_path/'L.fasta').write_text('>known\n'+natural+'\n')
    contexts = recruitment_templates([template.locus], tmp_path, {'L': template})
    assert any(t.sequence(8) == natural for t in contexts.values())
    assert len(contexts) == 2


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
