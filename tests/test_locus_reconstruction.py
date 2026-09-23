from dataclasses import replace

import numpy as np
import pytest

from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.sequence import revcomp
from mlvamaps.locus_products import LocusProduct, genotype_product
from mlvamaps.locus_reconstruction import recover_long_reads, reconstruct_short_products
from mlvamaps.repeat_likelihood import (
    panel_template, classify_pair, infer_repeat, InsertDistribution, estimate_insert_distribution,
)


@pytest.fixture
def template():
    return panel_template(Locus('X', forward_primer='ACGTCAGTACGATCG', reverse_primer='TCGTAGCACTGATCG',
        left_flank_sequence='GCTAGCTACGTTACAGTGCA', right_flank_sequence='CGACTGATGCATCGTACGAC',
        repeat_motif='AATG', repeat_unit_length_bp=4, expected_min_repeats=3, expected_max_repeats=10))


def pair(a, b=None, name='m'):
    return ReadPair(name, ReadRecord(name+'/1', a, 'I'*len(a)), ReadRecord(name+'/2', b, 'I'*len(b)) if b else None)


def test_evidence_classes(template):
    t = template
    seq = t.sequence(12)
    assert 'E' in classify_pair(pair(seq), t).classes
    assert 'F' in classify_pair(pair(t.left + 'AATG'*5), t).classes
    assert 'F' in classify_pair(pair('AATG'*5 + t.right), t).classes
    assert 'S' in classify_pair(pair(t.left, revcomp(t.right)), t).classes
    assert 'FRR' in classify_pair(pair(t.left, revcomp('AATG'*8)), t).classes
    assert classify_pair(pair('AATG'*8), t) is None
    assert classify_pair(pair('CCCCCCCCCCCCCCCCCCCC'), t) is None
    assert 'discordant' in classify_pair(pair(t.left, t.right), t).classes


@pytest.mark.parametrize('repeat', [2, 8, 12, 12.5, 40])
def test_novel_and_outside_range(template, repeat):
    items = [classify_pair(pair(template.sequence(repeat), name=str(i)), template) for i in range(4)]
    result = infer_repeat(items, template)
    assert result.best == repeat
    assert result.identifiable
    products = reconstruct_short_products(items, template, result, 'sample')
    assert genotype_product(products[0], template.locus).repeat_count == repeat


def test_spanning_pair_length(template):
    t = template
    items = [classify_pair(pair(t.left, revcomp(t.right), str(i)), t) for i in range(40)]
    insert = InsertDistribution(len(t.sequence(30)), 4, 40, 'override')
    result = infer_repeat(items, t, insert)
    assert result.best == 30
    assert result.identifiable
    assert not infer_repeat(items, t).identifiable


def test_lower_bounds_never_force_exact_call(template):
    items = [classify_pair(pair(template.left + 'AATG'*25), template)] * 20
    result = infer_repeat(items, template, maximum=30)
    assert not result.identifiable
    assert result.interval[0] >= 24
    assert result.limit_reached


def test_insert_robustness_and_override():
    result = estimate_insert_distribution([299, 300, 301]*20 + [9000, 10000])
    assert result.mean == 300
    assert result.pairs_used == 60
    assert estimate_insert_distribution([], 500, 20).mean == 500
    assert estimate_insert_distribution([300]*2) is None
    with pytest.raises(ValueError):
        estimate_insert_distribution([], 300, None)


@pytest.mark.parametrize('fraction', [.9, .7, .5])
@pytest.mark.parametrize('change_repeat', [True, False])
def test_mixtures_and_snp_haplotypes(template, fraction, change_repeat):
    first = template.sequence(8)
    second = template.sequence(12 if change_repeat else 8)
    second = second[:20] + ('A' if second[20] != 'A' else 'C') + second[21:]
    sequences = [first]*round(100*fraction) + [second]*round(100*(1-fraction))
    items = [classify_pair(pair(seq, name=str(i)), template) for i,seq in enumerate(sequences)]
    inference = infer_repeat(items, template)
    products = reconstruct_short_products(items, template, inference, 'sample')
    assert len(products) == 2
    assert sorted(p.estimated_fraction for p in products) == pytest.approx(sorted([fraction, 1-fraction]), abs=.01)
    assert {p.sequence for p in products} == {first, second}
    assert len({genotype_product(p, template.locus).haplotype_id for p in products}) == 2


def test_cross_modality_product_equivalence(template):
    seq = template.sequence(9)
    seq = seq[:20] + 'A' + seq[21:]
    assembly = genotype_product(LocusProduct('s', 'X', 'assembly', 'a', seq), template.locus)
    pairs = [pair(seq if i % 2 else revcomp(seq), name=str(i)) for i in range(6)]
    lr, audit, _ = recover_long_reads(pairs, [template.locus], 's', max_anchor_edits=1)
    items = [classify_pair(p, template) for p in pairs]
    sr = reconstruct_short_products(items, template, infer_repeat(items, template), 's')
    for product in lr + sr:
        observed = genotype_product(product, template.locus)
        assert observed.repeat_count == assembly.repeat_count == 9
        assert observed.snp_sequence == assembly.snp_sequence
        assert observed.combined_marker == assembly.combined_marker
    assert {r['orientation'] for r in audit} == {'+', '-'}


def test_spanning_mixture(template):
    # Two fragment geometries imply repeat states 20 and 40 under one library.
    from mlvamaps.repeat_likelihood import MoleculeEvidence
    items = [MoleculeEvidence(str(i), 'X', ('S',), (), (), fragment_offset=500-r*4)
             for i,r in enumerate([20]*70+[40]*30)]
    result = infer_repeat(items, template, InsertDistribution(500, 4, 100, 'override'))
    assert result.fractions[20] == pytest.approx(.7, abs=.01)
    assert result.fractions[40] == pytest.approx(.3, abs=.01)
    assert result.fractions.get(26, 0) < .01


def write_panel(tmp_path, template):
    from dataclasses import asdict
    import csv
    panel = tmp_path / 'panel.tsv'
    row = asdict(template.locus)
    with panel.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(row), delimiter='\t')
        writer.writeheader()
        writer.writerow(row)
    return panel


def test_end_to_end_three_modes(tmp_path, template):
    from mlvamaps.assembly_call import run_assembly_call
    from mlvamaps.pipeline import run_call
    from mlvamaps.short_reads import run_short_read_call
    from mlvamaps.io import write_fastq, read_profiles
    seq = template.sequence(9)
    seq = seq[:20] + 'A' + seq[21:]
    panel = write_panel(tmp_path, template)
    fasta = tmp_path / 'truth.fa'
    fasta.write_text('>truth\n'+seq+'\n')
    reads = tmp_path / 'reads.fq'
    write_fastq([ReadRecord(str(i), seq, 'I'*len(seq)) for i in range(8)], reads)
    mate = tmp_path / 'mate.fq'
    write_fastq([ReadRecord(str(i), revcomp(seq), 'I'*len(seq)) for i in range(8)], mate)
    outputs = [
        run_assembly_call(str(fasta), str(panel), str(tmp_path/'assembly'), 's'),
        run_call(str(reads), str(panel), str(tmp_path/'lr'), 's'),
        run_short_read_call(str(reads), str(mate), str(panel), str(tmp_path/'sr'), 's', show_progress=False),
    ]
    calls = [read_profiles(o['calls'])[0] for o in outputs]
    assert {float(row['repeat_count']) for row in calls} == {9}
    products = [read_profiles(o['reconstructed_locus_variants'])[0] for o in outputs]
    assert len({p['combined_marker'] for p in products}) == 1
    assert len({p['snp_sequence'] for p in products}) == 1
    assert 'repeat_likelihoods' in outputs[2]


def test_database_does_not_change_rich_panel_novel_call(tmp_path, template):
    from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
    from mlvamaps.io import write_fastq
    database = tmp_path/'db'
    database.mkdir()
    (database/'X.fasta').write_text('>old\n'+template.sequence(8)+'\n')
    reads = tmp_path/'reads.fq'
    seq = template.sequence(17)
    write_fastq([ReadRecord(str(i), seq, 'I'*len(seq)) for i in range(6)], reads)
    calls = []
    for name, db in [('without', None), ('with', database)]:
        rows, _, _, _ = run_reconstructed_fastq_inference(reads1=reads, reads2=None, loci=[template.locus],
            database_path=db, outdir=tmp_path/name, sample_id='s', technology='illumina',
            minimum_molecules=3, minimum_probability=.8)
        calls.append(rows[0]['repeat_count'])
    assert calls == [17, 17]


def test_partial_snp_consensus_and_uncovered_bases(template):
    seq = template.sequence(30)
    items = [classify_pair(pair(seq[:55], revcomp(seq[-55:]), str(i)), template) for i in range(40)]
    result = infer_repeat(items, template, InsertDistribution(len(seq), 4, 40, 'override'))
    products = reconstruct_short_products(items, template, result, 's')
    assert len(products) == 1
    assert products[0].sequence == seq
    assert genotype_product(products[0], template.locus).repeat_count == 30


def test_partial_paired_snp_mixture(template):
    original = template.sequence(30)
    variant = original[:20] + 'A' + original[21:]
    truth = [original]*70 + [variant]*30
    items = [classify_pair(pair(seq[:55], revcomp(seq[-55:]), str(i)), template) for i,seq in enumerate(truth)]
    inference = infer_repeat(items, template, InsertDistribution(len(original), 4, 100, 'override'))
    products = reconstruct_short_products(items, template, inference, 's')
    assert {p.sequence for p in products} == {original, variant}
    assert sorted(p.estimated_fraction for p in products) == pytest.approx([.3, .7])


def test_low_quality_snp_is_not_reported_as_supported(template):
    seq = template.sequence(30)
    items = []
    for i in range(40):
        p = pair(seq[:55], revcomp(seq[-55:]), str(i))
        p = replace(p, read1=replace(p.read1, quality='I'*20+'!'+ 'I'*34))
        items.append(classify_pair(p, template))
    inference = infer_repeat(items, template, InsertDistribution(len(seq), 4, 40, 'override'))
    products = reconstruct_short_products(items, template, inference, 's')
    assert products[0].sequence[20] == 'N'


def test_shared_flanks_are_not_double_counted(tmp_path, template):
    from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
    from mlvamaps.io import write_fastq
    reads = tmp_path/'reads.fq'
    seq = template.sequence(8)
    write_fastq([ReadRecord('r', seq, 'I'*len(seq))], reads)
    rows, audit, _, _ = run_reconstructed_fastq_inference(reads1=reads, reads2=None,
        loci=[template.locus, replace(template.locus, locus_id='Y')], database_path=None,
        outdir=tmp_path/'out', sample_id='s', technology='illumina', minimum_molecules=1, minimum_probability=.8)
    assert all(row['repeat_count'] == '' for row in rows)
    assert all(row['classes'] == 'ambiguous_locus' for row in audit)


def test_cli_engine_and_insert_options():
    from mlvamaps.cli import build_parser
    args = build_parser().parse_args(['call', '-i', 'reads.fq', '-p', 'primers.tsv',
                                     '--sr-engine', 'competitive', '--lr-engine', 'competitive',
                                     '--insert-mean', '400', '--insert-sd', '30'])
    assert args.insert_mean == 400 and args.insert_sd == 30
    assert args.sr_engine == args.lr_engine == 'competitive'


def test_benchmark_missing_calls_and_tied_ranks(tmp_path):
    from scripts.benchmark_cross_mode import compare_runs, ranks
    for mode, count in [('assembly', '9'), ('sr', '10'), ('lr', '')]:
        directory = tmp_path/mode
        directory.mkdir()
        (directory/'calls.tsv').write_text('locus_id\trepeat_count\tstatus\nX\t'+count+'\t'+('PASS' if count else 'NOT_FOUND')+'\n')
    result, rows = compare_runs([{'sample_id':'s', 'mode':mode, 'outdir':tmp_path/mode} for mode in ['assembly', 'sr', 'lr']])
    assert result['assembly:sr']['within_one_repeat_concordance'] == 1
    assert result['assembly:sr']['mean_absolute_repeat_difference'] == 1
    assert result['assembly:lr']['exact_repeat_concordance'] is None
    assert ranks([1, 1, 3]).tolist() == [.5, .5, 2]


def test_mask_excludes_repeat_length_differences(template):
    products = [genotype_product(LocusProduct('s', 'X', 'short_read', str(r), template.sequence(r)), template.locus)
                for r in [3, 8, 12.5, 40]]
    assert len({p.snp_sequence for p in products}) == 1
    assert len({p.haplotype_id for p in products}) == 1


def test_singleton_secondary_sequence_remains_trace(template):
    seq = template.sequence(8)
    alt = seq[:20] + 'A' + seq[21:]
    items = [classify_pair(pair(s, name=str(i)), template) for i,s in enumerate([seq]*20+[alt])]
    products = reconstruct_short_products(items, template, infer_repeat(items, template), 's')
    assert {p.sequence for p in products} == {seq, alt}
    assert next(p for p in products if p.sequence == alt).evidence['meaningful'] == 'no'


def test_flank_insert_geometry(template):
    from mlvamaps.repeat_likelihood import flank_insert_length
    import random
    rng = random.Random(12)
    flank = ''.join(rng.choices('ACGT', k=600))
    t = replace(template, left=flank)
    assert flank_insert_length(pair(flank[20:120], revcomp(flank[220:320])), t) == 300
    assert flank_insert_length(pair(t.left[-100:], revcomp(t.right)), t) is None


def test_overlapping_partial_reads_retain_nonrepeat_indel(template):
    # Product-size calibration is independently applied after sequence recovery.
    truth = template.sequence(30)
    truth = truth[:20] + 'G' + truth[20:]
    items = [classify_pair(pair(truth[:55], revcomp(truth[-55:]), str(i)), template) for i in range(40)]
    inference = infer_repeat(items, template, InsertDistribution(len(truth), 4, 40, 'override'))
    products = reconstruct_short_products(items, template, inference, 's')
    assert products[0].sequence == truth
    assert genotype_product(products[0], template.locus).repeat_count == genotype_product(
        LocusProduct('s', 'X', 'assembly', 'a', truth), template.locus).repeat_count


def test_shared_snp_alignment_ignores_repeat_and_retains_snp(tmp_path, template):
    from mlvamaps.locus_products import write_products
    from mlvamaps.io import read_profiles
    primary = template.sequence(8)
    secondary = template.sequence(12)
    products = [LocusProduct('s', 'X', 'short_read', 'a', primary, effective_depth=70, estimated_fraction=.7),
                LocusProduct('s', 'X', 'short_read', 'b', secondary, effective_depth=30, estimated_fraction=.3)]
    paths = write_products(products, [template.locus], tmp_path)
    assert read_profiles(paths['mapping_snps']) == []
    products[1] = replace(products[1], sequence=secondary[:20]+'A'+secondary[21:])
    paths = write_products(products, [template.locus], tmp_path)
    rows = read_profiles(paths['mapping_snps'])
    assert len(rows) == 1
    assert float(rows[0]['alternate_frequency']) == .3
    assert rows[0]['coordinate_system'] == 'repeat_masked_product'


def test_custom_rounding_is_retained_in_canonical_output(tmp_path, template):
    from mlvamaps.locus_products import write_products
    from mlvamaps.io import read_profiles
    seq = template.sequence(8) + 'A'
    product = LocusProduct('s', 'X', 'assembly', 'v', seq, round_tolerance=.4)
    assert genotype_product(product, template.locus).repeat_count == 8
    paths = write_products([product], [template.locus], tmp_path)
    assert float(read_profiles(paths['reconstructed_locus_variants'])[0]['repeat_count']) == 8


def test_candidate_memory_safety_guard(template):
    t = replace(template, locus=replace(template.locus, expected_max_repeats=20000))
    with pytest.raises(ValueError, match='memory safety limit'):
        infer_repeat([], t, maximum=20000)
