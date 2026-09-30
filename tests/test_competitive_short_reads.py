"""Competitive mapping must retain assembly calibration and novel observations."""
from dataclasses import asdict, replace
import random

import pytest

from mlvamaps.assembly_call import run_assembly_call
from mlvamaps.io import read_profiles, write_fastq, write_tsv
from mlvamaps.models import Locus, ReadRecord
from mlvamaps.sequence import revcomp
from mlvamaps.short_reads import run_short_read_call


@pytest.mark.parametrize('change,mode', [('none', 'single'), ('insertion', 'single'),
    ('deletion', 'single'), ('none', 'paired'), ('none', 'orphan'), ('none', 'cached'),
    ('none', 'primers'), ('insertion', 'primers'), ('deletion', 'primers')])
def test_competitive_single_read_matches_assembly_outside_candidate_grid(tmp_path, change, mode, monkeypatch):
    minimap = '/nonexistent/minimap2'
    rng = random.Random(372)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    left, right, motif = dna(70), dna(70), 'AGTC'
    # Calibration differs from the supplied nonrepeat flank length by half a unit.
    locus = Locus('L', forward_primer=left[:20], reverse_primer=revcomp(right[-20:]),
        left_flank_sequence=left[20:], right_flank_sequence=right[:-20],
        repeat_motif=motif, repeat_unit_length_bp=4, expected_min_repeats=2,
        expected_max_repeats=5, expected_product_size_bp=150, nominal_repeat_units=3)
    sequence = left + motif*9 + right
    if change == 'insertion':
        sequence = sequence[:8] + 'A' + sequence[8:]
    elif change == 'deletion':
        sequence = sequence[:8] + sequence[9:]
    if mode == 'primers':
        locus = replace(locus, left_flank_sequence='', right_flank_sequence='', repeat_motif='NNNN')
    panel, reads, assembly = tmp_path/'panel.tsv', tmp_path/'reads.fq', tmp_path/'assembly.fa'
    write_tsv([asdict(locus)], panel, list(asdict(locus)))
    write_fastq([ReadRecord('one', sequence, 'I'*len(sequence))], reads)
    assembly.write_text('>truth\n'+sequence+'\n')
    mate = None
    database = None
    if mode in {'paired', 'orphan'}:
        mate = tmp_path/'mate.fq'
        write_fastq([ReadRecord('one/1', sequence if mode == 'orphan' else sequence[:140],
                               'I'*(len(sequence) if mode == 'orphan' else 140))], reads)
        second = revcomp(sequence[-140:])
        write_fastq([ReadRecord('one/2', second, ('!' if mode == 'orphan' else 'I')*len(second))], mate)
    if mode == 'cached':
        import shutil
        minimap = shutil.which('minimap2')
        if minimap is None:
            pytest.skip('minimap2 is required for database calling')
        from mlvamaps.candidate_contexts import generate_candidate_contexts, write_candidate_contexts
        database = tmp_path/'database'
        paths = write_candidate_contexts(generate_candidate_contexts([locus]), database/'competitive_mapping')
        (database/'competitive_mapping'/'short.mmi').write_bytes(b'unused reference index')
        monkeypatch.setattr('mlvamaps.minimap_mapping.minimap2_version', lambda _: 'not_used_in_this_check')
        monkeypatch.setattr('mlvamaps.mapping_classification.run_mapping_classification', lambda **kwargs: {})
    reference = run_assembly_call(str(assembly), str(panel), str(tmp_path/'assembly'), 's')
    result = run_short_read_call(str(reads), str(mate) if mate else None, str(panel), str(tmp_path/'sr'), 's',
        database_path=str(database) if database else None,
        minimap2_bin=minimap, threads=1, min_depth=100,
        show_progress=False)
    truth = read_profiles(reference['calls'])[0]
    actual = read_profiles(result['calls'])[0]
    assert float(actual['repeat_count']) == float(truth['repeat_count']) == 9.5
    assert float(actual['product_size_bp']) == float(truth['product_size_bp'])
    assert actual['status'] == 'PASS'
    assert actual['read_depth'] == '1'

    if database:
        assert 'candidate_metadata' not in result
    variants = read_profiles(result['reconstructed_locus_variants'])
    assert variants and float(variants[0]['repeat_count']) == 9.5
    assembly_variants = read_profiles(reference['reconstructed_locus_variants'])
    assert variants[0]['combined_marker'] == assembly_variants[0]['combined_marker']
    assert variants[0]['snp_sequence'] == assembly_variants[0]['snp_sequence']


def test_boundary_reads_do_not_select_short_candidates(tmp_path):
    minimap = '/nonexistent/minimap2'
    from scripts.benchmark_reconstruction import make_inputs
    first, second, panel = make_inputs(tmp_path, 120, 'unresolved')
    result = run_short_read_call(str(first), str(second), str(panel), str(tmp_path/'result'), 's',
        minimap2_bin=minimap, threads=1, show_progress=False)
    calls = read_profiles(result['calls'])
    assert len(calls) == 6
    assert all(row['present'] == 'yes' and row['status'] != 'PASS' for row in calls)
    assert all(row['repeat_count'] != '15' for row in calls)


def test_partial_half_repeat_preserves_snp_mask_boundaries(tmp_path):
    from mlvamaps.locus_products import LocusProduct, genotype_product
    from mlvamaps.models import ReadPair
    from mlvamaps.short_read_evidence import classify_pair, panel_template
    from mlvamaps.targeted_reconstruction import direct_products
    from scripts.benchmark_cross_mode import snp_comparison

    rng = random.Random(418)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    left, right, motif = dna(40), dna(40), dna(12)
    locus = Locus('L', forward_primer=left[:20], reverse_primer=revcomp(right[-20:]),
        left_flank_sequence=left[20:], right_flank_sequence=right[:-20],
        repeat_motif=motif, repeat_unit_length_bp=12)
    sequence = left+motif*5+motif[:6]+right
    # The repeat is fully spanned, but only the first ten bases of the right
    # flank are observed. A motif-only mask would retain the half unit as SNPs.
    partial = sequence[:-30]
    item = classify_pair(ReadPair('one', ReadRecord('one', partial, 'I'*len(partial))), panel_template(locus))
    products = direct_products([item], panel_template(locus), 's', .01, 0)
    actual = genotype_product(products[0], locus)
    expected = genotype_product(LocusProduct('s', 'L', 'assembly', 'truth', sequence), locus)
    assert actual.repeat_count == expected.repeat_count == 5.5
    assert actual.repeat_sequence == expected.repeat_sequence
    assert actual.snp_sequence == expected.snp_sequence[:-30]+'N'*30
    assert snp_comparison(actual.snp_sequence, expected.snp_sequence) == (50, 50)


@pytest.mark.parametrize('mode', ['novel', 'nested_database', 'paired', 'orphan', 'inferred', 'unknown_insert',
                                  'boundary', 'competing_locus'])
def test_reference_calling_measures_reads_missing_both_primers(tmp_path, mode):
    import shutil
    from mlvamaps.unified_fastq import run_unified_fastq_inference
    minimap = shutil.which('minimap2')
    if minimap is None:
        pytest.skip('minimap2 is required for reference rescue integration')
    rng = random.Random(1947)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    left, right, motif = dna(180), dna(180), dna(9)
    locus = Locus('L', forward_primer=left[:20], reverse_primer=revcomp(right[-20:]),
        repeat_motif='N'*9, repeat_unit_length_bp=9,
        expected_product_size_bp=432, nominal_repeat_units=8)
    loci = [locus]
    if mode == 'competing_locus':
        loci.append(replace(locus, locus_id='other'))
    database = tmp_path/'references'
    resource = database/'database' if mode == 'nested_database' else database
    resource.mkdir(parents=True)
    for target in loci:
        (resource/f'{target.locus_id}.fasta').write_text('>reference\n'+left+motif*8+right+'\n')
    # The observed allele is outside the reference's expanded 7..9 grid.
    sequence = left+motif*11+right
    read = left[-60:]+motif*7 if mode == 'boundary' else sequence[100:-100]
    reads = tmp_path/'reads.fq'
    mate = None
    if mode in {'paired', 'inferred', 'unknown_insert'}:
        mate = tmp_path/'mate.fq'
        second = dna(150) if mode == 'paired' else revcomp(right[:80])
        if mode != 'paired':
            read = left[-80:]
        write_fastq([ReadRecord('one/2', second, 'I'*len(second))], mate)
    write_fastq([ReadRecord('one/1' if mate else 'one', read, 'I'*len(read))], reads)
    kwargs = dict(reads1=reads, reads2=mate, loci=loci, sample_id='s', technology='illumina',
        minimap2_bin=minimap, threads=1, minimum_molecules=0, minimum_probability=.8,
        maximum_candidate_repeat_count=10)
    if mode == 'orphan':
        orphan = tmp_path/'orphans.fq'
        reads.rename(orphan)
        reads.write_text('')
        kwargs['orphan_path'] = orphan
    engine = run_unified_fastq_inference
    if mode == 'inferred':
        from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
        engine = run_reconstructed_fastq_inference
        kwargs.update(insert_mean=160+11*9, insert_sd=.1, maximum_candidate_repeat_count=100)
    without, _, _, _ = engine(**kwargs, database_path=None, outdir=tmp_path/'free')
    assert all(row['repeat_count'] == '' for row in without)
    calls, _, _, paths = engine(**kwargs, database_path=database,
        outdir=tmp_path/'rescued', keep_alignments=True)
    if mode in {'novel', 'nested_database', 'paired', 'orphan', 'inferred'}:
        assert calls[0]['repeat_count'] == 11
        assert calls[0]['product_size_bp'] == len(sequence)
        assert calls[0]['status'] == 'called'
        assert 'reference_assisted:reference' in calls[0]['reason']
        variants = read_profiles(paths['reconstructed_locus_variants'])
        expected = ('N'*100+left[-80:]+'N'*99+right[:80]+'N'*100
                    if mode == 'inferred' else 'N'*100+read+'N'*100)
        assert variants[0]['sequence'] == expected
        summary = read_profiles(paths['reference_calling'])[0]
        assert summary['status'] == summary['call_status'] == 'called'
        if mode == 'inferred':
            assert float(summary['insert_mean']) == 259
            assert float(summary['insert_sd']) == .1
    else:
        assert all(row['repeat_count'] == '' for row in calls)
        assert all(row['status'] != 'called' for row in calls)
    assert paths['reference_calling_alignments_0'].is_file()


def test_tied_fragment_lengths_remain_uncalled_at_low_probability_threshold():
    from mlvamaps.competitive_short_reads import candidate_likelihood
    from mlvamaps.short_read_evidence import InsertDistribution, MoleculeEvidence, RepeatTemplate
    locus = Locus('L', repeat_unit_length_bp=4)
    template = RepeatTemplate(locus, 'A'*40, 'C'*40, 'AGTC', 4)
    item = MoleculeEvidence('one', 'L', ('FLANK_PAIR',), (), (), fragment_offset=80)
    result = candidate_likelihood([item], template, InsertDistribution(121, .1, 1, 'override'), 100, .3)
    assert not result.identifiable and result.best is None
    assert result.reason == 'minimum_likelihood_tie'
