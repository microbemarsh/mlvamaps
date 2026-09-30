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
