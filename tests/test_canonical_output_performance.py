"""Scientific and output equivalence for native canonical genotyping."""
import csv
from dataclasses import replace
import io
import random

import parasail
import pytest

from mlvamaps.combined_marker_phylogeny import (
    _IUPAC, _iupac_find, _longest_motif_run, _longest_motif_run_scalar,
    decompose_marker_sequence,
)
from mlvamaps.locus_products import LocusProduct, genotype_product, write_product_alignments
from mlvamaps.models import Locus


def scalar_find(pattern, sequence, start=0):
    pattern, sequence = pattern.upper(), sequence.upper()
    if not pattern:
        return -1
    for index in range(start, len(sequence) - len(pattern) + 1):
        if all(base in _IUPAC.get(code, frozenset(code)) for code, base in zip(pattern, sequence[index:])):
            return index
    return -1


def scalar_run(sequence, motif):
    motif = motif.upper()
    return _longest_motif_run_scalar(sequence, motif) if motif and set(motif) != {'N'} else None


def test_native_iupac_search_preserves_literal_sequence_ambiguity():
    rng = random.Random(362)
    cases = [('RY', 'AGCTNRY'), ('N', 'NNNACGT'), ('[]', 'AC[]GT'), ('é', 'ACÉGT'), ('', 'ACGT')]
    cases += [(''.join(rng.choices('ACGTNRYSW', k=rng.randrange(1, 12))),
               ''.join(rng.choices('ACGTNRYSW-', k=60))) for _ in range(200)]
    for pattern, sequence in cases:
        for start in (-1, 0, 5, len(sequence), len(sequence) + 2):
            assert _iupac_find(pattern, sequence, start) == scalar_find(pattern, sequence, start)


def test_native_motif_runs_match_scalar_with_errors_phases_and_ties():
    rng = random.Random(919)
    cases = [('', 'A'), ('A', ''), ('AAAA', 'N'), ('AAAA', 'AAAAA'),
             ('AAAACTAAAA', 'AA'), ('ACGTACGT', 'RYSW'), ('NACNNACN', 'NAC'),
             ('éééAééé', 'é'), ('ACGTNNNNACGT', 'ACGT'), ('aaaAAAA', 'a')]
    for unit in (1, 4, 8, 12, 60, 128):
        motif = ''.join(rng.choices('ACGT', k=unit))
        for phase in (0, unit // 2):
            sequence = (motif * 7)[phase:]
            for errors in (0, unit // 8, unit // 8 + 1):
                bases = list(sequence)
                for index in rng.sample(range(len(bases)), errors):
                    bases[index] = next(b for b in 'ACGT' if b != bases[index])
                cases.append((''.join(bases), motif))
        cases.append((''.join(rng.choices('ACGTNR-', k=500)), motif))
    # Cross the native block boundary without expensive scalar repeat expansion.
    cases.append((''.join(rng.choices('ACGT', k=18000)), 'ACGT'*16))
    for sequence, motif in cases:
        assert _longest_motif_run(sequence, motif) == scalar_run(sequence, motif)


def test_marker_components_match_scalar_in_both_orientations(monkeypatch):
    import mlvamaps.combined_marker_phylogeny as marker
    from mlvamaps.sequence import revcomp
    locus = Locus('L', forward_primer='ACGTCGACTG', reverse_primer='TAGCGTACGA',
                  repeat_motif='AATG', repeat_unit_length_bp=4)
    native_find, native_run = marker._iupac_find, marker._longest_motif_run
    for with_flanks in (False, True):
        template = replace(locus, left_flank_sequence='GCACTGA' if with_flanks else '',
                           right_flank_sequence='CATGACG' if with_flanks else '')
        sequence = locus.forward_primer + 'GCACTGA' + 'AATG'*40 + 'CATGACG' + revcomp(locus.reverse_primer)
        for seq in (sequence, revcomp(sequence), sequence.replace('GCACTGA', 'NCACAGA')):
            monkeypatch.setattr(marker, '_iupac_find', scalar_find)
            monkeypatch.setattr(marker, '_longest_motif_run', scalar_run)
            expected = decompose_marker_sequence(template, seq)
            monkeypatch.setattr(marker, '_iupac_find', native_find)
            monkeypatch.setattr(marker, '_longest_motif_run', native_run)
            assert decompose_marker_sequence(template, seq) == expected


def test_canonical_identity_alignment_matches_parasail(tmp_path, monkeypatch):
    locus = Locus('L')
    sequences = ['ACGTNNACGTTGCA', 'ACGTNNACGTTGCA', 'ACGTNAACGTTGCA', 'ACGTNACGTTGCA', 'ACGTRNACGTTGCA']
    genotypes = [genotype_product(LocusProduct('s', 'L', 'short_read', f'v{i}', seq,
                 estimated_fraction=1/(i+1)), locus) for i, seq in enumerate(sequences)]
    native = parasail.nw_trace_striped_32
    called = []
    def tracked(query, reference, *args):
        called.append(query)
        return native(query, reference, *args)
    monkeypatch.setattr(parasail, 'nw_trace_striped_32', tracked)
    paths = write_product_alignments(genotypes, tmp_path)
    with paths['canonical_alignments'].open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert called == sequences[2:]
    matrix = parasail.matrix_create('ACGTN', 2, -4)
    for row, sequence in zip(rows, sequences):
        expected = native(sequence, sequences[0], 5, 1, matrix)
        assert row['aligned_query'] == expected.traceback.query
        assert row['aligned_reference'] == expected.traceback.ref
        assert row['cigar'] == expected.cigar.decode.decode()
        assert int(row['alignment_score']) == expected.score


def test_streamed_membership_and_prediction_files_match_dict_writer(tmp_path, monkeypatch):
    from mlvamaps.locus_reconstruction import write_compatibility_products, _common_call
    from mlvamaps.pipeline import ASV_MEMBERSHIP_FIELDS, PREDICTION_FIELDS
    import mlvamaps.locus_products as products_module
    locus = Locus('L', repeat_motif='AATG', repeat_unit_length_bp=4,
                  expected_product_size_bp=44, nominal_repeat_units=10)
    product = LocusProduct('s', 'L', 'short_read', 'v', 'AATG'*12, support_count=9003,
                           effective_depth=9003, evidence={'molecule_ids':
                           [f'm{i}' for i in range(9000)] + ['tab\tid', 'new\nline', 'quote"id']})
    expected_genotype = genotype_product(product, locus)
    def unused(*args):
        pytest.fail('repeat-only output must not compute SNP masking')
    monkeypatch.setattr(products_module, 'decompose_marker_sequence', unused)
    paths = write_compatibility_products([product], [locus], tmp_path)
    call = _common_call(locus, [product], product.support_count, 's', 'illumina', 1, .8)
    assert call['repeat_count'] == expected_genotype.repeat_count
    for key, fields in [('mapped_read_memberships', ASV_MEMBERSHIP_FIELDS), ('read_predictions', PREDICTION_FIELDS)]:
        expected = io.StringIO(newline='')
        writer = csv.DictWriter(expected, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for molecule in product.evidence['molecule_ids']:
            writer.writerow({'sample_id': 's', 'read_id': molecule, 'locus_id': 'L', 'variant_id': 'v',
                             'repeat_count': expected_genotype.repeat_count,
                             'predicted_repeat_count': expected_genotype.repeat_count,
                             'raw_repeat_count_estimate': expected_genotype.repeat_count_raw,
                             'probability': product.reconstruction_confidence, 'evidence_weight': 1})
        assert paths[key].read_bytes() == expected.getvalue().encode()
