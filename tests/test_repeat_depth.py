"""Coverage estimates are reportable, but never acquire unsupported confidence."""
from dataclasses import replace
import csv
import json
import random

import numpy as np
import pytest

from mlvamaps.locus_reconstruction import locus_templates, run_reconstructed_fastq_inference
from mlvamaps.models import Locus, ReadPair, ReadRecord
from mlvamaps.repeat_depth import estimate_graph_lengths
from mlvamaps.sequence import revcomp
from mlvamaps.short_read_evidence import MoleculeEvidence
from mlvamaps.targeted_reconstruction import LocusRecovery


def reads(count, depth=1, noise=False, unit=9, seed=189):
    rng = random.Random(seed)
    def dna(n):
        return ''.join(rng.choices('ACGT', k=n))
    left, right, motif = dna(90), dna(90), dna(unit)
    sequence = left+motif*count+right
    locus = Locus('L', forward_primer=left[:20], reverse_primer=revcomp(right[-20:]),
        repeat_unit_length_bp=unit, expected_product_size_bp=180+8*unit, nominal_repeat_units=8)
    world = dna(180)+sequence+dna(180)
    pairs, items = [], []
    for copy in range(depth):
        for start in range(len(world)-149):
            seq = world[start:start+150]
            if noise and copy % 3:
                pos = rng.randrange(150)
                seq = seq[:pos]+rng.choice([b for b in 'ACGT' if b != seq[pos]])+seq[pos+1:]
            name = f'{copy}:{start}'
            # Alternate input strands; retained evidence is locus oriented.
            pairs.append(ReadPair(name, ReadRecord(name, revcomp(seq) if start % 2 else seq, 'I'*len(seq))))
            if start < 180+90 or start+150 > 180+90+count*unit:
                items.append(MoleculeEvidence(name, 'L', ('UNINFORMATIVE',), (seq,), ('+',)))
    return locus, sequence, pairs, items


@pytest.mark.parametrize('count,depth,noise,unit', [(30, 1, False, 9), (70, 1, False, 9),
    (70, 5, True, 9), (30, 2, True, 39), (140, 1, False, 9)])
def test_cycle_multiplicity_uses_all_reads_not_anchor_selected_depth(count, depth, noise, unit):
    locus, sequence, pairs, items = reads(count, depth, noise, unit)
    recovery = LocusRecovery(method='AMBIGUOUS', reason='assembly_node_limit', limit_reached=True)
    diagnostics = estimate_graph_lengths([recovery], [locus], locus_templates([locus]), {'L': items}, lambda: iter(pairs))
    estimate = diagnostics['L']
    assert abs(estimate['amplicon_length']-len(sequence)) <= 9
    assert abs(estimate['repeat_count']-count) <= 1
    assert estimate['repeat_interval'][0] <= count <= estimate['repeat_interval'][1]
    assert estimate['status'] == 'diagnostic_only'
    assert recovery.best is None and recovery.product_size_bp is None
    assert recovery.method == 'AMBIGUOUS' and not recovery.identifiable
    assert recovery.limit_reached and not recovery.products
    assert diagnostics['L']['arm_depths'][0] > 0
    assert 'assembly_node_limit' in recovery.reason
    assert 'depth_only_not_length_identifying' in recovery.reason


def test_graph_ignores_panel_nominal_count_and_unconnected_background():
    locus, sequence, pairs, items = reads(60)
    # Change the nominal allele while preserving only the length calibration.
    locus = replace(locus, expected_product_size_bp=1080, nominal_repeat_units=100)
    rng = random.Random(711)
    for i in range(2100):
        seq = ''.join(rng.choices('ACGT', k=150))
        items.append(MoleculeEvidence(f'noise{i}', 'L', ('UNINFORMATIVE',), (seq,), ('+',)))
    recovery = LocusRecovery()
    diagnostics = estimate_graph_lengths([recovery], [locus], locus_templates([locus]), {'L': items}, lambda: iter(pairs))
    assert diagnostics['L']['repeat_count'] == 60
    assert diagnostics['L']['amplicon_length'] == len(sequence)
    assert recovery.best is None and recovery.product_size_bp is None


def test_matching_background_inflates_depth_without_establishing_a_locus_count():
    from mlvamaps.locus_reconstruction import _common_call
    from mlvamaps.unified_fastq import common_calls_to_compatibility
    locus, sequence, pairs, items = reads(70)
    templates = locus_templates([locus])
    clean = estimate_graph_lengths([LocusRecovery()], [locus], templates,
                                  {'L': items}, lambda: iter(pairs))['L']
    # Identical repeat sequence can occur outside this locus. It changes neither
    # the target allele nor its flank coverage, but contributes to the ratio.
    repeated = sequence[90:99]*17
    background = [ReadPair(f'background{i}', ReadRecord(f'background{i}', repeated, 'I'*len(repeated)))
                  for i in range(5000)]
    recovery = LocusRecovery()
    inflated = estimate_graph_lengths([recovery], [locus], templates,
        {'L': items}, lambda: iter(pairs+background))['L']
    assert inflated['repeat_count'] > 5*clean['repeat_count']
    assert inflated['arm_depths'] == clean['arm_depths']
    assert inflated['shared_edges'] == 0  # A panel-only sharing check misses it.
    assert inflated['status'] == 'diagnostic_only'
    assert recovery.best is None and recovery.interval is None and recovery.product_size_bp is None
    call = _common_call(locus, [], len(items), 's', 'illumina', 0, .8, recovery, inflated)
    assert call['status'] == 'estimated'
    assert call['repeat_count'] == inflated['repeat_count']
    assert call['confidence'] == 0 and call['confidence_kind'] == 'unvalidated_depth'
    exported = common_calls_to_compatibility([call])[0]
    assert exported['status'] == 'ESTIMATED' and exported['allele_confidence'] == 0


def test_disconnected_signal_does_not_fabricate_a_graph_length():
    locus, _, pairs, items = reads(60)
    items = [item for item in items if locus.forward_primer in item.sequences[0]]
    recovery = LocusRecovery()
    def forbidden_replay():
        pytest.fail('No connected graph should not trigger an input rescan')
    diagnostics = estimate_graph_lengths([recovery], [locus], locus_templates([locus]), {'L': items}, forbidden_replay)
    assert recovery.best is None and recovery.product_size_bp is None
    assert diagnostics['L']['reason'] == 'no_primer_connected_kmer_graph'


def test_multiple_loci_share_one_quality_filtered_input_scan():
    loci, recoveries, evidence, pairs = [], [], {}, []
    for i, count in enumerate((30, 70)):
        locus, _, sample_pairs, items = reads(count, seed=189+i)
        locus = replace(locus, locus_id=f'L{i}')
        loci.append(locus)
        recoveries.append(LocusRecovery())
        evidence[locus.locus_id] = [replace(item, locus_id=locus.locus_id) for item in items]
        pairs.extend(sample_pairs)
        pairs.extend(replace(pair, read1=replace(pair.read1, quality='!'*150)) for pair in sample_pairs)
    scans = []
    def replay():
        scans.append(True)
        return iter(pairs)
    diagnostics = estimate_graph_lengths(recoveries, loci, locus_templates(loci), evidence, replay)
    assert len(scans) == 1
    assert [diagnostics[locus.locus_id]['repeat_count'] for locus in loci] == [30, 70]
    assert all(result.best is None for result in recoveries)


def test_repeat_phase_in_learned_flank_does_not_create_false_one_copy_span():
    from mlvamaps.short_read_evidence import RepeatTemplate, classify_pair
    # A gapped flank alignment used to match the periodic tail alone and
    # report a 152 bp product, although these reads come from a 305 bp locus.
    left = 'TATGCTACTAACTGCGAGGGGATCTCGTAGCCGCGTTAATAGTCTTTGTCTTGCAACCTTCCGGCTGTGAGAGTCTTT'
    right = 'AACTAACTCGGGCGTTCCCCGGAGTCCAGCATGATTGCCGATATATTATACAACATGACTGGTCTACGTTT'
    motif = 'AACTAACTCGGGCGTTACCTTCCGGCTGTGAGAGTCTTT'
    locus = Locus('L', forward_primer=left[:20], reverse_primer=revcomp(right[-20:]), repeat_unit_length_bp=39)
    template = RepeatTemplate(locus, left, right, motif, 39)
    first = (left[22:]+motif*3)[:150]
    mate = motif[-12:]+right
    pair = ReadPair('p', ReadRecord('p/1', first), ReadRecord('p/2', revcomp(mate)))
    item = classify_pair(pair, template)
    assert item is not None and 'FULL_SPAN' not in item.classes
    assert item.observed_repeat is None
    assert item.alignment['0']['right'] is None


def test_depth_only_counts_reach_exports_and_report_with_zero_confidence(tmp_path, monkeypatch):
    from mlvamaps import targeted_reconstruction
    from mlvamaps.unified_fastq import common_calls_to_compatibility
    from mlvamaps.report import _repeat_count_svg
    locus, sequence, pairs, _ = reads(70)
    monkeypatch.setattr(targeted_reconstruction, 'microassemble',
        lambda *args, **kwargs: ('', [], 'assembly_node_limit'))
    calls, _, _, paths = run_reconstructed_fastq_inference(reads1=None, reads2=None,
        pairs=iter(pairs), loci=[locus], database_path=None, outdir=tmp_path,
        sample_id='s', technology='illumina', minimum_molecules=2, minimum_probability=.8)
    call = calls[0]
    assert call['repeat_count'] == 70 and call['product_size_bp'] == len(sequence)
    assert call['status'] == 'estimated'
    assert call['confidence'] == call['best_probability'] == 0
    assert call['confidence_kind'] == 'unvalidated_depth'
    assert call['repeat_count_interval_kind'] == 'coverage_sensitivity'
    assert 'depth_only_not_length_identifying' in call['reason']
    with paths['short_read_repeat_evidence'].open() as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    assert row['best_repeat'] == '70' and float(row['amplicon_length']) == len(sequence)
    assert float(row['confidence']) == 0
    assert row['limit_reached'] == 'True'
    with paths['common_locus_calls'].open() as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    assert row['repeat_count'] == '70' and float(row['product_size_bp']) == len(sequence)
    diagnostic = json.loads(paths['repeat_length_estimates'].read_text())['L']
    assert diagnostic['amplicon_length'] == len(sequence) and diagnostic['repeat_count'] == 70
    assert diagnostic['status'] == 'diagnostic_only'
    exported = common_calls_to_compatibility(calls)
    assert exported[0]['repeat_count'] == 70 and exported[0]['allele_confidence'] == 0
    chart = _repeat_count_svg(exported, assembly=True)
    assert '70 repeats' in chart and 'ESTIMATED' in chart and 'confidence 0' in chart


def test_rescued_reconstruction_precedes_depth_estimation(tmp_path, monkeypatch):
    locus, sequence, pairs, _ = reads(70)
    attempts = []
    def assemble(items, template, *args, **kwargs):
        attempts.append(len(items))
        assert not template.primer_only  # Anchor learning precedes the single assembly attempt.
        return sequence, [item.molecule_id for item in items], ''
    monkeypatch.setattr('mlvamaps.targeted_reconstruction.microassemble', assemble)
    calls, _, _, paths = run_reconstructed_fastq_inference(reads1=None, reads2=None,
        pairs=iter(pairs), loci=[locus], database_path=None, outdir=tmp_path,
        sample_id='s', technology='illumina', minimum_molecules=2, minimum_probability=.8)
    assert len(attempts) == 1
    assert calls[0]['inference_method'] == 'RECONSTRUCTED'
    assert calls[0]['status'] == 'called' and calls[0]['repeat_count'] == 70
    assert 'repeat_length_estimates' not in paths


def test_depth_fallback_retains_calibrated_likelihood_uncertainty(tmp_path, monkeypatch):
    locus, sequence, pairs, _ = reads(70)
    # Known flanks make likelihood coordinates identical to panel coordinates.
    locus = replace(locus, left_flank_sequence=sequence[20:90],
                    right_flank_sequence=sequence[-90:-20], repeat_motif=sequence[90:99])
    monkeypatch.setattr('mlvamaps.targeted_reconstruction.microassemble',
                        lambda *args, **kwargs: ('', [], 'assembly_node_limit'))
    monkeypatch.setattr('mlvamaps.targeted_reconstruction.candidate_likelihood',
        lambda *args, **kwargs: LocusRecovery(method='AMBIGUOUS', best=None,
            interval=(0, 100), states=np.array([0., 100.]),
            posterior=np.array([.5, .5]), log_likelihoods=np.zeros(2),
            reason='minimum_likelihood_tie', limit_reached=True))
    calls, _, _, paths = run_reconstructed_fastq_inference(reads1=None, reads2=None,
        pairs=iter(pairs), loci=[locus], database_path=None, outdir=tmp_path,
        sample_id='s', technology='illumina', minimum_molecules=2, minimum_probability=.8)
    assert calls[0]['inference_method'] == 'KMER_DEPTH'
    assert calls[0]['repeat_count'] == 70 and calls[0]['confidence'] == 0
    assert calls[0]['repeat_count_interval_kind'] == 'coverage_sensitivity'
    with paths['repeat_likelihoods'].open() as handle:
        likelihoods = list(csv.DictReader(handle, delimiter='\t'))
    assert [float(row['repeat_count']) for row in likelihoods] == [0, 100]
    assert [float(row['posterior']) for row in likelihoods] == [.5, .5]
    assert 'minimum_likelihood_tie' in calls[0]['reason']
    assert json.loads(paths['repeat_length_estimates'].read_text())['L']['status'] == 'diagnostic_only'


def test_one_covered_flank_cannot_supply_a_repeat_call():
    locus, _, pairs, items = reads(70)
    # Missing coverage remains visible in diagnostics, without a repeat call.
    left_reads = pairs[:180+90]
    recovery = LocusRecovery()
    diagnostics = estimate_graph_lengths([recovery], [locus], locus_templates([locus]),
                                         {'L': items}, lambda: iter(left_reads))
    assert recovery.product_size_bp is None and recovery.best is None
    assert diagnostics['L']['arm_depths'][1] == 0
    assert diagnostics['L']['status'] == 'diagnostic_only' and recovery.confidence == 0


def test_single_pair_depth_diagnostic_cannot_supply_a_repeat_call():
    locus, sequence, _, _ = reads(5)
    first, second = sequence[:150], sequence[-150:]
    pair = ReadPair('one', ReadRecord('one/1', first, 'I'*150),
                    ReadRecord('one/2', revcomp(second), 'I'*150))
    item = MoleculeEvidence('one', 'L', ('UNINFORMATIVE',), (first, second), ('+', '-'))
    recovery = LocusRecovery()
    diagnostics = estimate_graph_lengths([recovery], [locus], locus_templates([locus]),
                                         {'L': [item]}, lambda: iter([pair]))
    assert recovery.best is None and recovery.product_size_bp is None
    assert diagnostics['L']['amplicon_length'] > 0 and diagnostics['L']['status'] == 'diagnostic_only'
    assert diagnostics['L']['read_start_positions'] == [1, 1]
    assert recovery.confidence == 0 and not recovery.products


def test_native_depth_counts_match_scalar_across_boundaries_and_quality(tmp_path):
    from collections import Counter
    from mlvamaps.repeat_depth import _count_depth_chunk, _segments, _depth_targets, _depth_index
    rng = random.Random(912)
    sequences = [''.join(rng.choices('ACGT', k=150)) for _ in range(40)]
    sequences += ['A'*21, 'C'*14, 'G'*15, 'T'*21, '', 'N'*40,
                  'ACGT'*12+'N'+'TGCA'*12, 'acgt'*30]
    pairs, wanted = [], {k: set() for k in (15, 21)}
    for i, seq in enumerate(sequences):
        # Include absent targets, overlapping hits, missing mates and qualities.
        quality = ''.join('!' if j % 47 == 0 else 'I' for j in range(len(seq))) if i % 2 else None
        pairs.append(ReadPair(str(i), ReadRecord(str(i), seq, quality),
                              ReadRecord(str(i)+'b', revcomp(seq), quality[::-1] if quality else None)))
        for k, keys in wanted.items():
            keys.update(seq.upper()[j:j+k] for j in range(len(seq)-k+1)
                        if set(seq.upper()[j:j+k]) <= set('ACGT'))
    wanted = {k: tuple(sorted(keys | {'A'*k, 'T'*k})) for k, keys in wanted.items()}
    expected = Counter()
    pairs.append(ReadPair('orphan', ReadRecord('orphan', 'A'*70, 'I'*70)))
    for pair in pairs:
        for read in (pair.read1, pair.read2):
            if read is not None:
                for k, keys in wanted.items():
                    for segment in _segments(read.sequence.upper(), read.quality, k):
                        expected.update(segment[j:j+k] for j in range(len(segment)-k+1)
                                        if segment[j:j+k] in keys)
    observed = Counter()
    paths = tuple((k, str(tmp_path/f'{k}.npy')) for k in wanted)
    for k, path in paths:
        np.save(path, _depth_targets(wanted[k], k))
    for start in range(0, len(pairs), 7):
        size, counts = _count_depth_chunk((pairs[start:start+7], paths))
        assert size == len(pairs[start:start+7])
        for k, (positions, multiplicities) in counts.items():
            observed.update({wanted[k][i]: int(n) for i, n in zip(positions, multiplicities)})
    assert observed == expected
    size, counts = _count_depth_chunk(([], paths))
    assert size == 0 and all(not len(positions) for positions, _ in counts.values())
    assert all(isinstance(index, np.memmap) and not index.flags.writeable
               for index in _depth_index(paths).values())


def test_process_depth_graphs_and_counts_match_serial():
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing
    locus, _, pairs, items = reads(30)
    templates, evidence = locus_templates([locus]), {'L': items}
    serial, parallel = LocusRecovery(), LocusRecovery()
    expected = estimate_graph_lengths([serial], [locus], templates, evidence, lambda: iter(pairs))
    with ProcessPoolExecutor(max_workers=2, mp_context=multiprocessing.get_context('spawn')) as executor:
        observed = estimate_graph_lengths([parallel], [locus], templates, evidence, lambda: iter(pairs),
                                         locus_executor=executor, threads=2)
    assert observed == expected
    assert (parallel.best, parallel.interval, parallel.product_size_bp, parallel.reason) == (
        serial.best, serial.interval, serial.product_size_bp, serial.reason)
