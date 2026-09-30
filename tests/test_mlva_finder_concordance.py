"""Compare FASTQ calls with an independent, pinned MLVA_finder oracle."""
import csv
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys

import pytest

from mlvamaps.io import write_fastq
from mlvamaps.models import ReadRecord
from mlvamaps.sequence import revcomp
from scripts.benchmark_cross_mode import compare_runs, load_mlva_finder


ORACLE = Path(__file__).parent/'data'/'mlva_finder_oracle'/'fastq_output.csv'


def make_oracle_inputs(root):
    """Deterministic, arbitrary DNA; reference alleles differ from every read allele."""
    rng = random.Random(554938)
    dna = lambda n: ''.join(rng.choices('ACGT', k=n))
    assembly = root/'assemblies'
    database = root/'database'
    assembly.mkdir(parents=True)
    database.mkdir()
    primer_rows, fasta, single, first, second = [], [], [], [], []
    for i, (count, extra, change) in enumerate([
            (7, 0, ''), (7, 1, ''), (7, 3, ''), (8, 0, 'substitution'),
            (6, 0, 'insertion'), (9, 0, 'reverse')], 1):
        left, right, motif = dna(60), dna(60), dna(4)
        name = f'T{i}_4bp_140bp_5U'
        primer_rows.append(f'{name}\t{left[:20]}\t{revcomp(right[-20:])}\n')
        (database/f'{name}.fasta').write_text(f'>reference\n{left}{motif*5}{right}\n')
        sequence = left+motif*count+motif[:extra]+right
        if change == 'substitution':
            sequence = sequence[:8]+next(b for b in 'ACGT' if b != sequence[8])+sequence[9:]
        elif change == 'insertion':
            sequence = sequence[:8]+'A'+sequence[8:]
        elif change == 'reverse':
            sequence = revcomp(sequence)
        fasta.append(f'>{name}\n{sequence}\n')
        for j in range(4):
            read_id = f'{name}_{j}'
            single.append(ReadRecord(read_id, sequence, 'I'*len(sequence)))
            first.append(ReadRecord(read_id+'/1', sequence[:100], 'I'*100))
            second.append(ReadRecord(read_id+'/2', revcomp(sequence[-100:]), 'I'*100))
    (assembly/'oracle.fasta').write_text(''.join(fasta))
    primers = root/'primers.tsv'
    primers.write_text(''.join(primer_rows))
    for name, records in [('single.fq', single), ('r1.fq', first), ('r2.fq', second)]:
        write_fastq(records, root/name)
    return primers, assembly/'oracle.fasta', database


@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('tolerance', [.25, 0])
def test_fastq_and_assembly_match_independent_mlva_finder_calls(tmp_path, paired, tolerance):
    from mlvamaps.assembly_call import run_assembly_call
    from mlvamaps.locus_reconstruction import run_reconstructed_fastq_inference
    from mlvamaps.primers import read_primer_pairs
    if shutil.which('minimap2') is None:
        pytest.skip('minimap2 is required for reference-backed FASTQ calling')
    primers, assembly, database = make_oracle_inputs(tmp_path)
    oracle_path = ORACLE if tolerance else ORACLE.with_name('fastq_unrounded_output.csv')
    expected = load_mlva_finder(oracle_path, 'oracle.fasta')
    assembly_paths = run_assembly_call(str(assembly), None, str(tmp_path/'assembly_calls'),
        's', primers_path=str(primers), threads=1, show_progress=False, assembly_round_tolerance=tolerance)
    with assembly_paths['calls'].open() as handle:
        assembly_calls = {row['locus_id']: row for row in csv.DictReader(handle, delimiter='\t')}
    calls, _, _, _ = run_reconstructed_fastq_inference(
        reads1=tmp_path/('r1.fq' if paired else 'single.fq'), reads2=tmp_path/'r2.fq' if paired else None,
        loci=read_primer_pairs(primers), database_path=database, outdir=tmp_path/'sr', sample_id='s',
        technology='illumina', minimum_molecules=0, minimum_probability=.8, threads=1,
        show_progress=False, round_tolerance=tolerance)
    assert len(calls) == len(expected) == 6
    for call in calls:
        oracle, assembled = expected[call['locus']], assembly_calls[call['locus']]
        assert call['status'] == 'called', call
        assert call['repeat_count'] == float(assembled['repeat_count']) == oracle['repeat_count']
        assert call['product_size_bp'] == float(assembled['product_size_bp']) == float(oracle['amplicon_length'])


@pytest.mark.parametrize('change', ['exact', 'wrong', 'missing', 'estimated'])
def test_external_comparison_requires_recovery_of_all_oracle_calls(tmp_path, change):
    expected = load_mlva_finder(ORACLE, 'oracle.fasta')
    sr = tmp_path/'sr'
    sr.mkdir()
    with (sr/'calls.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['locus_id', 'repeat_count', 'status'])
        for i, (locus, row) in enumerate(expected.items()):
            count, status = row['repeat_count'], 'PASS'
            if i == 0:
                if change == 'wrong': count += 1
                if change == 'missing': count, status = '', 'PRESENT_COUNT_UNKNOWN'
                if change == 'estimated': status = 'ESTIMATED'
            writer.writerow([locus, count, status])
    manifest = [{'sample_id': 's', 'mode': 'mlva_finder', 'outdir': ORACLE, 'source_sample': 'oracle.fasta'},
                {'sample_id': 's', 'mode': 'sr', 'outdir': sr}]
    summary, rows = compare_runs(manifest)
    comparison = summary['mlva_finder:sr']
    assert comparison['strict_repeat_match'] == (change == 'exact')
    assert comparison['exact_repeat_recovery_by_mode']['mlva_finder'] == (1 if change == 'exact' else 5/6)
    assert len(rows) == 6
    if change in {'missing', 'estimated'}:
        assert comparison['exact_repeat_concordance'] == 1  # Cannot hide the missing sixth call.
        assert comparison['missing_call_loci_by_mode']['sr'] == 1
    manifest_path = tmp_path/'runs.tsv'
    with manifest_path.open('w') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=['sample_id', 'mode', 'outdir', 'source_sample'])
        writer.writeheader()
        writer.writerows(manifest)
    output = tmp_path/'comparison.json'
    result = subprocess.run([sys.executable, 'scripts/benchmark_cross_mode.py', str(manifest_path),
                             '--output', str(output), '--require-exact'], capture_output=True, text=True)
    assert result.returncode == (0 if change == 'exact' else 1), result.stderr
    assert json.loads(output.read_text())['comparisons']['mlva_finder:sr'] == comparison
    assert output.with_suffix('.loci.tsv').is_file()


def test_external_comparison_rejects_wrong_sample_and_wide_fingerprint(tmp_path):
    with pytest.raises(ValueError, match='source_sample'):
        load_mlva_finder(ORACLE, 'wrong_sample')
    wide = tmp_path/'MLVA_analysis.csv'
    wide.write_text('key,Access_number,T1\n001,s,7\n')
    with pytest.raises(ValueError, match='detailed'):
        load_mlva_finder(wide, 's')


def test_external_comparison_retains_zero_and_missing_and_rejects_invalid_rows(tmp_path):
    output = tmp_path/'external_output.csv'
    output.write_text('strain,primer,allele,size\ns,zero,0,120\ns,missing,,\nother,zero,9,156\n')
    details = load_mlva_finder(output, 's')
    assert details['zero']['repeat_count'] == 0 and details['zero']['status'] == 'PASS'
    assert details['missing']['repeat_count'] is None and details['missing']['status'] == 'NOT_FOUND'
    for rows, message in [('s,L,7,148\ns,L,8,152\n', 'duplicate'),
                          ('s,L,NaN,148\n', 'invalid allele'), ('s,L,-1,148\n', 'invalid allele'),
                          ('s,L\n', 'incomplete')]:
        output.write_text('strain,primer,allele,size\n'+rows)
        with pytest.raises(ValueError, match=message):
            load_mlva_finder(output, 's')


def test_strict_comparison_detects_an_entire_missing_sample_run(tmp_path):
    run = tmp_path/'sr'
    run.mkdir()
    expected = load_mlva_finder(ORACLE, 'oracle.fasta')
    (run/'calls.tsv').write_text('locus_id\trepeat_count\tstatus\n'+''.join(
        f'{locus}\t{row["repeat_count"]}\tPASS\n' for locus, row in expected.items()))
    summary, _ = compare_runs([
        {'sample_id': 's', 'mode': 'sr', 'outdir': run},
        {'sample_id': 's', 'mode': 'mlva_finder', 'outdir': ORACLE, 'source_sample': 'oracle.fasta'},
        {'sample_id': 'missing', 'mode': 'mlva_finder', 'outdir': ORACLE, 'source_sample': 'oracle.fasta'}])
    result = summary['mlva_finder:sr']
    assert result['full_profile_concordance'] == 1
    assert result['missing_runs_by_mode']['sr'] == ['missing']
    assert not result['strict_repeat_match']
