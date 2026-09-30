#!/usr/bin/env python3
"""Compare existing mlvamaps runs; no dataset or organism assumptions.

Manifest TSV: sample_id, mode (assembly/sr/lr), outdir.
Run: python scripts/benchmark_cross_mode.py manifest.tsv --output comparison.json
Optional: --distance-matrix assembly=assembly_distances.tsv --distance-matrix sr=sr_distances.tsv
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np
import parasail

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mlvamaps.alignment import MASKED_DNA_MATRIX


def read_table(path):
    if not Path(path).is_file():
        return []
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def load_run(directory):
    if not (Path(directory) / "calls.tsv").is_file():
        raise ValueError(f"missing calls.tsv in {directory}")
    calls = {r['locus_id']: number(r.get('repeat_count')) for r in read_table(Path(directory)/'calls.tsv')
             if r.get('status') in {'PASS', 'LOW_DEPTH', 'MULTIPLE_VARIANTS', 'PRESENT'}}
    calls = {k: v for k, v in calls.items() if v is not None}
    variants = {}
    for row in read_table(Path(directory)/'reconstructed_locus_variants.tsv'):
        evidence = json.loads(row.get('evidence') or '{}')
        if evidence.get('meaningful') == 'no':
            continue
        variants.setdefault(row['locus_id'], []).append(row)
    return calls, variants


def snp_comparison(first, second):
    if not first or not second:
        return 0, 0
    alignment = parasail.nw_trace_striped_32(first, second, 5, 1, MASKED_DNA_MATRIX)
    pairs = [(a,b) for a,b in zip(alignment.traceback.query, alignment.traceback.ref)
             if a in 'ACGT-' and b in 'ACGT-' and (a,b) != ('-', '-')]
    return sum(a == b for a,b in pairs), len(pairs)


def compare_runs(manifest):
    runs, details = {}, {}
    for row in manifest:
        key = (row['sample_id'], row['mode'])
        if key in runs:
            raise ValueError(f'duplicate manifest entry: {key}')
        runs[key] = load_run(row['outdir'])
        directory = Path(row['outdir'])
        details[key] = {r['locus_id']: dict(r) for r in read_table(directory/'calls.tsv')}
        for qc in read_table(directory/'short_read_repeat_evidence.tsv'):
            details[key].setdefault(qc['locus_id'], {}).update(qc)
    modes = sorted({mode for _,mode in runs})
    samples = sorted({sample for sample,_ in runs})
    result, locus_rows = {}, []
    for a,b in itertools.combinations(modes, 2):
        differences, estimate_differences, exact_profiles, components = [], [], [], []
        snp_matches = snp_sites = 0
        repeat_matches = repeat_sites = 0
        callable_by_mode = {a: 0, b: 0}
        shared_by_sample = {}
        fraction_distances = []
        total_by_mode = {a: 0, b: 0}
        ambiguous_by_mode = {a: 0, b: 0}
        per_locus = {}
        for sample in samples:
            if (sample,a) not in runs or (sample,b) not in runs:
                continue
            ca, va = runs[sample,a]
            cb, vb = runs[sample,b]
            callable_by_mode[a] += len(ca)
            callable_by_mode[b] += len(cb)
            for mode in (a, b):
                total_by_mode[mode] += len(details[sample, mode])
                ambiguous_by_mode[mode] += sum(r.get('call_method') in {'AMBIGUOUS', 'NO_CALL'} or
                    r.get('status') in {'NOT_FOUND', 'UNRESOLVED', 'PRESENT_UNTYPED', 'PRESENT_COUNT_UNKNOWN', 'AMBIGUOUS', 'LOCUS_DROPOUT'}
                    for r in details[sample, mode].values())
            shared = ca.keys() & cb.keys()
            shared_by_sample[sample] = len(shared)
            exact_profiles.append(bool(ca) and ca == cb)
            for locus in sorted(details[sample,a].keys() | details[sample,b].keys()):
                row = {'sample_id': sample, 'mode_a': a, 'mode_b': b, 'locus_id': locus,
                       'repeat_a': ca.get(locus, ''), 'repeat_b': cb.get(locus, '')}
                for mode, suffix in ((a, 'a'), (b, 'b')):
                    info = details[sample, mode].get(locus, {})
                    row['best_estimate_' + suffix] = number(info.get('repeat_count'))
                    row['status_' + suffix] = info.get('status', '')
                if row['best_estimate_a'] is not None and row['best_estimate_b'] is not None:
                    delta = abs(row['best_estimate_a'] - row['best_estimate_b'])
                    row['best_estimate_absolute_difference'] = delta
                    estimate_differences.append(delta)
                rate = per_locus.setdefault(locus, {a: [0, 0], b: [0, 0]})
                for mode, calls, suffix in ((a, ca, 'a'), (b, cb, 'b')):
                    rate[mode][0] += int(locus in calls)
                    rate[mode][1] += int(locus in details[sample,mode])
                    info = details[sample,mode].get(locus, {})
                    for field in ('call_method', 'confidence', 'effective_depth', 'read_length', 'motif_length',
                                  'vntr_length', 'amplicon_length', 'amplicon_insert_ratio',
                                  'n_full_span', 'n_left_boundary', 'n_right_boundary', 'n_flank_pairs',
                                  'n_repeat_rich', 'n_anchored_repeat', 'n_softclip_left', 'n_softclip_right'):
                        row[field+'_'+suffix] = info.get(field, '')
                    read_length, vntr_length = number(info.get('read_length')), number(info.get('vntr_length'))
                    row['vntr_read_ratio_'+suffix] = vntr_length/read_length if read_length and vntr_length is not None else ''
                if locus not in shared:
                    row['comparison'] = 'missing_in_both' if locus not in ca and locus not in cb else 'missing_in_a' if locus not in ca else 'missing_in_b'
                    locus_rows.append(row)
                    continue
                delta = abs(ca[locus]-cb[locus])
                differences.append(delta)
                row.update(absolute_difference=delta, comparison='exact' if delta == 0 else 'discordant')
                if va.get(locus) and vb.get(locus):
                    pa = max(va[locus], key=lambda r: float(r['estimated_fraction']))
                    pb = max(vb[locus], key=lambda r: float(r['estimated_fraction']))
                    match, sites = snp_comparison(pa['snp_sequence'], pb['snp_sequence'])
                    snp_matches += match
                    snp_sites += sites
                    aa = {r['combined_marker']: float(r['estimated_fraction']) for r in va[locus]}
                    bb = {r['combined_marker']: float(r['estimated_fraction']) for r in vb[locus]}
                    components.append(aa.keys() == bb.keys())
                    fraction_distances.append(sum(abs(aa.get(k,0)-bb.get(k,0)) for k in aa.keys() | bb.keys()) / 2)
                    row['snp_matches'], row['snp_comparable_sites'] = match, sites
                    if delta == 0:
                        rm, rs = snp_comparison(pa.get('repeat_sequence', ''), pb.get('repeat_sequence', ''))
                        repeat_matches += rm
                        repeat_sites += rs
                        row['repeat_sequence_matches'], row['repeat_sequence_sites'] = rm, rs
                locus_rows.append(row)
        exact_matches = sum(delta == 0 for delta in differences)
        result[f'{a}:{b}'] = {
            'best_estimates_including_ambiguous': {
                'comparable_loci': len(estimate_differences),
                'exact_repeat_concordance': float(np.mean(np.asarray(estimate_differences)==0)) if estimate_differences else None,
                'within_one_repeat_concordance': float(np.mean(np.asarray(estimate_differences)<=1)) if estimate_differences else None,
            },
            'call_rate_by_mode': {m: callable_by_mode[m]/n if n else None for m,n in total_by_mode.items()},
            'ambiguous_no_call_rate_by_mode': {m: ambiguous_by_mode[m]/n if n else None for m,n in total_by_mode.items()},
            'incorrect_call_rate_among_comparable': (len(differences)-exact_matches)/len(differences) if differences else None,
            'per_locus_call_rate': {l: {m: called/total if total else None for m,(called,total) in values.items()}
                                    for l,values in per_locus.items()},
            'shared_callable_loci': len(differences), 'shared_callable_loci_by_sample': shared_by_sample,
            'callable_loci_by_mode': callable_by_mode,
            'exact_repeat_matches': exact_matches,
            'discordant_repeat_loci': len(differences) - exact_matches,
            # Each denominator includes calls missing from the other mode.
            # This prevents improved shared-call accuracy from hiding dropout.
            'exact_repeat_recovery_by_mode': {
                mode: exact_matches / count if count else None for mode, count in callable_by_mode.items()},
            'missing_call_loci_by_mode': {
                a: callable_by_mode[b] - len(differences), b: callable_by_mode[a] - len(differences)},
            'exact_repeat_concordance': float(np.mean(np.asarray(differences)==0)) if differences else None,
            'within_one_repeat_concordance': float(np.mean(np.asarray(differences)<=1)) if differences else None,
            'mean_absolute_repeat_difference': float(np.mean(differences)) if differences else None,
            'snp_concordance': snp_matches/snp_sites if snp_sites else None,
            'snp_comparable_sites': snp_sites,
            'repeat_sequence_concordance_at_equal_counts': repeat_matches/repeat_sites if repeat_sites else None,
            'repeat_sequence_comparable_sites': repeat_sites,
            'full_profile_concordance': float(np.mean(exact_profiles)) if exact_profiles else None,
            'component_set_concordance': float(np.mean(components)) if components else None,
            'mean_component_fraction_total_variation': float(np.mean(fraction_distances)) if fraction_distances else None,
        }
    return result, locus_rows


def ranks(values):
    values = np.asarray(values)
    result = np.empty(len(values))
    order = np.argsort(values, kind='stable')
    for value in np.unique(values):
        positions = np.flatnonzero(values[order] == value)
        result[order[positions]] = positions.mean()
    return result


def distance_correlations(paths):
    matrices = {}
    for item in paths:
        mode, path = item.split('=', 1)
        rows = read_table(path)
        matrix = {}
        for row in rows:
            id_field = next(iter(row))
            sample = row[id_field]
            for other, value in row.items():
                if other != id_field and sample < other and number(value) is not None:
                    matrix[sample,other] = number(value)
        matrices[mode] = matrix
    result = {}
    for a,b in itertools.combinations(sorted(matrices), 2):
        shared = sorted(matrices[a].keys() & matrices[b].keys())
        x = ranks([matrices[a][k] for k in shared])
        y = ranks([matrices[b][k] for k in shared])
        rho = float(np.corrcoef(x,y)[0,1]) if len(shared)>1 and np.std(x)>0 and np.std(y)>0 else None
        result[f'{a}:{b}'] = {'shared_distances': len(shared), 'spearman_rho': rho}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--distance-matrix', action='append', default=[], metavar='MODE=TSV')
    args = parser.parse_args()
    manifest = read_table(args.manifest)
    if not manifest or not {'sample_id','mode','outdir'} <= manifest[0].keys():
        parser.error('manifest needs sample_id, mode, outdir columns and at least one row')
    for row in manifest:
        path = Path(row['outdir'])
        row['outdir'] = str(path if path.is_absolute() else args.manifest.parent/path)
    summary, rows = compare_runs(manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({'comparisons': summary, 'distance_correlations': distance_correlations(args.distance_matrix)}, indent=2, allow_nan=False)+'\n')
    with args.output.with_suffix('.loci.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, delimiter='\t', fieldnames=list(dict.fromkeys(['sample_id','mode_a','mode_b','locus_id',
            'repeat_a','repeat_b','comparison','absolute_difference','snp_matches','snp_comparable_sites'] +
            [field for row in rows for field in row])))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
