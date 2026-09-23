#!/usr/bin/env python3
"""Time short-read recruitment on a bounded subset of real FASTQs.

Example (run in the installed mlvamaps environment):
  python scripts/benchmark_short_read_recruitment.py --reads1 R1.fq.gz \
    --reads2 R2.fq.gz --primers panel.tsv --database reference_build \
    --pairs 10000 --threads 1 4 8 --output recruitment_timing.json

This benchmarks recruitment, not QC, repeat inference, classification or reports.
Every worker count must produce the same retained evidence fingerprint. Run this
inside an allocation large enough for the largest requested worker count.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
from itertools import islice
import json
from pathlib import Path
import time
import sys

# Direct execution should benchmark this checkout, including spawned workers,
# rather than require an editable install or accidentally use an older package.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mlvamaps.io import read_fastq_pairs
from mlvamaps.locus_reconstruction import locus_templates
from mlvamaps.primers import read_loci_or_primers
from mlvamaps.repeat_likelihood import _cyclic_search_text, _flank_profile, _motif_phases, _primitive_motif, _read_profile, flank_hit, repetitive
from mlvamaps.short_read_recruitment import recruit_short_reads
from mlvamaps.short_reads import qc_read_pairs


def positive_int(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError('must be positive')
    return parsed


def benchmark(pairs, templates, threads):
    # Prevent an earlier serial run from warming the next one's caches.
    for function in (_cyclic_search_text, _motif_phases, _primitive_motif, _read_profile, _flank_profile, flank_hit, repetitive):
        function.cache_clear()
    digest = hashlib.sha256()
    statistics = {}
    unique = ambiguous = 0
    started = time.perf_counter()
    for chunk in recruit_short_reads(iter(pairs), templates, 'benchmark', threads, statistics=statistics):
        unique += len(chunk.evidence)
        ambiguous += len(chunk.ambiguous)
        # Hash one record at a time to keep audit memory bounded.
        for kind, rows in (('unique', (asdict(e) for e in chunk.evidence)), ('ambiguous', chunk.ambiguous)):
            for row in rows:
                digest.update(kind.encode())
                digest.update(json.dumps(row, sort_keys=True).encode())
        for length in chunk.insert_lengths:
            digest.update(repr(length).encode())
    elapsed = time.perf_counter() - started
    return {'requested_threads': threads, 'workers': statistics['workers'],
            'pairs': len(pairs), 'seconds': elapsed, 'pairs_per_second': len(pairs)/elapsed,
            'uniquely_recruited_pairs': unique, 'ambiguous_locus_rows': ambiguous,
            'evidence_sha256': digest.hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reads1', required=True)
    parser.add_argument('--reads2')
    panel = parser.add_mutually_exclusive_group(required=True)
    panel.add_argument('--primers')
    panel.add_argument('--loci')
    parser.add_argument('--database')
    parser.add_argument('--pairs', type=positive_int, default=10000)
    parser.add_argument('--threads', type=positive_int, nargs='+', default=[1, 4])
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    loci = read_loci_or_primers(args.loci, args.primers)
    templates = locus_templates(loci, args.database)
    pairs, qc = qc_read_pairs(list(islice(read_fastq_pairs(args.reads1, args.reads2), args.pairs)),
                             min_length=40, min_mean_quality=15, trim_quality=0, min_pair_retention=.5)
    results = []
    for threads in args.threads:
        row = benchmark(pairs, templates, threads)
        print(json.dumps(row), flush=True)
        if results and row['evidence_sha256'] != results[0]['evidence_sha256']:
            raise RuntimeError('recruitment evidence changed with worker count')
        results.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({'qc': qc, 'recruitment_benchmarks': results}, indent=2)+'\n')


if __name__ == '__main__':
    main()
