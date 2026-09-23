#!/usr/bin/env python3
"""Compare native FASTQ input/QC/output backends on the same bounded subset.

Run inside your allocation, optionally using --scratch-dir "$SLURM_TMPDIR".
Includes strict pair validation, QC, gzip writes, and evidence hashing, but no
recruitment or genotyping. Temporary filtered files are removed after each run.
"""
from contextlib import ExitStack, closing
import argparse
import gzip
import hashlib
from itertools import islice
import json
from pathlib import Path
import sys
import tempfile
import time

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from mlvamaps.io import read_fastq_pairs
from mlvamaps.short_reads import qc_read_pairs


def positive_int(value):
    result = int(value)
    if result < 1:
        raise argparse.ArgumentTypeError('must be positive')
    return result


def benchmark(args, backend, writer):
    decoder_threads = min(4, max(1, (args.threads - 1) // (2 if args.reads2 else 1) - 1))
    decoders = (decoder_threads, decoder_threads if args.reads2 else 0) if backend == 'rapidgzip' else (0, 0)
    digest, counters = hashlib.sha256(), {}
    with tempfile.TemporaryDirectory(prefix='mlvamaps-io-', dir=args.scratch_dir) as temp:
        started = time.perf_counter()
        with ExitStack() as stack:
            handles = [stack.enter_context(writer.open(Path(temp)/f'{i}.fq.gz', 'wb', compresslevel=1))
                       for i in range(3)]
            pairs = stack.enter_context(closing(read_fastq_pairs(args.reads1, args.reads2,
                                                                 decompression_threads=decoders)))
            selected = islice(pairs, args.pairs)
            while chunk := list(islice(selected, 5000)):
                retained, metrics = qc_read_pairs(chunk, 40, 15, 0, .5)
                for key, value in metrics.items():
                    counters[key] = counters.get(key, 0) + value
                for pair in retained:
                    records = ((pair.read1, 0), (pair.read2, 1)) if pair.read2 else ((pair.read1, 2),)
                    for read, index in records:
                        data = f'@{read.read_id}\n{read.sequence}\n+\n{read.quality}\n'.encode('ascii')
                        handles[index].write(data)
                        digest.update(bytes([index]))
                        digest.update(data)
        elapsed = time.perf_counter() - started
    return {'input_backend': backend, 'gzip_writer': writer.__name__, 'decoder_threads': decoders,
            'seconds': elapsed, 'qc': counters, 'retained_sha256': digest.hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reads1', required=True)
    parser.add_argument('--reads2')
    parser.add_argument('--pairs', type=positive_int, default=100000)
    parser.add_argument('--threads', type=positive_int, default=8)
    parser.add_argument('--scratch-dir')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from isal import igzip
    import rapidgzip  # Fail early with a clear missing dependency.
    if args.threads < (5 if args.reads2 else 3):
        parser.error('allocate at least 5 threads for paired input, or 3 for single input')
    results = []
    for backend, writer in [('htslib', gzip), ('htslib', igzip), ('rapidgzip', igzip)]:
        row = benchmark(args, backend, writer)
        print(json.dumps(row), flush=True)
        if results and (row['qc'], row['retained_sha256']) != (results[0]['qc'], results[0]['retained_sha256']):
            raise RuntimeError('QC or retained reads differ between I/O backends')
        results.append(row)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({'rapidgzip_version': rapidgzip.__version__, 'benchmarks': results}, indent=2) + '\n')


if __name__ == '__main__':
    main()
