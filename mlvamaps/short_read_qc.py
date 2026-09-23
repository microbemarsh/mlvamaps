"""Stream validated, filtered molecules with optional reusable FASTQ files."""
from contextlib import ExitStack
from itertools import islice
import time

from .io import _gzip_writer, open_text, read_fastq, read_fastq_pairs
from .models import ReadPair
from .progress import ProgressReporter
from .short_reads import qc_read_pairs


def filtered_pairs(reads1, reads2, filtered1, filtered2, orphans, *,
                   materialize, counters, statistics, min_length,
                   min_mean_quality, trim_quality, min_pair_retention,
                   sample_id, show_progress=False, chunk_size=5000,
                   decompression_threads=(0, 0)):
    """Filter once and feed recruitment while optionally writing its inputs.

    Paired survivors precede orphans, matching the historical disk round trip
    and preserving accumulation order in downstream mixture calculations.
    Only orphans must be spooled when no later stage needs the FASTQ files.
    ``seconds`` counts input/QC/write work, excluding time yielded to callers.
    The caller must close the generator if downstream inference fails.
    """
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    progress = ProgressReporter(enabled=show_progress)
    statistics.update(seconds=0.0, materialized=bool(materialize), gzip_writer=_gzip_writer.__name__)
    with ExitStack() as stack:
        first = stack.enter_context(open_text(filtered1, "wt")) if materialize else None
        second = stack.enter_context(open_text(filtered2, "wt")) if materialize else None
        orphan = stack.enter_context(open_text(orphans, "wt")) if reads2 else None
        iterator = read_fastq_pairs(reads1, reads2, decompression_threads=decompression_threads)
        stack.callback(iterator.close)
        while True:
            started = time.perf_counter()
            chunk = list(islice(iterator, chunk_size))
            if not chunk:
                statistics["seconds"] += time.perf_counter() - started
                break
            retained, metrics = qc_read_pairs(chunk, min_length, min_mean_quality,
                                              trim_quality, min_pair_retention)
            for key, value in metrics.items():
                counters[key] = counters.get(key, 0) + value
            ready = []
            for pair in retained:
                if reads2 and pair.read2 is None:
                    records = ((pair.read1, orphan),)
                else:
                    ready.append(pair)
                    records = ((pair.read1, first), (pair.read2, second))
                for record, handle in records:
                    if record is not None and handle is not None:
                        quality = record.quality or "I" * len(record.sequence)
                        handle.write(f"@{record.read_id}\n{record.sequence}\n+\n{quality}\n")
            statistics["seconds"] += time.perf_counter() - started
            progress.count(f"[{sample_id}] Input pairs filtered/validated", counters["input_pairs"])
            yield from ready
        # Finish gzip members before either orphan reading or classification.
    if reads2:
        for read in read_fastq(orphans):
            yield ReadPair(read.read_id, read)
    progress.step(f"[{sample_id}] QC finished: {counters.get('input_pairs', 0):,} input pairs; "
                  f"{statistics['seconds']:.1f}s spent reading, filtering and writing")
