"""Optional native gzip acceleration, retaining HTSlib's FASTQ parser."""
from contextlib import contextmanager
from importlib.util import find_spec
import os
from pathlib import Path
import threading

import pysam

from .concurrency import resolve_threads


def short_read_io_plan(paths, threads, minimum_bytes=8 * 1024 * 1024):
    """Reserve CPUs for decoders, pipe feeders and QC before recruitment.

    Small gzip files use HTSlib directly to avoid decoder startup overhead.
    Explicit positive parallelism prevents rapidgzip using the entire host.
    """
    total = resolve_threads(threads)
    selected = [bool(path and Path(path).is_file() and Path(path).suffix == '.gz'
                     and Path(path).stat().st_size >= minimum_bytes) for path in paths]
    available = find_spec('rapidgzip') is not None and Path('/dev/fd').is_dir()
    # Two decoder threads and one feeding thread per input, plus the parent.
    reserved = 3 * sum(selected) + 1
    # Keep at least half the allocation for recruitment. In particular, an
    # eight-CPU batch sample must not lose seven CPUs to paired decompression.
    active = available and any(selected) and total >= max(8, 2 * reserved)
    decoders = [2 if active and select else 0 for select in selected]
    return {'input_backends': ['rapidgzip+htslib' if n else 'htslib' for n in decoders],
            'decoder_threads': decoders, 'feeder_threads': sum(n > 0 for n in decoders),
            'recruitment_workers': max(1, total - reserved) if active else total,
            'requested_threads': total}


@contextmanager
def rapidgzip_fastx(path, threads):
    """Pipe native parallel decoding into native parsing without a disk copy.

    A bounded OS pipe supplies backpressure. Closing the parser also unblocks
    the producer on early exit. Decoder failures propagate after pipe EOF;
    truncated input must not become a successful, shortened sample.
    """
    import rapidgzip
    if threads < 1:
        raise ValueError('rapidgzip decoder threads must be positive')
    read_fd, write_fd = os.pipe()
    failures = []

    def feed():
        try:
            with os.fdopen(write_fd, 'wb', buffering=0) as target:
                with rapidgzip.open(str(path), parallelization=threads) as source:
                    while block := source.read(1024 * 1024):
                        pending = memoryview(block)
                        while pending:
                            pending = pending[target.write(pending):]
        except BrokenPipeError:
            pass  # The caller deliberately stopped or rejected a FASTQ record.
        except Exception as exc:
            failures.append(exc)

    producer = threading.Thread(target=feed, name='mlvamaps-rapidgzip')
    producer.start()
    try:
        with pysam.FastxFile(f'/dev/fd/{read_fd}', persist=False) as reader:
            # FastxFile owns a separate descriptor. Keep no spare reader open:
            # it would prevent BrokenPipeError after an early parser close.
            os.close(read_fd)
            read_fd = None
            yield reader
        producer.join()
        if failures:
            raise ValueError(f'rapidgzip failed to decode {str(path)!r}: {failures[0]}') from failures[0]
    finally:
        if read_fd is not None:
            os.close(read_fd)
        producer.join()
