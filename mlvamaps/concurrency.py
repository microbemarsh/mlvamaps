from __future__ import annotations

import os

DEFAULT_THREADS = 32


def resolve_threads(threads: int | None) -> int:
    if threads is None or threads <= 0:
        return max(1, os.cpu_count() or 1)
    return max(1, threads)


def bounded_ordered_map(executor, function, iterable, max_pending: int):
    """Preserve order while bounding both queued tasks and completed results.

    Unlike Executor.map on Python <3.14, this never consumes the whole input.
    Keeping completed futures in the same bounded deque also prevents a slow
    first chunk from allowing an unbounded backlog of later completed chunks.
    """
    from collections import deque
    from itertools import islice

    if max_pending < 1:
        raise ValueError("max_pending must be positive")
    iterator = iter(iterable)
    pending = deque(executor.submit(function, item) for item in islice(iterator, max_pending))
    while pending:
        yield pending.popleft().result()
        for item in islice(iterator, 1):
            pending.append(executor.submit(function, item))
