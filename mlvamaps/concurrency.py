from __future__ import annotations

import os

DEFAULT_THREADS = 32


def resolve_threads(threads: int | None) -> int:
    if threads is None or threads <= 0:
        return max(1, os.cpu_count() or 1)
    return max(1, threads)
