from __future__ import annotations

import sys
import time
from typing import TextIO
from contextlib import contextmanager


class ProgressReporter:
    def __init__(self, enabled: bool = True, stream: TextIO | None = None, min_interval: float = 1.0):
        self.enabled = enabled
        self.stream = stream or sys.stderr
        self.min_interval = min_interval
        self.started = time.monotonic()
        self._last_update: dict[str, float] = {}

    @contextmanager
    def phase(self, name: str, detail: str = ""):
        """Report a phase boundary and its own elapsed wall time."""
        started = time.monotonic()
        self.step(f"Phase {name} started{': ' + detail if detail else ''}")
        try:
            yield
        finally:
            self.step(f"Phase {name} finished in {time.monotonic() - started:.1f}s")

    def step(self, message: str) -> None:
        if not self.enabled:
            return
        elapsed = time.monotonic() - self.started
        print(f"[{elapsed:6.1f}s] {message}", file=self.stream, flush=True)

    def count(self, label: str, current: int, total: int | None = None, force: bool = False) -> None:
        if not self.enabled:
            return
        now = time.monotonic()
        last = self._last_update.get(label, 0.0)
        if not force and now - last < self.min_interval:
            return
        self._last_update[label] = now
        if total:
            pct = min(100.0, current / total * 100)
            self.step(f"{label}: {current:,}/{total:,} ({pct:.1f}%)")
        else:
            self.step(f"{label}: {current:,}")
