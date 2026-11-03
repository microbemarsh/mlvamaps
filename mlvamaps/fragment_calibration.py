"""Bounded exact-overlap fragment calibration, conditional on overlap detection."""

from __future__ import annotations

from collections import Counter
import hashlib
import heapq
import statistics

from .models import ReadPair
from .short_reads import merge_read_pair


class OverlapFragmentCollector:
    """Sample original molecules without retaining the complete FASTQ stream.

    A narrow overlapping component does not establish the whole library's
    insert distribution: longer, nonoverlapping fragments remain unobserved.
    """

    def __init__(self, max_pairs: int = 5000):
        if isinstance(max_pairs, bool) or not isinstance(max_pairs, int) or max_pairs < 101:
            raise ValueError("overlap calibration requires a sample limit of at least 101 pairs")
        self.max_pairs = max_pairs
        self._pairs: dict[str, ReadPair | None] = {}
        self._heap: list[tuple[int, str]] = []
        self._records = 0
        self._limited = False

    def add(self, pair: ReadPair) -> None:
        self._records += 1
        if pair.read2 is None:
            return
        identifier = pair.molecule_id
        if identifier in self._pairs:
            previous = self._pairs[identifier]
            if previous is not None and (
                previous.read1.sequence != pair.read1.sequence
                or previous.read2.sequence != pair.read2.sequence
            ):
                self._pairs[identifier] = None
            return
        # ponytail: a stable bottom-hash sample bounds overlap work and storage
        # at max_pairs; larger libraries improve precision by raising this cap.
        rank = int.from_bytes(hashlib.blake2b(identifier.encode(), digest_size=16).digest(), "big")
        if len(self._pairs) == self.max_pairs:
            self._limited = True
            if rank >= -self._heap[0][0]:
                return
            _, removed = heapq.heapreplace(self._heap, (-rank, identifier))
            del self._pairs[removed]
        else:
            heapq.heappush(self._heap, (-rank, identifier))
        self._pairs[identifier] = pair

    def finish(self) -> tuple[dict | None, dict]:
        """Return a guarded empirical component and explicit selection limits."""
        diagnostics = dict(
            calibration_status="insufficient_exact_overlaps", source="exact_overlap_component",
            whole_library=False, unobserved_tail_excluded=False,
            input_pair_records=self._records, sampled_fragments=len(self._pairs),
            sample_limit=self.max_pairs, sample_limited=self._limited,
            required_fragments=101, eligible_fragments=0, retained_fragments=0,
            conflicting_molecule_ids=sum(pair is None for pair in self._pairs.values()),
            minimum_overlap_bp=12,
        )
        observations = []
        for identifier in sorted(self._pairs):
            pair = self._pairs[identifier]
            if pair is None:
                continue
            merged = merge_read_pair(pair, min_overlap=12, max_mismatch_fraction=0, require_unique=True)
            if merged is None:
                continue
            first, second = len(pair.read1.sequence), len(pair.read2.sequence)
            length = len(merged.sequence)
            overlap = pair.read1.sequence[-(first + second - length):]
            if (set(overlap) - set("ACGT") or len(set(overlap)) < 3
                    or max(Counter(overlap).values()) > .8 * len(overlap)
                    or any(all(base == overlap[index % period] for index, base in enumerate(overlap))
                           for period in range(1, min(10, len(overlap) // 2) + 1))):
                continue
            observations.append((length, max(first, second), first + second - 12))
        lengths = sorted(row[0] for row in observations)
        diagnostics["eligible_fragments"] = len(lengths)
        if len(lengths) < 101:
            return None, diagnostics

        center = statistics.median(lengths)
        mad = statistics.median(abs(length - center) for length in lengths)
        retained = [row for row in observations if abs(row[0] - center) <= max(1, 4.5 * mad)]
        diagnostics.update(retained_fragments=len(retained), median_bp=center,
                           trimmed_fraction=1 - len(retained) / len(lengths))
        # A small error tail may be trimmed; a second library component must
        # not be discarded to manufacture a narrow calibration.
        if len(retained) < 101 or len(retained) / len(lengths) < .95:
            diagnostics["calibration_status"] = "heterogeneous_overlap_lengths"
            return None, diagnostics
        kept_lengths = [row[0] for row in retained]
        mean = statistics.fmean(kept_lengths)
        deviation = max(1.0, statistics.stdev(kept_lengths))
        diagnostics.update(mean_bp=mean, sd_bp=deviation, sd_floor_bp=1.0)
        if deviation > max(3, .03 * center):
            diagnostics["calibration_status"] = "broad_overlap_lengths"
            return None, diagnostics
        tail_size = max(1, int(.05 * len(lengths)))
        if any(lengths[index] - lengths[index - 1] > max(3, deviation)
               for index in range(tail_size, len(lengths) - tail_size + 1)):
            diagnostics["calibration_status"] = "multimodal_overlap_lengths"
            return None, diagnostics
        # Read lengths bound both detectable fragment ends. Four SD of room
        # is required, including the 1 bp error floor for a zero-MAD library.
        margin = min(min(mean - lower, upper - mean) for _, lower, upper in retained)
        diagnostics["detection_boundary_margin_bp"] = margin
        if margin < 4 * deviation:
            diagnostics["calibration_status"] = "overlap_edge_censoring"
            return None, diagnostics

        diagnostics["calibration_status"] = "calibrated_overlap_component"
        calibration = dict(
            mean_bp=mean, sd_bp=deviation, retained_fragments=len(retained),
            source="exact_overlap_component", whole_library=False,
            unobserved_tail_excluded=False, observed_min_bp=min(kept_lengths),
            observed_max_bp=max(kept_lengths),
        )
        return calibration, diagnostics
