from __future__ import annotations

from collections import Counter
from functools import lru_cache
from typing import Optional, Tuple

try:
    import sassy
except ImportError:  # pragma: no cover - depends on optional Rust extension
    sassy = None

try:
    import edlib
except ImportError:  # pragma: no cover - depends on optional C extension
    edlib = None


_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(sequence: str) -> str:
    return sequence.translate(_RC)[::-1].upper()


def mean_qscore(quality: Optional[str]) -> float:
    if not quality:
        return 0.0
    return sum(ord(ch) - 33 for ch in quality) / len(quality)


def hamming_distance(a: str, b: str) -> int:
    length = min(len(a), len(b))
    return sum(1 for i in range(length) if a[i] != b[i]) + abs(len(a) - len(b))


@lru_cache(maxsize=2)
def _sassy_searcher(alphabet: str):
    if sassy is None:
        return None
    try:
        return sassy.Searcher(alphabet, rc=False)
    except TypeError:
        searcher = sassy.Searcher(alphabet)
        return searcher


def _find_best_sassy(pattern: str, sequence: str, max_mismatches: int) -> Tuple[Optional[int], Optional[int]]:
    searcher = _sassy_searcher("iupac")
    if searcher is None:
        return None, None
    matches = searcher.search(pattern.encode(), sequence.encode(), k=max_mismatches)
    if not matches:
        return None, None
    best = min(matches, key=lambda match: (match.cost, match.text_start))
    return int(best.text_start), int(best.cost)


def _find_best_edlib(pattern: str, sequence: str, max_mismatches: int) -> Tuple[Optional[int], Optional[int]]:
    if edlib is None:
        return None, None
    result = edlib.align(pattern, sequence, mode="HW", task="locations", k=max_mismatches)
    distance = result.get("editDistance", -1)
    locations = result.get("locations") or []
    if distance >= 0 and locations:
        return locations[0][0], distance
    return None, None


def find_best(pattern: str, sequence: str, max_mismatches: int) -> Tuple[Optional[int], Optional[int]]:
    pattern = pattern.upper()
    sequence = sequence.upper()
    if not pattern or len(pattern) > len(sequence):
        return None, None
    pos, distance = _find_best_sassy(pattern, sequence, max_mismatches)
    if pos is not None:
        return pos, distance
    pos, distance = _find_best_edlib(pattern, sequence, max_mismatches)
    if pos is not None:
        return pos, distance
    best_pos = None
    best_mismatches = len(pattern) + 1
    window = len(pattern)
    for pos in range(0, len(sequence) - window + 1):
        mismatches = hamming_distance(pattern, sequence[pos : pos + window])
        if mismatches < best_mismatches:
            best_pos = pos
            best_mismatches = mismatches
            if mismatches == 0:
                break
    if best_mismatches <= max_mismatches:
        return best_pos, best_mismatches
    return None, None


def phred_error_probability(qscore: float) -> float:
    if qscore <= 0:
        return 0.2
    return min(0.5, 10 ** (-qscore / 10))


def majority_consensus(sequences: list[str]) -> str:
    if not sequences:
        return ""
    max_len = max(len(seq) for seq in sequences)
    consensus = []
    for idx in range(max_len):
        bases = [seq[idx] for seq in sequences if idx < len(seq)]
        if not bases:
            continue
        consensus.append(Counter(bases).most_common(1)[0][0])
    return "".join(consensus)


def identity(a: str, b: str) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    shared = min(len(a), len(b))
    matches = sum(1 for i in range(shared) if a[i] == b[i])
    return matches / max(len(a), len(b))
