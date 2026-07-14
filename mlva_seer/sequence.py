from __future__ import annotations

from collections import Counter
from threading import local
from typing import Optional, Tuple

import sassy


_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_SASSY_LOCAL = local()


def revcomp(sequence: str) -> str:
    return sequence.translate(_RC)[::-1].upper()


def mean_qscore(quality: Optional[str]) -> float:
    if not quality:
        return 0.0
    return sum(ord(ch) - 33 for ch in quality) / len(quality)


def hamming_distance(a: str, b: str) -> int:
    length = min(len(a), len(b))
    return sum(1 for i in range(length) if a[i] != b[i]) + abs(len(a) - len(b))


def _sassy_searcher(alphabet: str, rc: bool = False):
    searchers = getattr(_SASSY_LOCAL, "searchers", None)
    if searchers is None:
        searchers = {}
        _SASSY_LOCAL.searchers = searchers
    cache_key = (alphabet, rc, id(sassy))
    if cache_key in searchers:
        return searchers[cache_key]
    try:
        searcher = sassy.Searcher(alphabet, rc=rc)
    except TypeError:
        searcher = sassy.Searcher(alphabet)
    searchers[cache_key] = searcher
    return searcher


def _clear_sassy_searchers() -> None:
    """Clear searchers cached by the calling thread (primarily for tests)."""
    _SASSY_LOCAL.searchers = {}


def _find_best_sassy(pattern: str, sequence: str, max_mismatches: int) -> Tuple[Optional[int], Optional[int]]:
    searcher = _sassy_searcher("iupac")
    matches = searcher.search(pattern.encode(), sequence.encode(), k=max_mismatches)
    matches = [
        match
        for match in matches
        if all(
            base in "ACGT"
            for base in sequence[
                int(match.text_start) : int(
                    getattr(match, "text_end", int(match.text_start) + len(pattern))
                )
            ]
        )
    ]
    if not matches:
        return None, None
    best = min(matches, key=lambda match: (match.cost, match.text_start))
    return int(best.text_start), int(best.cost)


def find_best(pattern: str, sequence: str, max_mismatches: int) -> Tuple[Optional[int], Optional[int]]:
    pattern = pattern.upper()
    sequence = sequence.upper()
    if not pattern or len(pattern) > len(sequence):
        return None, None
    pos, distance = _find_best_sassy(pattern, sequence, max_mismatches)
    return pos, distance


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
