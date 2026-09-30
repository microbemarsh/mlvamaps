from __future__ import annotations

from collections import Counter
from threading import local
from typing import Optional, Tuple

import numpy as np
from .sassy_cli import Searcher


_RC = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")
_SASSY_LOCAL = local()


def revcomp(sequence: str) -> str:
    return sequence.translate(_RC)[::-1].upper()


def mean_qscore(quality: Optional[str]) -> float:
    if not quality:
        return 0.0
    values = np.frombuffer(quality.encode("ascii"), dtype=np.uint8)
    return float(values.mean(dtype=np.float64) - 33.0)


def hamming_distance(a: str, b: str) -> int:
    length = min(len(a), len(b))
    if not length:
        return abs(len(a) - len(b))
    left = np.frombuffer(a[:length].encode("ascii"), dtype=np.uint8)
    right = np.frombuffer(b[:length].encode("ascii"), dtype=np.uint8)
    return int(np.count_nonzero(left != right)) + abs(len(a) - len(b))


def repeat_motif_statistics(
    sequence: str, motif: str, motif_length: int
) -> tuple[list[str], int, int]:
    """Compute motif comparisons in NumPy's native loops, retaining labels."""
    full_chunks = len(sequence) // motif_length
    full_length = full_chunks * motif_length
    if full_chunks:
        bases = np.frombuffer(
            sequence[:full_length].encode("ascii"), dtype=np.uint8
        ).reshape(full_chunks, motif_length)
        motif_bases = np.frombuffer(motif.encode("ascii"), dtype=np.uint8)
        shared_length = min(motif_length, len(motif))
        equal = bases[:, :shared_length] == motif_bases[:shared_length]
        exact_chunks = (
            np.all(equal, axis=1)
            if motif_length == len(motif)
            else np.zeros(full_chunks, dtype=np.bool_)
        )
        mismatches = int(equal.size - np.count_nonzero(equal))
        mismatches += full_chunks * abs(motif_length - len(motif))
        motif_kmers = int(np.count_nonzero(exact_chunks))
        pattern_parts = [
            motif if bool(exact_chunks[idx]) else sequence[idx * motif_length : (idx + 1) * motif_length]
            for idx in range(full_chunks)
        ]
    else:
        pattern_parts = []
        mismatches = 0
        motif_kmers = 0
    if full_length < len(sequence):
        pattern_parts.append(f"{sequence[full_length:]}:partial")
    return pattern_parts, mismatches, motif_kmers


def _sassy_searcher(alphabet: str, rc: bool = False):
    searchers = getattr(_SASSY_LOCAL, "searchers", None)
    if searchers is None:
        searchers = {}
        _SASSY_LOCAL.searchers = searchers
    cache_key = (alphabet, rc, id(Searcher))
    if cache_key in searchers:
        return searchers[cache_key]
    try:
        searcher = Searcher(alphabet, rc=rc)
    except TypeError:
        searcher = Searcher(alphabet)
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
