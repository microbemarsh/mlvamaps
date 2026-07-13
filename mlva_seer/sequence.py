from __future__ import annotations

from collections import Counter
from threading import local
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
    if sassy is None:
        return None
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
    if searcher is None:
        return None, None
    matches = searcher.search(pattern.encode(), sequence.encode(), k=max_mismatches)
    if not matches:
        return None, None
    best = min(matches, key=lambda match: (match.cost, match.text_start))
    return int(best.text_start), int(best.cost)


def find_best_many(
    patterns: list[str],
    texts: list[str],
    max_mismatches: int,
    threads: int,
) -> dict[tuple[int, int, str], tuple[int, int]] | None:
    """Search all patterns and texts in Sassy using Rust-owned threading.

    The result maps ``(pattern index, text index, orientation)`` to the best
    ``(position, edit distance)``. A return value of ``None`` means the
    installed Sassy binding does not provide the batched API.
    """
    searcher = _sassy_searcher("iupac", rc=True)
    if searcher is None or not hasattr(searcher, "search_many"):
        return None
    encoded_patterns = [pattern.upper().encode() for pattern in patterns]
    encoded_texts = [text.upper().encode() for text in texts]
    matches = searcher.search_many(
        encoded_patterns,
        encoded_texts,
        k=max_mismatches,
        threads=max(1, threads),
        mode="batch_texts",
    )
    best_matches: dict[tuple[int, int, str], tuple[int, int]] = {}
    for match in matches:
        orientation = "reverse" if match.strand == "-" else "forward"
        if orientation == "reverse":
            position = len(encoded_texts[match.text_idx]) - int(match.text_end)
        else:
            position = int(match.text_start)
        key = (int(match.pattern_idx), int(match.text_idx), orientation)
        candidate = (position, int(match.cost))
        current = best_matches.get(key)
        if current is None or (candidate[1], candidate[0]) < (current[1], current[0]):
            best_matches[key] = candidate
    return best_matches


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
