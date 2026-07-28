from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Sequence

from .calling import allele_grid, estimate_repeat_count_from_spanning_read, gaussian_allele_probabilities, legacy_round_repeat_count, repeat_unit_length
from .models import AnchorMeasurement, Locus, LocusMeasurement
from .sequence import _sassy_searcher, mean_qscore, repeat_motif_statistics, revcomp


@dataclass(frozen=True)
class ExtractedReadInterval:
    start: int
    end: int
    sequence: str
    quality: str | list[int] | None
    orientation: str


_IUPAC = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"),
    "T": frozenset("T"), "R": frozenset("AG"), "Y": frozenset("CT"),
    "S": frozenset("CG"), "W": frozenset("AT"), "K": frozenset("GT"),
    "M": frozenset("AC"), "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"), "N": frozenset("ACGT"),
}


def _compatible(pattern_base: str, observed_base: str) -> bool:
    return observed_base in _IUPAC.get(pattern_base, frozenset(pattern_base))


_CIGAR_OPERATION_RE = re.compile(r"(\d+)([=XIDM])")


def _sassy_edit_counts(cigar: str) -> tuple[int, int, int]:
    """Translate a Sassy pattern-to-text CIGAR into read-relative edits."""
    mismatches = insertions = deletions = 0
    for length_text, operation in _CIGAR_OPERATION_RE.findall(cigar):
        length = int(length_text)
        if operation == "X":
            mismatches += length
        elif operation == "D":
            insertions += length
        elif operation == "I":
            deletions += length
    return mismatches, insertions, deletions


def find_anchor(pattern: str, sequence: str, max_edits: int = 3, search_start: int = 0, search_end: int | None = None) -> AnchorMeasurement | None:
    """Find a complete anchor with Sassy's native SIMD edit-distance search."""
    pattern = pattern.upper()
    sequence = sequence.upper()
    search_end = len(sequence) if search_end is None else min(search_end, len(sequence))
    if not pattern or search_start >= search_end:
        return None
    exact_start = sequence.find(pattern, max(0, search_start), search_end) if set(pattern) <= set("ACGT") else -1
    if exact_start >= 0:
        return AnchorMeasurement(exact_start, exact_start + len(pattern), 1.0, 0, 0, 0, 0, True)
    # Short flanks cannot tolerate the same absolute error count as a 20--25 bp
    # primer without becoming non-specific.
    max_edits = min(max_edits, max(1, len(pattern) // 5))
    offset = max(0, search_start)
    window = sequence[offset:search_end]
    searcher = _sassy_searcher("iupac")
    matches = searcher.search(
        pattern.encode("ascii"),
        window.encode("ascii"),
        k=max_edits,
    )
    candidates = []
    for match in matches:
        if (
            int(getattr(match, "pattern_start", 0)) != 0
            or int(getattr(match, "pattern_end", len(pattern))) != len(pattern)
        ):
            continue
        relative_start = int(match.text_start)
        relative_end = int(match.text_end)
        if relative_start < 0 or relative_end <= relative_start:
            continue
        observed = window[relative_start:relative_end]
        # The previous measurement engine treated ambiguity in the configured
        # anchor as IUPAC but did not treat N in a read/assembly as a wildcard.
        if any(base not in "ACGT" for base in observed):
            continue
        distance = int(match.cost)
        identity = (len(pattern) - distance) / max(len(pattern), len(observed))
        start = offset + relative_start
        end = offset + relative_end
        candidates.append(
            (
                distance,
                -identity,
                start,
                abs(len(observed) - len(pattern)),
                end,
                str(match.cigar),
            )
        )
    if not candidates:
        return None
    distance, _identity, start, _delta, end, cigar = min(candidates)
    mm, ins, dels = _sassy_edit_counts(cigar)
    identity = max(0.0, 1.0 - distance / max(len(pattern), end - start, 1))
    return AnchorMeasurement(start, end, round(identity, 6), distance, mm, ins, dels, True)


def reference_interval_to_query_interval(alignment, ref_start: int, ref_end: int) -> tuple[int, int] | None:
    """Project a half-open reference interval through CIGAR-aligned pairs."""
    if ref_start < 0 or ref_end <= ref_start:
        return None
    by_reference = {int(r): int(q) for q, r in alignment.get_aligned_pairs(matches_only=False) if q is not None and r is not None}
    first = by_reference.get(ref_start)
    last = by_reference.get(ref_end - 1)
    if first is None or last is None:
        return None
    start, end = sorted((first, last))
    return start, end + 1


def original_read_in_locus_orientation(alignment, sequence: str, quality: str | Sequence[int] | None = None) -> tuple[str, str | list[int] | None]:
    """Orient an original FASTQ record consistently with the mapped locus."""
    if not alignment.is_reverse:
        return sequence.upper(), quality if isinstance(quality, str) else (list(quality) if quality is not None else None)
    oriented_quality = quality[::-1] if isinstance(quality, str) else list(reversed(quality)) if quality is not None else None
    return revcomp(sequence), oriented_quality


def extract_reference_interval_from_original_read(
    alignment,
    sequence: str,
    quality: str | Sequence[int] | None,
    ref_start: int,
    ref_end: int,
    padding: int = 0,
) -> ExtractedReadInterval | None:
    """CIGAR-project and extract an interval from the original FASTQ read."""
    projected = reference_interval_to_query_interval(alignment, ref_start, ref_end)
    if projected is None:
        return None
    oriented_sequence, oriented_quality = original_read_in_locus_orientation(
        alignment, sequence, quality
    )
    start = max(0, projected[0] - max(0, padding))
    end = min(len(oriented_sequence), projected[1] + max(0, padding))
    extracted_quality = (
        oriented_quality[start:end] if oriented_quality is not None else None
    )
    return ExtractedReadInterval(
        start,
        end,
        oriented_sequence[start:end],
        extracted_quality,
        "reverse" if alignment.is_reverse else "forward",
    )


def _quality_text(qualities: str | Sequence[int] | None) -> str | None:
    if qualities is None:
        return None
    if isinstance(qualities, str):
        return qualities
    return "".join(chr(min(93, max(0, int(value))) + 33) for value in qualities)


def _credible_set(ranked: list[tuple[int | float, float]], mass: float = 0.95) -> tuple[int | float, ...]:
    selected, cumulative = [], 0.0
    for allele, probability in ranked:
        selected.append(allele)
        cumulative += probability
        if cumulative >= mass:
            break
    return tuple(selected)


def measure_locus_product(sequence: str, locus_model: Locus, qualities: str | Sequence[int] | None = None, source: str = "assembly", sequence_id: str | None = None, max_anchor_edits: int = 3, round_tolerance: float = 0.25, calibrated_product_size_bp: int | None = None) -> LocusMeasurement:
    """Measure a locus from independently verified anchors in a sequence."""
    sequence = sequence.upper()
    quality = _quality_text(qualities)
    forward = find_anchor(locus_model.forward_primer, sequence, max_anchor_edits)
    reverse = find_anchor(revcomp(locus_model.reverse_primer), sequence, max_anchor_edits, forward.end) if forward else None
    product_start = forward.start if forward else None
    product_end = reverse.end if reverse else None
    inner_start = forward.end if forward else 0
    inner_end = reverse.start if reverse else len(sequence)
    left = find_anchor(locus_model.left_flank_sequence, sequence, max_anchor_edits, inner_start, inner_end) if locus_model.left_flank_sequence else None
    right_start = left.end if left else inner_start
    right = find_anchor(locus_model.right_flank_sequence, sequence, max_anchor_edits, right_start, inner_end) if locus_model.right_flank_sequence else None
    unit_length = repeat_unit_length(locus_model) or None
    flank_repeat_length = right.start - left.end if left and right else None
    plausible_minimum = max(
        1,
        (max(0, locus_model.expected_min_repeats - 1) * unit_length)
        if unit_length else 1,
    )
    plausible_maximum = (
        (locus_model.expected_max_repeats + 1) * unit_length
        if unit_length else len(sequence)
    )
    flanks_resolved = bool(
        left and right and flank_repeat_length is not None
        and plausible_minimum <= flank_repeat_length <= plausible_maximum
    )
    repeat_start = left.end if flanks_resolved else inner_start if forward and reverse else None
    repeat_end = right.start if flanks_resolved else inner_end if forward and reverse else None
    if repeat_start is None or repeat_end is None or repeat_end < repeat_start:
        repeat_start = repeat_end = None
    repeat_sequence = sequence[repeat_start:repeat_end] if repeat_start is not None and repeat_end is not None else ""
    repeat_length = len(repeat_sequence) if repeat_start is not None else None
    raw_count = None
    method = ""
    if forward and reverse and repeat_length is not None:
        raw_count, method = estimate_repeat_count_from_spanning_read(locus_model, calibrated_product_size_bp or (reverse.end - forward.start), repeat_length, flanks_resolved=flanks_resolved)
    elif flanks_resolved and unit_length:
        raw_count, method = repeat_length / unit_length, "flank_bounded_repeat_length"
    motif_identity = None
    if repeat_sequence and locus_model.repeat_motif and unit_length:
        _parts, mismatches, _kmers = repeat_motif_statistics(repeat_sequence, locus_model.repeat_motif, unit_length)
        motif_identity = max(0.0, 1.0 - mismatches / max(len(repeat_sequence), 1))
    called = legacy_round_repeat_count(raw_count, round_tolerance) if raw_count is not None else None
    likelihoods = None
    confidence = None
    credible: tuple[int | float, ...] = ()
    second_allele = None
    second_probability = None
    if raw_count is not None:
        candidates = allele_grid(locus_model, step=0.5)
        anchor_indels = sum(anchor.insertions + anchor.deletions for anchor in (forward, reverse, left, right) if anchor is not None)
        qscore = mean_qscore(quality[repeat_start:repeat_end] if quality is not None and repeat_start is not None and repeat_end is not None else quality)
        sigma = max(0.08, (0.35 + 0.12 * anchor_indels + max(0.0, 15.0 - qscore) / 30.0) / max(math.sqrt(unit_length or 1), 1.0))
        probabilities = gaussian_allele_probabilities(raw_count, candidates, sigma)
        ranked = sorted(zip(candidates, probabilities), key=lambda item: (-item[1], float(item[0])))
        likelihoods = {allele: round(probability, 8) for allele, probability in ranked}
        confidence = ranked[0][1]
        credible = _credible_set(ranked)
        if len(ranked) > 1:
            second_allele, second_probability = ranked[1]
    if forward and reverse and raw_count is not None:
        status, failure = "FULL_PRODUCT", None
    elif flanks_resolved and raw_count is not None:
        status, failure = "REPEAT_INFORMATIVE", None
    elif forward or reverse or left or right:
        status, failure = "PRESENCE_ONLY", "insufficient anchors to measure the repeat interval"
    else:
        status, failure = "ANCHOR_FAILURE", "no complete locus anchors were identified"
    product_sequence = sequence[product_start:product_end] if product_start is not None and product_end is not None else ""
    product_quality = quality[product_start:product_end] if quality is not None and product_start is not None and product_end is not None else None
    return LocusMeasurement(locus_model.locus_id, sequence_id, source, product_start, product_end, forward.start if forward else None, forward.end if forward else None, reverse.start if reverse else None, reverse.end if reverse else None, repeat_start, repeat_end, repeat_length, unit_length, raw_count, called, forward.identity if forward else None, reverse.identity if reverse else None, motif_identity, likelihoods, confidence, status, failure, forward, reverse, credible, second_allele, second_probability, product_sequence, product_quality, repeat_sequence, method)
