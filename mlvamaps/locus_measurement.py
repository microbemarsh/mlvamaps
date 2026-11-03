from __future__ import annotations

import math
import re
import regex
from dataclasses import dataclass
from functools import lru_cache
from typing import Sequence

from .calling import assembly_equivalent_product_allele, gaussian_allele_probabilities, legacy_round_repeat_count, repeat_unit_length
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
    concrete = set(pattern) <= set("ACGT")
    exact_start = sequence.find(pattern, max(0, search_start), search_end) if concrete else -1
    if exact_start >= 0:
        return AnchorMeasurement(exact_start, exact_start + len(pattern), 1.0, 0, 0, 0, 0, True,
                                 sequence.find(pattern, exact_start + 1, search_end) >= 0)
    # Short flanks cannot tolerate the same absolute error count as a 20--25 bp
    # primer without becoming non-specific.
    max_edits = min(max_edits, max(1, len(pattern) // 5))
    offset = max(0, search_start)
    window = sequence[offset:search_end]
    if len(window) < len(pattern) - max_edits:
        return None
    # With at most k edits, one of k+1 disjoint pattern pieces must survive
    # unchanged. This only rejects impossible matches; Sassy still resolves
    # every possible match, including indels and ambiguous primer bases.
    pieces = max_edits + 1
    if concrete and 0 < pieces <= len(pattern) and not any(
        pattern[index * len(pattern) // pieces:(index + 1) * len(pattern) // pieces] in window
        for index in range(pieces)
    ):
        return None
    # Short windows are cheaper to reject in-process than with a Sassy launch.
    # This checks existence only; Sassy still determines every accepted match
    # and its traceback. Larger windows retain the native SIMD search directly.
    if concrete and len(window) <= 512 and not regex.search(
        f"(?:{pattern}){{e<={max_edits}}}", window,
    ):
        return None
    searcher = _sassy_searcher("iupac")
    # Small products share one native-search cache entry across the ordered
    # primer/flank subwindows. Filter boundaries after the search.
    search_text, search_offset = (sequence, 0) if len(sequence) <= 512 else (window, offset)
    matches = searcher.search(
        pattern.encode("ascii"),
        search_text.encode("ascii"),
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
        if search_offset + relative_start < offset or search_offset + relative_end > search_end:
            continue
        observed = search_text[relative_start:relative_end]
        # The previous measurement engine treated ambiguity in the configured
        # anchor as IUPAC but did not treat N in a read/assembly as a wildcard.
        if any(base not in "ACGT" for base in observed):
            continue
        distance = int(match.cost)
        identity = (len(pattern) - distance) / max(len(pattern), len(observed))
        start = search_offset + relative_start
        end = search_offset + relative_end
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
    ambiguous = any(candidate[0] == distance and (candidate[2], candidate[4]) != (start, end)
                    for candidate in candidates)
    return AnchorMeasurement(start, end, round(identity, 6), distance, mm, ins, dels, True, ambiguous)


def reference_interval_to_query_interval(alignment, ref_start: int, ref_end: int) -> tuple[int, int] | None:
    """Project a half-open reference interval through CIGAR-aligned pairs."""
    if ref_start < 0 or ref_end <= ref_start:
        return None
    reference, query = alignment.reference_start, 0
    endpoints = {}
    for operation, length in alignment.cigartuples or ():
        if operation in (0, 7, 8):
            for boundary in (ref_start, ref_end - 1):
                if reference <= boundary < reference + length:
                    endpoints[boundary] = query + boundary - reference
        if operation in (0, 1, 4, 7, 8):
            query += length
        if operation in (0, 2, 3, 7, 8):
            reference += length
    first = endpoints.get(ref_start)
    last = endpoints.get(ref_end - 1)
    if first is None or last is None:
        return None
    return (first, last + 1) if first <= last else None


def mapped_locus_window(row, context, *, refine=True):
    """Locate observed locus sequence using CIGAR and locally verified flanks.

    Coordinates refer to the reference-oriented original read. CIGAR only
    locates search windows; it never supplies repeat length or missing bases.
    """
    return _mapped_window(row.query_sequence, row.query_quality, row.reverse,
                          row.reference_start, row.reference_end, row.query_start, row.query_end,
                          row.cigar, context.sequence, context.repeat_start, context.repeat_end, refine)


@lru_cache(maxsize=512)
def _mapped_window(original_sequence, original_quality, reverse, reference_start, reference_end,
                   query_start, query_end, cigar, reference_sequence, repeat_start, repeat_end, refine):
    # All sequence/quality/coordinate inputs participate in this bounded cache;
    # identical reads reuse alignment work while retaining distinct molecule IDs.
    sequence = revcomp(original_sequence) if reverse else original_sequence
    quality = original_quality[::-1] if reverse and original_quality else original_quality
    padding = 24
    start = max(0, query_start - reference_start - padding)
    end = min(len(sequence), query_end + len(reference_sequence) - reference_end + padding)
    if not refine:
        return sequence[start:end], quality[start:end] if quality else None, start, False
    # Project boundaries in O(CIGAR operations), without allocating aligned pairs.
    reference, query = reference_start, 0
    positions = {}
    for size, operation in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
        size = int(size)
        if operation in "M=X":
            for boundary in (repeat_start, repeat_end):
                if reference <= boundary <= reference + size:
                    positions.setdefault(boundary, query + boundary - reference)
        if operation in "MIS=X":
            query += size
        if operation in "MDN=X":
            reference += size
    anchors = []
    for left, boundary in ((True, repeat_start), (False, repeat_end)):
        position = positions.get(boundary)
        flank = reference_sequence[:boundary] if left else reference_sequence[boundary:]
        anchor = None
        if position is not None and len(flank) >= 21:
            # ponytail: three bounded native searches, up to 75 flank bases.
            # Wider genomic context would require an extended database schema.
            for length in sorted({min(len(flank), size) for size in (75, 35, 21)}, reverse=True):
                pattern = flank[-length:] if left else flank[:length]
                anchor = find_anchor(pattern, sequence, max(3, length // 10),
                                     max(start, position - (length if left else 0) - padding),
                                     min(end, position + (0 if left else length) + padding))
                if anchor:
                    break
        anchors.append(anchor)
    left, right = anchors
    spanning = False
    if left and right and left.end <= right.start:
        product_start = left.start - (repeat_start - (left.end - left.start))
        product_end = right.end + len(reference_sequence) - repeat_end - (right.end - right.start)
        spanning = product_start >= 0 and product_end <= len(sequence) and "H" not in cigar
        start, end = max(start, product_start - padding), min(end, product_end + padding)
    return sequence[start:end], quality[start:end] if quality else None, start, spanning


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
    if start >= end or projected[1] > len(oriented_sequence):
        return None
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


@lru_cache(maxsize=512)
def _reference_measurement(sequence: str, locus: Locus) -> LocusMeasurement:
    return measure_locus_product(sequence, locus, round_tolerance=0)


def reference_product_allele(sequence: str, locus: Locus, round_tolerance: float = 0.25):
    """Use primer-adjusted reference length on the assembly allele scale.

    Only an unambiguous product covering the complete reference permits primer
    correction. Other reference records retain their whole-length convention.
    Reference anchor searches are shared across calibration and inference.
    """
    measured = _reference_measurement(sequence, locus)
    if (measured.status == "FULL_PRODUCT" and measured.product_start == 0
            and measured.product_end == len(sequence)
            and not measured.forward_anchor.ambiguous and not measured.reverse_anchor.ambiguous):
        raw = measured.raw_repeat_count
        return raw, legacy_round_repeat_count(raw, round_tolerance) if raw is not None else None
    return assembly_equivalent_product_allele(locus, len(sequence), round_tolerance)


def reference_nonrepeat_lengths(loci, contexts):
    """Calibrate flank spans only when reference products agree on nonrepeat size.

    Configured flanks can be short boundary anchors, so their lengths alone
    cannot reconstruct the primer-to-primer nonrepeat length. Reference repeat
    counts and copy numbers do not determine this calibration.
    """
    models = {locus.locus_id: locus for locus in loci
              if repeat_unit_length(locus) and locus.left_flank_sequence and locus.right_flank_sequence}
    lengths = {name: set() for name in models}
    for context in contexts:
        locus = models.get(context.locus_id)
        if locus is None:
            continue
        measured = _reference_measurement(context.sequence, locus)
        anchors = (measured.forward_anchor, measured.reverse_anchor)
        if (measured.status != "FULL_PRODUCT" or not measured.metadata.get("flanks_resolved")
                or measured.metadata.get("boundary_ambiguous")
                or any(anchor is None or anchor.ambiguous for anchor in anchors)):
            lengths[locus.locus_id].add(None)
            continue
        lengths[locus.locus_id].add(measured.reverse_anchor.start - measured.forward_anchor.end
            + len(locus.forward_primer) + len(locus.reverse_primer) - measured.repeat_length_bp)
    return {name: tuple(sorted(values)) if None not in values else () for name, values in lengths.items()}


def measure_locus_product(sequence: str, locus_model: Locus, qualities: str | Sequence[int] | None = None, source: str = "assembly", sequence_id: str | None = None, max_anchor_edits: int = 3, round_tolerance: float = 0.25, calibrated_product_size_bp: int | None = None, reference_nonrepeat_bp: tuple[int, ...] | None = None) -> LocusMeasurement:
    """Measure a locus from independently verified anchors in a sequence."""
    sequence = sequence.upper()
    quality = _quality_text(qualities)
    if quality is not None and len(quality) != len(sequence):
        raise ValueError("Locus sequence and quality lengths differ")
    if len(sequence) <= 512:
        # Batch only plausible non-exact anchors; exact products launch no
        # subprocess. This reuses the existing Sassy batch/cache implementation.
        pending = {}
        for pattern in (locus_model.forward_primer, revcomp(locus_model.reverse_primer),
                        locus_model.left_flank_sequence, locus_model.right_flank_sequence):
            edits = min(max_anchor_edits, max(1, len(pattern) // 5))
            if (pattern and set(pattern) <= set("ACGT") and pattern not in sequence
                    and regex.search(f"(?:{pattern}){{e<={edits}}}", sequence)):
                pending.setdefault(edits, []).append(pattern.encode("ascii"))
        if pending:
            prime = getattr(_sassy_searcher("iupac"), "prime_many", None)
            if prime:
                for edits, patterns in pending.items():
                    prime(patterns, [sequence.encode("ascii")], edits)
    forward = find_anchor(locus_model.forward_primer, sequence, max_anchor_edits)
    reverse = find_anchor(revcomp(locus_model.reverse_primer), sequence, max_anchor_edits, forward.end) if forward else None
    if forward and forward.ambiguous and reverse:
        # An alternative start beyond the selected end cannot form this pair.
        forward = find_anchor(locus_model.forward_primer, sequence, max_anchor_edits, 0, reverse.start)
    product_start = forward.start if forward else None
    product_end = reverse.end if reverse else None
    inner_start = forward.end if forward else 0
    inner_end = reverse.start if reverse else len(sequence)
    left_pattern, right_pattern = locus_model.left_flank_sequence, locus_model.right_flank_sequence
    left = find_anchor(left_pattern, sequence, max_anchor_edits, inner_start, inner_end) if left_pattern else None
    right_start = left.end if left else inner_start
    right = find_anchor(right_pattern, sequence, max_anchor_edits, right_start, inner_end) if right_pattern else None
    partial_flanks = False
    if reference_nonrepeat_bp is not None and len(reference_nonrepeat_bp) == 1:
        # A read endpoint can truncate a long configured flank while retaining
        # its repeat boundary. Require independently calibrated nonrepeat length
        # and retain the complete-flank search whenever it succeeds.
        if left is None and forward is None and len(left_pattern) > 21:
            candidate = find_anchor(left_pattern[-21:], sequence, min(2, max_anchor_edits),
                                    inner_start, right.start if right else inner_end)
            if candidate and candidate.end - inner_start < len(left_pattern):
                left, left_pattern, partial_flanks = candidate, left_pattern[-21:], True
        if right is None and len(right_pattern) > 21:
            candidate = find_anchor(right_pattern[:21], sequence, min(2, max_anchor_edits),
                                    left.end if left else inner_start, inner_end)
            if (candidate and inner_end - candidate.start < len(right_pattern)
                    and not find_anchor(revcomp(locus_model.reverse_primer), sequence,
                                        max_anchor_edits, candidate.end, inner_end)):
                right, right_pattern, partial_flanks = candidate, right_pattern[:21], True
    if left and left.ambiguous and right:
        left = find_anchor(left_pattern, sequence, min(2, max_anchor_edits) if left_pattern != locus_model.left_flank_sequence else max_anchor_edits,
                           inner_start, right.start)
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
    strong_flanks = all(anchor and not anchor.ambiguous and anchor.identity >= .9
                        and len(pattern) >= 21 and len(set(pattern)) >= 3
                        and not set(pattern) - set("ACGT")
                        for anchor, pattern in ((left, left_pattern), (right, right_pattern)))
    # Long, unique high-identity anchors also establish novel/zero alleles.
    # Weak fuzzy anchors still require the expected range; truncated flanks
    # always require the stronger check, including for in-range measurements.
    flanks_resolved = bool(
        left and right and flank_repeat_length is not None
        and flank_repeat_length >= 0
        and (not partial_flanks or strong_flanks)
        and (
            left.edit_distance == right.edit_distance == 0
            or strong_flanks
            or plausible_minimum <= flank_repeat_length <= plausible_maximum
        )
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
        # Match assembly PCR/MLVA_finder: primer indels do not change the
        # calibrated product size, which uses the configured primer lengths.
        product_size = calibrated_product_size_bp
        if product_size is None:
            product_size = (
                reverse.start - forward.end
                + len(locus_model.forward_primer) + len(locus_model.reverse_primer)
            )
        raw_count, _ = assembly_equivalent_product_allele(locus_model, product_size, round_tolerance)
        method = "assembly_product_length"
    elif flanks_resolved and unit_length:
        if reference_nonrepeat_bp is None:
            raw_count, method = repeat_length / unit_length, "flank_bounded_repeat_length"
        elif len(reference_nonrepeat_bp) == 1:
            raw_count, _ = assembly_equivalent_product_allele(
                locus_model, reference_nonrepeat_bp[0] + repeat_length, round_tolerance)
            method = "assembly_calibrated_flank_length"
        else:
            method = "unresolved_flank_calibration"
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
    boundary_ambiguous = any(anchor.ambiguous for anchor in (left, right) if anchor)
    if raw_count is not None:
        # Measurement support is local to the observed sequence. Catalog bounds
        # must not allocate an unlimited allele grid for a single observation.
        center = math.floor(raw_count * 2) / 2
        candidates = sorted({center + delta / 2 for delta in range(-4, 6) if center + delta / 2 >= 0} | {called})
        anchor_indels = sum(anchor.insertions + anchor.deletions for anchor in (forward, reverse, left, right) if anchor is not None)
        qscore = mean_qscore(quality[repeat_start:repeat_end]
                            if quality is not None and repeat_start is not None and repeat_end is not None and repeat_end > repeat_start
                            else quality)
        sigma = max(0.08, (0.35 + 0.12 * anchor_indels + max(0.0, 15.0 - qscore) / 30.0) / max(math.sqrt(unit_length or 1), 1.0))
        if quality is not None:
            start, end = (forward.end, reverse.start) if forward and reverse else (repeat_start, repeat_end)
            sizing_quality = mean_qscore(quality[start:end])
            # Longer/noisier products have more opportunities for sizing errors.
            sigma = max(sigma, math.sqrt(max(0, end - start) * 10 ** (-sizing_quality / 10)) / (unit_length or 1))
        if boundary_ambiguous:
            sigma = max(sigma, 1.0)
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
        status = "PRESENCE_ONLY"
        failure = ("reference nonrepeat length is unavailable or ambiguous" if method == "unresolved_flank_calibration"
                   else "insufficient anchors to measure the repeat interval")
    else:
        status, failure = "ANCHOR_FAILURE", "no complete locus anchors were identified"
    product_sequence = sequence[product_start:product_end] if product_start is not None and product_end is not None else ""
    product_quality = quality[product_start:product_end] if quality is not None and product_start is not None and product_end is not None else None
    return LocusMeasurement(locus_model.locus_id, sequence_id, source, product_start, product_end, forward.start if forward else None, forward.end if forward else None, reverse.start if reverse else None, reverse.end if reverse else None, repeat_start, repeat_end, repeat_length, unit_length, raw_count, called, forward.identity if forward else None, reverse.identity if reverse else None, motif_identity, likelihoods, confidence, status, failure, forward, reverse, credible, second_allele, second_probability, product_sequence, product_quality, repeat_sequence, method,
                            {"repeat_sigma": sigma,
                             "boundary_ambiguous": boundary_ambiguous,
                             "flanks_resolved": flanks_resolved} if raw_count is not None else {})
