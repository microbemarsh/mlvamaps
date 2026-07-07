from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

from .concurrency import resolve_threads
from .models import Assignment, Locus, ReadRecord
from .sequence import find_best, revcomp


def _score_locus(sequence: str, locus: Locus, max_primer_mismatches: int):
    reverse_site = revcomp(locus.reverse_primer)
    f_pos, f_mm = find_best(locus.forward_primer, sequence, max_primer_mismatches)
    r_pos, r_mm = find_best(reverse_site, sequence, max_primer_mismatches)
    detected = int(f_pos is not None) + int(r_pos is not None)
    if detected == 0:
        return None
    primer_len = max(len(locus.forward_primer) + len(reverse_site), 1)
    mismatch_penalty = ((f_mm or 0) + (r_mm or 0)) / primer_len
    length_bonus = 0.0
    if f_pos is not None and r_pos is not None and r_pos > f_pos:
        amplicon_len = r_pos + len(reverse_site) - f_pos
        if locus.expected_amplicon_min_bp <= amplicon_len <= locus.expected_amplicon_max_bp:
            length_bonus = 0.2
    score = detected / 2 - mismatch_penalty + length_bonus
    return score, f_pos, r_pos, f_mm, r_mm


def assign_reads(
    reads: list[ReadRecord],
    loci: list[Locus],
    sample_id: str,
    max_primer_mismatches: int = 3,
    min_assignment_score: float = 0.55,
    threads: int = 0,
) -> list[Assignment]:
    thread_count = resolve_threads(threads)

    def assign_one(read: ReadRecord) -> Assignment:
        candidates = []
        for orientation, sequence, quality in (
            ("forward", read.sequence, read.quality),
            ("reverse", revcomp(read.sequence), read.quality[::-1] if read.quality else None),
        ):
            for locus in loci:
                scored = _score_locus(sequence, locus, max_primer_mismatches)
                if scored is None:
                    continue
                score, f_pos, r_pos, f_mm, r_mm = scored
                candidates.append((score, orientation, sequence, quality, locus, f_pos, r_pos, f_mm, r_mm))
        if not candidates:
            return Assignment(
                read.read_id,
                sample_id,
                "UNASSIGNED",
                0.0,
                "unknown",
                False,
                False,
                False,
                read.sequence,
                read.quality,
                None,
                None,
                None,
                None,
                None,
                None,
            )
        else:
            score, orientation, sequence, quality, locus, f_pos, r_pos, f_mm, r_mm = max(
                candidates, key=lambda item: item[0]
            )
            reverse_len = len(locus.reverse_primer)
            primers_ordered = f_pos is not None and r_pos is not None and r_pos > f_pos
            return Assignment(
                read.read_id,
                sample_id,
                locus.locus_id,
                round(score, 4),
                orientation,
                f_pos is not None,
                r_pos is not None,
                score >= min_assignment_score and primers_ordered,
                sequence,
                quality,
                f_pos,
                f_pos + len(locus.forward_primer) if f_pos is not None else None,
                r_pos,
                r_pos + reverse_len if r_pos is not None else None,
                f_mm,
                r_mm,
            )

    if thread_count == 1 or len(reads) <= 1:
        return [assign_one(read) for read in reads]
    with ThreadPoolExecutor(max_workers=thread_count) as executor:
        return list(executor.map(assign_one, reads))
