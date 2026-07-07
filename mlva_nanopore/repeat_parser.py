from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

from .calling import estimate_repeat_count_from_inner_length, repeat_unit_length
from .concurrency import resolve_threads
from .models import Assignment, Locus, RepeatFeature
from .sequence import find_best, hamming_distance, mean_qscore


def _bounded_find(pattern: str, sequence: str, start: int, end: int, max_mismatches: int):
    if not pattern or start >= end:
        return None, None
    pos, mismatches = find_best(pattern, sequence[start:end], max_mismatches)
    if pos is None:
        return None, None
    return start + pos, mismatches


def extract_repeat_features(assignments: list[Assignment], loci: list[Locus], threads: int = 0) -> list[RepeatFeature]:
    by_locus = {locus.locus_id: locus for locus in loci}

    def extract_one(assignment: Assignment) -> RepeatFeature | None:
        if not assignment.passes_assignment_qc or assignment.assigned_locus not in by_locus:
            return None
        locus = by_locus[assignment.assigned_locus]
        sequence = assignment.oriented_sequence
        repeat_start = assignment.forward_end or 0
        repeat_end = assignment.reverse_start if assignment.reverse_start is not None else len(sequence)
        left_flank_score = 0.0
        right_flank_score = 0.0
        if locus.left_flank_sequence:
            pos, mm = _bounded_find(locus.left_flank_sequence, sequence, repeat_start, repeat_end, 3)
            if pos is not None:
                repeat_start = pos + len(locus.left_flank_sequence)
                left_flank_score = 1 - (mm or 0) / max(len(locus.left_flank_sequence), 1)
        if locus.right_flank_sequence:
            pos, mm = _bounded_find(locus.right_flank_sequence, sequence, repeat_start, repeat_end, 3)
            if pos is not None:
                repeat_end = pos
                right_flank_score = 1 - (mm or 0) / max(len(locus.right_flank_sequence), 1)
        if repeat_end < repeat_start:
            return None
        repeat_sequence = sequence[repeat_start:repeat_end]
        motif = locus.repeat_motif or "N"
        motif_len = repeat_unit_length(locus) or max(len(motif), 1)
        raw_count = estimate_repeat_count_from_inner_length(locus, len(repeat_sequence))
        if raw_count is None:
            raw_count = len(repeat_sequence) / motif_len
        nearest = round(raw_count)
        pattern_parts = []
        mismatches = 0
        motif_kmers = 0
        for idx in range(0, len(repeat_sequence), motif_len):
            chunk = repeat_sequence[idx : idx + motif_len]
            if len(chunk) < motif_len:
                pattern_parts.append(f"{chunk}:partial")
                continue
            distance = hamming_distance(chunk, motif)
            if distance == 0:
                pattern_parts.append(motif)
                motif_kmers += 1
            else:
                pattern_parts.append(chunk)
                mismatches += distance
        indels = len(repeat_sequence) % motif_len
        left_score = 1 - ((assignment.forward_mismatches or 0) / max(len(locus.forward_primer), 1))
        right_score = 1 - ((assignment.reverse_mismatches or 0) / max(len(locus.reverse_primer), 1))
        flank_scores = [left_score, right_score]
        if locus.left_flank_sequence:
            flank_scores.append(left_flank_score)
        if locus.right_flank_sequence:
            flank_scores.append(right_flank_score)
        flank_quality = sum(flank_scores) / len(flank_scores)
        return RepeatFeature(
            assignment.read_id,
            locus.locus_id,
            repeat_start,
            repeat_end,
            len(repeat_sequence),
            motif,
            round(raw_count, 3),
            nearest,
            round(flank_quality, 4),
            "-".join(pattern_parts),
            repeat_sequence,
            round(mean_qscore(assignment.oriented_quality), 3),
            indels,
            mismatches,
            motif_kmers,
            round(left_score, 4),
            round(right_score, 4),
            round(left_flank_score, 4),
            round(right_flank_score, 4),
        )

    thread_count = resolve_threads(threads)
    if thread_count == 1 or len(assignments) <= 1:
        results = [extract_one(assignment) for assignment in assignments]
    else:
        with ThreadPoolExecutor(max_workers=thread_count) as executor:
            results = list(executor.map(extract_one, assignments))
    return [feature for feature in results if feature is not None]
