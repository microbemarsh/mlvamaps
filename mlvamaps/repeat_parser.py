from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

from .calling import estimate_repeat_count_from_inner_length, repeat_unit_length
from .concurrency import DEFAULT_THREADS, resolve_threads
from .models import Assignment, Locus, RepeatFeature
from .progress import ProgressReporter
from .sequence import find_best, mean_qscore, repeat_motif_statistics


def _bounded_find(pattern: str, sequence: str, start: int, end: int, max_mismatches: int):
    if not pattern or start >= end:
        return None, None
    pos, mismatches = find_best(pattern, sequence[start:end], max_mismatches)
    if pos is None:
        return None, None
    return start + pos, mismatches


def extract_repeat_features(
    assignments: list[Assignment],
    loci: list[Locus],
    threads: int = DEFAULT_THREADS,
    progress: ProgressReporter | None = None,
) -> list[RepeatFeature]:
    by_locus = {locus.locus_id: locus for locus in loci}
    progress = progress or ProgressReporter(enabled=False)

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
        pattern_parts, mismatches, motif_kmers = repeat_motif_statistics(
            repeat_sequence, motif, motif_len
        )
        left_score = 1 - ((assignment.forward_mismatches or 0) / max(len(locus.forward_primer), 1))
        right_score = 1 - ((assignment.reverse_mismatches or 0) / max(len(locus.reverse_primer), 1))
        flank_scores = [left_score, right_score]
        if locus.left_flank_sequence:
            flank_scores.append(left_flank_score)
        if locus.right_flank_sequence:
            flank_scores.append(right_flank_score)
        flank_quality = sum(flank_scores) / len(flank_scores)
        amplicon_start = assignment.forward_start or 0
        amplicon_end = assignment.reverse_end if assignment.reverse_end is not None else len(sequence)
        amplicon_sequence = sequence[amplicon_start:amplicon_end]
        amplicon_quality = (
            assignment.oriented_quality[amplicon_start:amplicon_end]
            if assignment.oriented_quality is not None
            else None
        )
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
            mismatches,
            motif_kmers,
            round(left_score, 4),
            round(right_score, 4),
            round(left_flank_score, 4),
            round(right_flank_score, 4),
            amplicon_sequence,
            amplicon_quality,
            min(len(amplicon_sequence), max(0, repeat_start - amplicon_start)),
            min(len(amplicon_sequence), max(0, repeat_end - amplicon_start)),
        )

    thread_count = resolve_threads(threads)
    total = len(assignments)
    progress.step(f"Extracting repeat features from {total:,} assignments with {thread_count} worker(s)")
    if thread_count == 1 or len(assignments) <= 1:
        results = []
        for idx, assignment in enumerate(assignments, start=1):
            results.append(extract_one(assignment))
            progress.count("Parsed repeat features", idx, total)
    else:
        with ThreadPoolExecutor(max_workers=thread_count) as executor:
            results = []
            for idx, result in enumerate(executor.map(extract_one, assignments), start=1):
                results.append(result)
                progress.count("Parsed repeat features", idx, total)
    progress.count("Parsed repeat features", total, total, force=True)
    return [feature for feature in results if feature is not None]
