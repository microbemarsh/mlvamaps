from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

from .concurrency import resolve_threads
from .models import Assignment, Locus, ReadRecord
from .progress import ProgressReporter
from . import sequence as sequence_module
from .sequence import find_best, revcomp


_WORKER_LOCI: list[Locus] = []
_WORKER_SAMPLE_ID = ""
_WORKER_MAX_PRIMER_MISMATCHES = 3
_WORKER_MIN_ASSIGNMENT_SCORE = 0.55


def _init_process_worker(
    loci: list[Locus],
    sample_id: str,
    max_primer_mismatches: int,
    min_assignment_score: float,
) -> None:
    global _WORKER_LOCI, _WORKER_SAMPLE_ID, _WORKER_MAX_PRIMER_MISMATCHES, _WORKER_MIN_ASSIGNMENT_SCORE
    _WORKER_LOCI = loci
    _WORKER_SAMPLE_ID = sample_id
    _WORKER_MAX_PRIMER_MISMATCHES = max_primer_mismatches
    _WORKER_MIN_ASSIGNMENT_SCORE = min_assignment_score


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
    progress: ProgressReporter | None = None,
) -> list[Assignment]:
    thread_count = resolve_threads(threads)
    progress = progress or ProgressReporter(enabled=False)

    total = len(reads)
    progress.step(f"Assigning {total:,} reads to {len(loci):,} loci with {thread_count} worker(s)")
    chunk_size = max(1, min(500, total // max(thread_count * 8, 1) or 1))

    if thread_count == 1 or len(reads) <= 1:
        results = []
        for idx, read in enumerate(reads, start=1):
            results.append(_assign_read(read, loci, sample_id, max_primer_mismatches, min_assignment_score))
            progress.count("Assigned reads", idx, total)
        progress.count("Assigned reads", total, total, force=True)
        return results

    executor_cls = ThreadPoolExecutor if (sequence_module.sassy or sequence_module.edlib) else ProcessPoolExecutor
    if executor_cls is ThreadPoolExecutor:
        worker_args = [(read, loci, sample_id, max_primer_mismatches, min_assignment_score) for read in reads]
        results = []
        with executor_cls(max_workers=thread_count) as executor:
            for idx, assignment in enumerate(executor.map(_assign_one, worker_args, chunksize=chunk_size), start=1):
                results.append(assignment)
                progress.count("Assigned reads", idx, total)
    else:
        progress.step("Using process workers for pure-Python primer matching")
        results = []
        with executor_cls(
            max_workers=thread_count,
            initializer=_init_process_worker,
            initargs=(loci, sample_id, max_primer_mismatches, min_assignment_score),
        ) as executor:
            for idx, assignment in enumerate(executor.map(_assign_read_from_process, reads, chunksize=chunk_size), start=1):
                results.append(assignment)
                progress.count("Assigned reads", idx, total)
    progress.count("Assigned reads", total, total, force=True)
    return results


def _assign_one(args: tuple[ReadRecord, list[Locus], str, int, float]) -> Assignment:
    read, loci, sample_id, max_primer_mismatches, min_assignment_score = args
    return _assign_read(read, loci, sample_id, max_primer_mismatches, min_assignment_score)


def _assign_read_from_process(read: ReadRecord) -> Assignment:
    return _assign_read(
        read,
        _WORKER_LOCI,
        _WORKER_SAMPLE_ID,
        _WORKER_MAX_PRIMER_MISMATCHES,
        _WORKER_MIN_ASSIGNMENT_SCORE,
    )


def _assign_read(
    read: ReadRecord,
    loci: list[Locus],
    sample_id: str,
    max_primer_mismatches: int,
    min_assignment_score: float,
) -> Assignment:
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

    score, orientation, sequence, quality, locus, f_pos, r_pos, f_mm, r_mm = max(candidates, key=lambda item: item[0])
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
