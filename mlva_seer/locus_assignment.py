from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

from .concurrency import DEFAULT_THREADS, resolve_threads
from .models import Assignment, Locus, ReadRecord
from .progress import ProgressReporter
from . import sequence as sequence_module
from .sequence import find_best, find_best_many, revcomp


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
    return _score_primer_hits(locus, f_pos, r_pos, f_mm, r_mm)


def _score_primer_hits(locus: Locus, f_pos, r_pos, f_mm, r_mm):
    detected = int(f_pos is not None) + int(r_pos is not None)
    if detected == 0:
        return None
    reverse_len = len(locus.reverse_primer)
    primer_len = max(len(locus.forward_primer) + reverse_len, 1)
    mismatch_penalty = ((f_mm or 0) + (r_mm or 0)) / primer_len
    length_bonus = 0.0
    if f_pos is not None and r_pos is not None and r_pos > f_pos:
        amplicon_len = r_pos + reverse_len - f_pos
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
    threads: int = DEFAULT_THREADS,
    progress: ProgressReporter | None = None,
) -> list[Assignment]:
    thread_count = resolve_threads(threads)
    progress = progress or ProgressReporter(enabled=False)

    total = len(reads)
    progress.step(f"Assigning {total:,} reads to {len(loci):,} loci with {thread_count} worker(s)")
    chunk_size = max(1, min(500, total // max(thread_count * 8, 1) or 1))

    if sequence_module.sassy and len(reads) > 1:
        batched_results = _assign_reads_with_sassy(
            reads,
            loci,
            sample_id,
            max_primer_mismatches,
            min_assignment_score,
            thread_count,
            progress,
        )
        if batched_results is not None:
            return batched_results

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


def _assign_reads_with_sassy(
    reads: list[ReadRecord],
    loci: list[Locus],
    sample_id: str,
    max_primer_mismatches: int,
    min_assignment_score: float,
    thread_count: int,
    progress: ProgressReporter,
) -> list[Assignment] | None:
    patterns = []
    for locus in loci:
        patterns.extend((locus.forward_primer, revcomp(locus.reverse_primer)))
    progress.step(f"Searching all primer pairs with Sassy using {thread_count} Rust thread(s)")
    hits = find_best_many(
        patterns,
        [read.sequence for read in reads],
        max_primer_mismatches,
        thread_count,
    )
    if hits is None:
        progress.step("Installed Sassy lacks search_many; using per-read workers")
        return None

    candidate_keys: dict[int, set[tuple[int, int]]] = {}
    for pattern_idx, read_idx, orientation in hits:
        orientation_idx = 1 if orientation == "reverse" else 0
        candidate_keys.setdefault(read_idx, set()).add((orientation_idx, pattern_idx // 2))

    results = []
    total = len(reads)
    for read_idx, read in enumerate(reads):
        candidates = []
        for orientation_idx, locus_idx in sorted(candidate_keys.get(read_idx, ())):
            orientation = "reverse" if orientation_idx else "forward"
            sequence = revcomp(read.sequence) if orientation_idx else read.sequence
            quality = read.quality[::-1] if orientation_idx and read.quality else read.quality
            locus = loci[locus_idx]
            f_hit = hits.get((locus_idx * 2, read_idx, orientation))
            r_hit = hits.get((locus_idx * 2 + 1, read_idx, orientation))
            f_pos, f_mm = f_hit if f_hit is not None else (None, None)
            r_pos, r_mm = r_hit if r_hit is not None else (None, None)
            scored = _score_primer_hits(locus, f_pos, r_pos, f_mm, r_mm)
            if scored is None:
                continue
            score, f_pos, r_pos, f_mm, r_mm = scored
            candidates.append((score, orientation, sequence, quality, locus, f_pos, r_pos, f_mm, r_mm))
        results.append(_assignment_from_candidates(read, sample_id, candidates, min_assignment_score))
        progress.count("Assigned reads", read_idx + 1, total)
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
    return _assignment_from_candidates(read, sample_id, candidates, min_assignment_score)


def _assignment_from_candidates(
    read: ReadRecord,
    sample_id: str,
    candidates: list,
    min_assignment_score: float,
) -> Assignment:
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
