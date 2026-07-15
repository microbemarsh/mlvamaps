from __future__ import annotations

from .models import ReadRecord
from .progress import ProgressReporter
from .sequence import mean_qscore


def filter_reads(
    reads: list[ReadRecord],
    min_length: int = 50,
    max_length: int = 100000,
    min_qscore: float = 0.0,
    progress: ProgressReporter | None = None,
) -> tuple[list[ReadRecord], list[dict]]:
    progress = progress or ProgressReporter(enabled=False)
    kept = []
    failed_length = 0
    failed_qscore = 0
    total = len(reads)
    for idx, read in enumerate(reads, start=1):
        length = len(read.sequence)
        qscore = mean_qscore(read.quality)
        if length < min_length or length > max_length:
            failed_length += 1
            progress.count("QC-filtered reads", idx, total)
            continue
        if qscore < min_qscore:
            failed_qscore += 1
            progress.count("QC-filtered reads", idx, total)
            continue
        kept.append(read)
        progress.count("QC-filtered reads", idx, total)
    progress.count("QC-filtered reads", total, total, force=True)
    rows = [
        {"metric": "input_reads", "value": len(reads)},
        {"metric": "passed_reads", "value": len(kept)},
        {"metric": "failed_length", "value": failed_length},
        {"metric": "failed_qscore", "value": failed_qscore},
    ]
    return kept, rows
