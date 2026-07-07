from __future__ import annotations

from .models import ReadRecord
from .sequence import mean_qscore


def filter_reads(
    reads: list[ReadRecord],
    min_length: int = 50,
    max_length: int = 100000,
    min_qscore: float = 0.0,
) -> tuple[list[ReadRecord], list[dict]]:
    kept = []
    failed_length = 0
    failed_qscore = 0
    for read in reads:
        length = len(read.sequence)
        qscore = mean_qscore(read.quality)
        if length < min_length or length > max_length:
            failed_length += 1
            continue
        if qscore < min_qscore:
            failed_qscore += 1
            continue
        kept.append(read)
    rows = [
        {"metric": "input_reads", "value": len(reads)},
        {"metric": "passed_reads", "value": len(kept)},
        {"metric": "failed_length", "value": failed_length},
        {"metric": "failed_qscore", "value": failed_qscore},
    ]
    return kept, rows
