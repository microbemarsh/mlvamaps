from __future__ import annotations

from .models import Assignment, Locus, ReadRecord
from .progress import ProgressReporter
from .sequence import revcomp


def assignments_from_pcr(
    reads: list[ReadRecord],
    loci: list[Locus],
    rows: list[dict[str, str | int]],
    sample_id: str,
    min_assignment_score: float = 0.55,
    progress: ProgressReporter | None = None,
) -> list[Assignment]:
    """Build read assignments from native paired-primer products."""
    progress = progress or ProgressReporter(enabled=False)
    locus_by_id = {locus.locus_id: locus for locus in loci}
    read_by_id = {read.read_id: read for read in reads}
    candidates_by_read: dict[str, list] = {}
    for row in rows:
        read_id = str(row.get("reference_id", "")).split()[0]
        locus = locus_by_id.get(str(row.get("primer_name", "")))
        if locus is None:
            continue
        product_size = int(row["full_len"])
        min_len = locus.expected_amplicon_min_bp or len(locus.forward_primer) + len(locus.reverse_primer)
        max_len = locus.expected_amplicon_max_bp or 100000
        if not min_len <= product_size <= max_len:
            continue
        read = read_by_id.get(read_id)
        if read is None:
            continue
        orientation = "reverse" if row.get("strand") == "-" else "forward"
        sequence = revcomp(read.sequence) if orientation == "reverse" else read.sequence
        quality = read.quality[::-1] if orientation == "reverse" and read.quality else read.quality
        f_mm = int(row["fwd_mismatches"])
        r_mm = int(row["rev_mismatches"])
        primer_len = max(len(locus.forward_primer) + len(locus.reverse_primer), 1)
        score = 1.2 - ((f_mm + r_mm) / primer_len)
        candidates_by_read.setdefault(read_id, []).append(
            (
                score,
                orientation,
                sequence,
                quality,
                locus,
                int(row["fwd_start"]),
                int(row["fwd_end"]),
                int(row["rev_start"]),
                int(row["rev_end"]),
                f_mm,
                r_mm,
            )
        )

    results = []
    total = len(reads)
    for idx, read in enumerate(reads, start=1):
        results.append(
            _assignment_from_candidates(
                read,
                sample_id,
                candidates_by_read.get(read.read_id, []),
                min_assignment_score,
            )
        )
        progress.count("Assigned reads", idx, total)
    progress.count("Assigned reads", total, total, force=True)
    return results


# Source compatibility for the MLVAMaps 0.1 API.
assignments_from_amplirust = assignments_from_pcr


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

    score, orientation, sequence, quality, locus, f_pos, f_end, r_pos, r_end, f_mm, r_mm = max(
        candidates, key=lambda item: item[0]
    )
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
        f_end,
        r_pos,
        r_end,
        f_mm,
        r_mm,
    )
