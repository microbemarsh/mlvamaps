from __future__ import annotations

from dataclasses import replace

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
    retain_unassigned_sequence: bool = True,
    allow_multi_locus: bool = False,
) -> list[Assignment]:
    """Build read assignments from native paired-primer products."""
    progress = progress or ProgressReporter(enabled=False)
    locus_by_id = {locus.locus_id: locus for locus in loci}
    read_by_id = reads if hasattr(reads, "get") else {read.read_id: read for read in reads}
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
                product_size,
            )
        )

    results = []
    total = len(reads)
    for idx, read in enumerate(reads, start=1):
        candidates = candidates_by_read.get(read.read_id, [])
        groups = [candidates]
        if allow_multi_locus and candidates:
            from .mapped_reconstruction import query_interval_groups
            def interval(candidate):
                start, end = candidate[5], candidate[8]
                return (len(read.sequence) - end, len(read.sequence) - start) if candidate[1] == "reverse" else (start, end)
            groups = query_interval_groups(candidates, interval)
        for candidates in groups:
            assignment = _assignment_from_candidates(
                read,
                sample_id,
                candidates,
                min_assignment_score,
            )
            if allow_multi_locus and assignment.passes_assignment_qc:
                # Measure this product, not the first primer pair on a read
                # that can independently cover several loci.
                start, end = assignment.forward_start, assignment.reverse_end
                assignment = replace(assignment,
                    oriented_sequence=assignment.oriented_sequence[start:end],
                    oriented_quality=assignment.oriented_quality[start:end] if assignment.oriented_quality else None,
                    forward_start=0, forward_end=assignment.forward_end - start,
                    reverse_start=assignment.reverse_start - start, reverse_end=end - start)
            if not retain_unassigned_sequence and not assignment.passes_assignment_qc:
                assignment = replace(assignment, oriented_sequence="", oriented_quality=None)
            results.append(assignment)
        progress.count("Assigned reads", idx, total)
    progress.count("Assigned reads", total, total, force=True)
    return results


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

    (
        score,
        orientation,
        sequence,
        quality,
        locus,
        f_pos,
        f_end,
        r_pos,
        r_end,
        f_mm,
        r_mm,
        product_size,
    ) = max(
        candidates, key=lambda item: item[0]
    )
    primers_ordered = f_pos is not None and r_pos is not None and r_pos > f_pos
    unambiguous = len({
        (candidate[4].locus_id, candidate[11])
        for candidate in candidates if candidate[0] == score
    }) == 1
    return Assignment(
        read.read_id,
        sample_id,
        locus.locus_id,
        round(score, 4),
        orientation,
        f_pos is not None,
        r_pos is not None,
        score >= min_assignment_score and primers_ordered and unambiguous,
        sequence,
        quality,
        f_pos,
        f_end,
        r_pos,
        r_end,
        f_mm,
        r_mm,
        product_size,
    )
