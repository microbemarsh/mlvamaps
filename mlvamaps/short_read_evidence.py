"""Illumina-specific extraction into the shared candidate evidence model."""

from __future__ import annotations

import math
import re
from collections import defaultdict

from .alignment_evidence import CandidateAlignment, CandidateEvidence
from .candidate_contexts import CandidateContext
from .locus_measurement import measure_locus_product
from .models import Locus, ReadPair, ReadRecord
from .sequence import revcomp
from .short_reads import _merge_overlap_is_repeat_only, merge_read_pair


_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def _repeat_indel_delta(alignment: CandidateAlignment, context: CandidateContext) -> float:
    reference = alignment.reference_start
    delta = 0
    unit = max(context.repeat_unit_length, 1)
    for length_text, operation in _CIGAR.findall(alignment.cigar):
        length = int(length_text)
        if operation in "M=X":
            reference += length
        elif operation == "D":
            if context.repeat_start - unit <= reference <= context.repeat_end + unit and length % unit == 0:
                delta -= length / unit
            reference += length
        elif operation == "I":
            if context.repeat_start - unit <= reference <= context.repeat_end + unit and length % unit == 0:
                delta += length / unit
    return delta


def extract_short_read_evidence(
    alignments: list[CandidateAlignment],
    contexts: list[CandidateContext],
    loci: list[Locus],
) -> list[CandidateEvidence]:
    context_by_id = {context.candidate_id: context for context in contexts}
    locus_by_id = {locus.locus_id: locus for locus in loci}
    grouped: dict[tuple[str, str, int | float], list[CandidateAlignment]] = defaultdict(list)
    for alignment in alignments:
        grouped[(alignment.molecule_id, alignment.locus_id, alignment.repeat_count)].append(alignment)
    evidence = []
    direct_by_molecule: dict[tuple[str, str], tuple[int | float, str]] = {}
    for (molecule, locus_id, _repeat), rows in grouped.items():
        if (molecule, locus_id) in direct_by_molecule:
            continue
        mates = {row.mate: row for row in rows if row.mate in (1, 2)}
        sequences: list[ReadRecord] = []
        if 1 in mates and 2 in mates:
            pair = ReadPair(
                molecule,
                ReadRecord(mates[1].read_id, mates[1].query_sequence, mates[1].query_quality),
                ReadRecord(mates[2].read_id, mates[2].query_sequence, mates[2].query_quality),
            )
            merged = merge_read_pair(pair)
            if merged is not None and not _merge_overlap_is_repeat_only(pair, merged, locus_by_id[locus_id]):
                sequences.append(merged)
        sequences.extend(ReadRecord(row.read_id, row.query_sequence, row.query_quality) for row in rows)
        measurements = []
        for sequence in sequences:
            measurements.extend((
                measure_locus_product(sequence.sequence, locus_by_id[locus_id], sequence.quality, source="illumina_molecule"),
                measure_locus_product(revcomp(sequence.sequence), locus_by_id[locus_id], sequence.quality[::-1] if sequence.quality else None, source="illumina_molecule"),
            ))
        measured = max(measurements, key=lambda item: (item.status == "FULL_PRODUCT", item.status == "REPEAT_INFORMATIVE", item.confidence or 0), default=None)
        if measured and measured.called_allele is not None and measured.status in {"FULL_PRODUCT", "REPEAT_INFORMATIVE"}:
            direct_by_molecule[(molecule, locus_id)] = (measured.called_allele, measured.status)

    for (molecule, locus_id, repeat), rows in grouped.items():
        best = max(rows, key=lambda row: (row.alignment_score, row.alignment_identity))
        scores_by_candidate: dict[str, float] = {}
        for row in rows:
            scores_by_candidate[row.candidate_id] = max(
                scores_by_candidate.get(row.candidate_id, -math.inf), row.alignment_score
            )
        alternatives = [
            score for candidate_id, score in scores_by_candidate.items()
            if candidate_id != best.candidate_id
        ]
        background_margin = (
            best.alignment_score - max(alternatives) if alternatives else math.inf
        )
        context = context_by_id[best.candidate_id]
        left = any(row.reference_start < context.repeat_start < row.reference_end for row in rows)
        right = any(row.reference_start < context.repeat_end < row.reference_end for row in rows)
        full = any(row.reference_start <= context.repeat_start and row.reference_end >= context.repeat_end for row in rows)
        indel_delta = _repeat_indel_delta(best, context)
        measured = direct_by_molecule.get((molecule, locus_id))
        # Paired geometry is candidate-relative: deviation from the aligned
        # candidate product is represented in repeat units without using MAPQ.
        pair_rows = [row for row in rows if row.template_length]
        geometry = 0.0
        if pair_rows and any(row.reference_start < context.repeat_start for row in rows) and any(row.reference_end > context.repeat_end for row in rows):
            geometry = 1.0
        evidence.append(CandidateEvidence(
            locus_id, repeat, molecule, best.alignment_score, best.alignment_identity,
            min(1.0, best.alignment_identity * best.query_coverage),
            bool(measured and measured[1] == "FULL_PRODUCT"),
            full or (measured is not None), left, right, -abs(indel_delta), geometry,
            "illumina", measured[0] if measured else None, best.candidate_id, best.reference_id,
            best.mapping_quality, best.query_coverage, best.reference_coverage,
            best.cigar, best.cs,
            "direct_product" if measured is not None else "complete_vntr_span" if full
            else "repeat_boundary" if left or right else "pair_geometry" if geometry
            else "generic_locus_mapping",
            {
                "repeat_indel_observed_delta": indel_delta,
                "background_alignment_margin": background_margin,
            },
        ))
    return evidence
