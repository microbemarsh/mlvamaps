"""Illumina-specific extraction into the shared candidate evidence model."""

from __future__ import annotations

import math
import re
from collections import defaultdict
from functools import lru_cache

import regex

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
    *, round_tolerance: float = 0.25,
    recruited_fragments: dict | None = None,
) -> list[CandidateEvidence]:
    context_by_id = {context.candidate_id: context for context in contexts}
    locus_by_id = {locus.locus_id: locus for locus in loci}
    span_patterns = {}
    if recruited_fragments is not None:
        for locus in loci:
            anchors = (locus.forward_primer, revcomp(locus.reverse_primer),
                       locus.left_flank_sequence, locus.right_flank_sequence)
            span_patterns[locus.locus_id] = [
                regex.compile(f"(?:{anchor}){{e<={min(3, max(1, len(anchor) // 5))}}}")
                if anchor and set(anchor) <= set("ACGT") else None
                for anchor in anchors
            ]
        assigned = {name: locus_id for locus_id, pool in recruited_fragments.items() for name in pool}
        # A repeat-only competitor cannot establish locus identity. For pairs
        # missed by exact seeds, retain the mapper's aligned unique-flank bases.
        alignments = [row for row in alignments if (
            assigned.get(row.molecule_id) == row.locus_id
            or (row.molecule_id not in assigned and row.alignment_identity >= 0.9 and max(
                min(row.reference_end, context_by_id[row.candidate_id].repeat_start) - row.reference_start,
                row.reference_end - max(row.reference_start, context_by_id[row.candidate_id].repeat_end),
            ) >= 21)
        )]

    @lru_cache(maxsize=512)
    def measure_sequence(locus_id: str, sequence: str, quality: str | None):
        # Cache only the compact outcome, including failed measurements.
        # The cache is bounded and local to this sample/rounding setting;
        # quality and locus remain part of the key, and molecules stay distinct.
        measurements = []
        for oriented, qualities in ((sequence, quality), (revcomp(sequence), quality[::-1] if quality else None)):
            patterns = span_patterns.get(locus_id)
            if patterns and all(pattern is not None for pattern in patterns) and not (
                (patterns[0].search(oriented) and patterns[1].search(oriented))
                or (patterns[2].search(oriented) and patterns[3].search(oriented))
            ):
                continue
            measurements.append(measure_locus_product(
                oriented, locus_by_id[locus_id], qualities, source="illumina_molecule", round_tolerance=round_tolerance,
            ))
        if not measurements:
            return None, "PRESENCE_ONLY", 0
        best = max(measurements, key=lambda item: (item.status == "FULL_PRODUCT", item.status == "REPEAT_INFORMATIVE", item.confidence or 0))
        return best.called_allele, best.status, best.confidence or 0

    grouped: dict[tuple[str, str, int | float], list[CandidateAlignment]] = defaultdict(list)
    for alignment in alignments:
        grouped[(alignment.molecule_id, alignment.locus_id, alignment.repeat_count)].append(alignment)
    evidence = []
    direct_by_molecule: dict[tuple[str, str], tuple[int | float, str]] = {}
    fragment_methods = {}
    for locus_id, fragments in (recruited_fragments or {}).items():
        for molecule, fragment in fragments.items():
            pair = fragment.pair
            reads = [pair.read1, pair.read2]
            measured = [measure_sequence(locus_id, read.sequence, read.quality) for read in reads if read]
            complete = [item for item in measured if item[0] is not None and item[1] in {"FULL_PRODUCT", "REPEAT_INFORMATIVE"}]
            method = "direct_spanning"
            if not complete:
                merged = merge_read_pair(pair, require_unique=True)
                if merged is not None and not _merge_overlap_is_repeat_only(pair, merged, locus_by_id[locus_id]):
                    item = measure_sequence(locus_id, merged.sequence, merged.quality)
                    if item[0] is not None and item[1] in {"FULL_PRODUCT", "REPEAT_INFORMATIVE"}:
                        complete.append(item)
                        method = "merged_fragment"
            if complete and len({item[0] for item in complete}) == 1:
                best = max(complete, key=lambda item: (item[1] == "FULL_PRODUCT", item[2]))
                direct_by_molecule[(molecule, locus_id)] = best[:2]
                fragment_methods[(molecule, locus_id)] = method
    for (molecule, locus_id, _repeat), rows in grouped.items():
        if molecule in (recruited_fragments or {}).get(locus_id, {}):
            continue
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
            merged = merge_read_pair(pair, require_unique=recruited_fragments is not None)
            if merged is not None and not _merge_overlap_is_repeat_only(pair, merged, locus_by_id[locus_id]):
                sequences.append(merged)
        sequences.extend(ReadRecord(row.read_id, row.query_sequence, row.query_quality) for row in rows)
        measurements = (measure_sequence(locus_id, read.sequence, read.quality) for read in sequences)
        measured = max(measurements, key=lambda item: (item[1] == "FULL_PRODUCT", item[1] == "REPEAT_INFORMATIVE", item[2]), default=None)
        if measured and measured[0] is not None and measured[1] in {"FULL_PRODUCT", "REPEAT_INFORMATIVE"}:
            direct_by_molecule[(molecule, locus_id)] = measured[:2]

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
    # Candidate-independent evidence also covers mates absent from every
    # competitive alignment. Shared inference groups by the original molecule
    # ID, so seed, mate, and mapping paths never multiply fragment support.
    for locus_id, fragments in (recruited_fragments or {}).items():
        for molecule, fragment in fragments.items():
            measured = direct_by_molecule.get((molecule, locus_id))
            evidence.append(CandidateEvidence(
                locus_id, measured[0] if measured else 0, molecule,
                locus_confidence=1.0,
                direct_product_measurement=bool(measured and measured[1] == "FULL_PRODUCT"),
                full_repeat_span=measured is not None,
                # Flank recruitment alone establishes presence, not repeat size.
                technology="illumina", measured_repeat_count=measured[0] if measured else None,
                evidence_tier=fragment_methods.get((molecule, locus_id), "flank_recruitment"),
                metadata={"left_anchor": fragment.left, "right_anchor": fragment.right,
                          "mate_rescue": fragment.rescued},
            ))
    return evidence
