"""Long-read extraction into the shared candidate evidence model."""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache

from .alignment_evidence import CandidateAlignment, CandidateEvidence
from .candidate_contexts import CandidateContext
from .locus_measurement import measure_locus_product
from .models import Locus
from .sequence import revcomp
from .short_read_evidence import _repeat_indel_delta


def extract_long_read_evidence(
    alignments: list[CandidateAlignment],
    contexts: list[CandidateContext],
    loci: list[Locus],
    technology: str,
    *, round_tolerance: float = 0.25,
) -> list[CandidateEvidence]:
    context_by_id = {context.candidate_id: context for context in contexts}
    locus_by_id = {locus.locus_id: locus for locus in loci}
    grouped: dict[tuple[str, str, int | float], list[CandidateAlignment]] = defaultdict(list)
    for alignment in alignments:
        grouped[(alignment.molecule_id, alignment.locus_id, alignment.repeat_count)].append(alignment)
    @lru_cache(maxsize=512)
    def measure(locus_id, sequence, quality):
        observations = (
            measure_locus_product(sequence, locus_by_id[locus_id], quality, source="long_read_molecule", round_tolerance=round_tolerance),
            measure_locus_product(revcomp(sequence), locus_by_id[locus_id], quality[::-1] if quality else None, source="long_read_molecule", round_tolerance=round_tolerance),
        )
        observation = max(observations, key=lambda item: (item.status == "FULL_PRODUCT", item.status == "REPEAT_INFORMATIVE", item.confidence or 0))
        if observation.called_allele is not None and observation.status in {"FULL_PRODUCT", "REPEAT_INFORMATIVE"}:
            return (observation.called_allele, observation.status, observation.raw_repeat_count,
                    (observation.metadata or {}).get("repeat_sigma"), observation.product_sequence)
        return None

    best_by_molecule = {}
    for row in alignments:
        key = (row.molecule_id, row.locus_id)
        if key not in best_by_molecule or row.alignment_score > best_by_molecule[key].alignment_score:
            best_by_molecule[key] = row
    measured = {key: measure(key[1], row.query_sequence, row.query_quality)
                for key, row in best_by_molecule.items()}
    evidence = []
    for (molecule, locus_id, repeat), rows in grouped.items():
        best = max(rows, key=lambda row: (row.alignment_score, row.alignment_identity))
        context = context_by_id[best.candidate_id]
        left = best.reference_start < context.repeat_start < best.reference_end
        right = best.reference_start < context.repeat_end < best.reference_end
        full = best.reference_start <= context.repeat_start and best.reference_end >= context.repeat_end
        direct = measured.get((molecule, locus_id))
        indel_delta = _repeat_indel_delta(best, context)
        evidence.append(CandidateEvidence(
            locus_id, repeat, molecule, best.alignment_score, best.alignment_identity,
            min(1.0, best.alignment_identity * best.query_coverage),
            bool(direct and direct[1] == "FULL_PRODUCT"), full or bool(direct), left, right,
            -abs(indel_delta), 0.0, technology, direct[0] if direct else None,
            best.candidate_id, best.reference_id, best.mapping_quality,
            best.query_coverage, best.reference_coverage, best.cigar, best.cs,
            "direct_product" if direct and direct[1] == "FULL_PRODUCT" else
            "complete_vntr_span" if full or direct else "repeat_boundary" if left or right else
            "repeat_indel" if indel_delta else "generic_locus_mapping",
            {"product_sequence": direct[4] if direct and direct[1] == "FULL_PRODUCT" else "",
             "repeat_indel_observed_delta": indel_delta,
             "spanning_repeat_indel_delta": indel_delta if left and right else None,
             "aligned_flank_bases": max(
                 0, min(best.reference_end, context.repeat_start) - best.reference_start,
                 best.reference_end - max(best.reference_start, context.repeat_end),
             )},
            raw_repeat_count=direct[2] if direct else None,
            repeat_sigma=direct[3] if direct else None,
        ))
    return evidence
