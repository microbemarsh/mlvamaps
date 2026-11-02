"""Long-read extraction into the shared candidate evidence model."""

from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
import math
import re

import parasail

from .alignment_evidence import CandidateAlignment, CandidateEvidence
from .candidate_contexts import CandidateContext
from .locus_measurement import mapped_locus_window, measure_locus_product
from .models import Locus
from .short_read_evidence import _repeat_indel_delta
from .allele_inference import InferenceThresholds


_MOTIF_MATRIX = parasail.matrix_create("ACGTN", 2, -3)


def flank_error_model(alignments, contexts):
    """Regularize shared ONT flank errors; repeat indels never train the model."""
    known = {row.candidate_id: row for row in contexts}
    selected = {}
    for row in alignments:
        if row.alignment_identity < .9 or not row.cs or row.secondary or row.supplementary:
            continue
        previous = selected.get(row.molecule_id)
        if previous is None or (row.alignment_score, row.candidate_id) > (previous.alignment_score, previous.candidate_id):
            selected[row.molecule_id] = row
    bases = substitutions = insertions = deletions = 0
    for row in selected.values():
        context = known[row.candidate_id]
        position = row.reference_start
        for token in re.findall(r"=[A-Za-z]+|:\d+|\*[A-Za-z]{2}|[+-][A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2}", row.cs):
            operation = token[0]
            length = int(token[1:]) if operation == ":" else 1 if operation == "*" else len(token) - 1
            if operation == "~":
                position += int(token[3:-2])
                continue
            if operation == "+":
                if position < context.repeat_start or position > context.repeat_end:
                    insertions += length
                continue
            outside = max(0, min(position + length, context.repeat_start) - position)
            outside += max(0, position + length - max(position, context.repeat_end))
            bases += outside
            substitutions += outside if operation == "*" else 0
            deletions += outside if operation == "-" else 0
            position += length
    fitted = len(selected) >= 20 and bases >= 2000
    # ponytail: pooled flank rates are a conservative error proxy, not a
    # basecaller calibration. Context-specific correction needs truth data.
    return {"source": "pooled_flanks" if fitted else "conservative_prior", "molecules": len(selected),
            "flank_bases": bases, **{name: (errors + prior * 2000) / (bases + 2000) if fitted else prior
            for name, errors, prior in (("substitution", substitutions, .002),
                                        ("insertion", insertions, .003), ("deletion", deletions, .003))}}


def refine_repeat_sequence(sequence, motif, limits):
    """Bounded native banded scoring; scores are not calibrated probabilities."""
    if not motif or set(motif) - set("ACGT") or not sequence:
        return {}, "unsupported_motif"
    unit = len(motif)
    center = round(len(sequence) / unit)
    band = max(16, 3 * unit)
    if (len(sequence) > limits.maximum_sequence_length or
            len(sequence) * (2 * band + 1) > limits.maximum_alignment_cells):
        return {}, "alignment_resource_limit"
    scores = {}
    for radius in (2, 4):
        candidates = range(max(0, center - radius), center + radius + 1)
        if len(candidates) > limits.maximum_candidates:
            return scores, "candidate_limit"
        for count in candidates:
            if count in scores:
                continue
            if count * unit > limits.maximum_sequence_length:
                return scores, "alignment_resource_limit"
            scores[count] = (parasail.nw_banded(sequence, motif * count, 4, 1, band, _MOTIF_MATRIX).score
                             if count else -4 - max(0, len(sequence) - 1))
        best = max(scores.values())
        winners = [count for count, score in scores.items() if score == best]
        if all(count not in (candidates.start, candidates.stop - 1) or count == 0 for count in winners):
            return scores, ""
    return scores, "candidate_boundary"


def long_read_measurement_metadata(observation, locus, error_model, limits, refine=True):
    """Apply the same uncertainty rules to mapped and primer-rescued products."""
    metadata = dict(observation.metadata or {})
    repeat_sequence = getattr(observation, "repeat_sequence", "")
    identity = getattr(observation, "repeat_motif_identity", None)
    if error_model and repeat_sequence:
        unit = max(1, locus.repeat_unit_length_bp or len(locus.repeat_motif))
        sigma = math.sqrt(len(repeat_sequence) * max(.001, error_model["insertion"] + error_model["deletion"])) / unit
        metadata["repeat_sigma"] = max(metadata.get("repeat_sigma", .08), sigma)
        metadata["error_model"] = error_model
    if identity is not None and identity < .9 and repeat_sequence and metadata.get("flanks_resolved"):
        scores, fallback = (refine_repeat_sequence(repeat_sequence, locus.repeat_motif, limits)
                            if refine else ({}, "refinement_read_limit"))
        if scores or fallback:
            metadata.update(repeat_alignment_scores=scores, refinement_fallback=fallback, inference_tier=2)
        # A noisy motif fit cannot uniquely distinguish an indel error
        # from a true partial unit. Preserve the measured length and
        # widen its uncertainty; consensus may subsequently resolve it.
        if scores:
            metadata["repeat_sigma"] = max(metadata.get("repeat_sigma", .08), .5)
    if metadata.get("boundary_ambiguous"):
        metadata["repeat_sigma"] = max(metadata.get("repeat_sigma", .08), 1.0)
    return metadata


def extract_long_read_evidence(
    alignments: list[CandidateAlignment],
    contexts: list[CandidateContext],
    loci: list[Locus],
    technology: str,
    *, round_tolerance: float = 0.25,
    limits: InferenceThresholds | None = None,
) -> list[CandidateEvidence]:
    limits = limits or InferenceThresholds()
    error_model = flank_error_model(alignments, contexts) if technology in {"ont", "ont-hq"} else None
    context_by_id = {context.candidate_id: context for context in contexts}
    locus_by_id = {locus.locus_id: locus for locus in loci}
    grouped: dict[tuple[str, str, int | float], list[CandidateAlignment]] = defaultdict(list)
    for alignment in alignments:
        grouped[(alignment.molecule_id, alignment.locus_id, alignment.repeat_count)].append(alignment)
    @lru_cache(maxsize=512)
    def measure(locus_id, sequence, quality, refine=True):
        observation = measure_locus_product(sequence, locus_by_id[locus_id], quality,
                                           source="long_read_molecule", round_tolerance=round_tolerance)
        if observation.called_allele is not None and observation.status in {"FULL_PRODUCT", "REPEAT_INFORMATIVE"}:
            metadata = long_read_measurement_metadata(observation, locus_by_id[locus_id], error_model, limits, refine)
            return (observation.called_allele, observation.status, observation.raw_repeat_count,
                    metadata.get("repeat_sigma"), observation.product_sequence, metadata)
        return None

    best_by_molecule = {}
    for row in alignments:
        key = (row.molecule_id, row.locus_id)
        if key not in best_by_molecule or row.alignment_score > best_by_molecule[key].alignment_score:
            best_by_molecule[key] = row
    windows = {key: mapped_locus_window(row, context_by_id[row.candidate_id])
               for key, row in best_by_molecule.items()}
    by_locus = defaultdict(list)
    for key in best_by_molecule:
        by_locus[key[1]].append(key)
    refinement_keys = {key for keys in by_locus.values() for key in sorted(keys, key=lambda key: (
        -best_by_molecule[key].alignment_identity, -best_by_molecule[key].mapping_quality, key[0]))[:limits.maximum_refinement_reads]}
    measured = {key: measure(key[1], window[0], window[1], key in refinement_keys) for key, window in windows.items()}
    evidence = []
    for (molecule, locus_id, repeat), rows in grouped.items():
        best = max(rows, key=lambda row: (row.alignment_score, row.alignment_identity))
        context = context_by_id[best.candidate_id]
        from .mapped_reconstruction import aligned_flanks
        left_bases, right_bases = aligned_flanks(best, context)
        left = left_bases >= max(1, min(21, context.repeat_start)) and best.reference_end > context.repeat_start
        right = right_bases >= max(1, min(21, len(context.sequence) - context.repeat_end)) and best.reference_start < context.repeat_end
        full = left and right
        direct = measured.get((molecule, locus_id))
        window = windows[(molecule, locus_id)]
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
             "query_offset": best.query_offset,
             **(direct[5] if direct else {}),
             "locus_sequence": window[0] if window[3] else "",
             "locus_quality": window[1] if window[3] else None,
             "repeat_indel_observed_delta": indel_delta,
             "spanning_repeat_indel_delta": indel_delta if left and right else None,
             "aligned_flank_bases": max(
                 0, min(best.reference_end, context.repeat_start) - best.reference_start,
                 best.reference_end - max(best.reference_start, context.repeat_end),
             )},
            raw_repeat_count=direct[2] if direct else None,
            repeat_sigma=direct[3] if direct else None,
            orientation="reverse" if best.reverse else "forward",
            quality_flags=tuple(flag for flag, present in (
                ("boundary_ambiguous", direct and direct[5].get("boundary_ambiguous")),
                ("supplementary", best.supplementary),
                ("refinement_limited", direct and direct[5].get("refinement_fallback")),
            ) if present),
        ))
    return evidence
