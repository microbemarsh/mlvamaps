"""Measure shared alleles without assigning indistinguishable loci separately."""

from collections import defaultdict
from dataclasses import replace
import hashlib
import json

from .allele_inference import InferenceThresholds, infer_alleles
from .mapped_reconstruction import recruit_mapping_fragments
from .short_read_evidence import extract_short_read_evidence


LOCUS_GROUP_FIELDS = [
    "sample", "group_id", "member_loci", "locus_assignment", "group_status",
    "repeat_count", "repeat_count_raw", "secondary_repeat", "allele_fractions",
    "candidate_distribution", "independent_support", "presence_support",
    "primary_molecule_support", "secondary_molecule_support", "full_span_support",
    "repeat_lower_bound", "model_support", "inference_method", "reason",
]


def _context_signature(context):
    return (context.sequence, context.repeat_start, context.repeat_end,
            context.repeat_unit_length, context.repeat_count,
            context.left_flank, context.right_flank, context.repeat_motif)


def infer_locus_groups(alignments, contexts, loci, ambiguous_assignments, sample_id,
                       thresholds=None, *, unmapped_mates=None, round_tolerance=.25,
                       insert_mean=None, insert_sd=None):
    """Reuse existing alignments to recover a group allele, never a locus identity.

    Identical catalog sequences alone are insufficient: every member must use
    identical boundaries, primer rules, and allele calibration. Only molecules
    already ambiguous entirely within a group participate. Group inference fits
    multiple alleles even for an isolate, because different member loci can
    carry different alleles. It never fills missing group counts from a prior.
    """
    thresholds = thresholds or InferenceThresholds()
    by_locus = defaultdict(list)
    for context in contexts:
        if context.observed_or_synthetic == "observed":
            by_locus[context.locus_id].append(context)
    groups = defaultdict(list)
    for locus in loci:
        if not by_locus[locus.locus_id]:
            continue
        calibration = tuple(getattr(locus, field) for field in (
            "forward_primer", "reverse_primer", "left_flank_sequence", "right_flank_sequence",
            "repeat_motif", "repeat_unit_length_bp", "expected_product_size_bp", "nominal_repeat_units",
            "expected_min_repeats", "expected_max_repeats",
            "expected_amplicon_min_bp", "expected_amplicon_max_bp",
        ))
        catalog = frozenset(_context_signature(context) for context in by_locus[locus.locus_id])
        groups[(calibration, catalog)].append(locus)
    by_molecule = defaultdict(list)
    for alignment in alignments:
        if alignment.molecule_id in ambiguous_assignments:
            by_molecule[alignment.molecule_id].append(alignment)
    rows = []
    for members in sorted(groups.values(), key=lambda group: sorted(locus.locus_id for locus in group)):
        if len(members) < 2:
            continue
        members = sorted(members, key=lambda locus: locus.locus_id)
        names = {locus.locus_id for locus in members}
        group_id = "group_" + hashlib.sha256(json.dumps(sorted(names)).encode()).hexdigest()[:16]
        group_locus = replace(members[0], locus_id=group_id)
        group_contexts, aliases = {}, {}
        for context in sorted((context for name in names for context in by_locus[name]),
                              key=lambda context: context.candidate_id):
            signature = _context_signature(context)
            canonical = group_contexts.setdefault(signature, replace(context, locus_id=group_id))
            aliases[context.candidate_id] = canonical.candidate_id
        eligible = {molecule for molecule, candidates in ambiguous_assignments.items()
                    if len(set(candidates)) >= 2 and set(candidates) <= names}
        selected = {}
        for alignment in (row for molecule in sorted(eligible) for row in by_molecule[molecule]):
            if alignment.candidate_id not in aliases:
                continue
            row = replace(alignment, locus_id=group_id, candidate_id=aliases[alignment.candidate_id])
            # Equivalent references and member loci can map the same mate to
            # the same coordinates. Preserve distinct mapping alternatives.
            key = (row.molecule_id, row.mate, row.candidate_id, row.reference_start, row.reference_end,
                   row.query_start, row.query_end, row.reverse, row.cigar, row.query_sequence, row.query_offset)
            previous = selected.get(key)
            if previous is None or (row.alignment_score, row.primary, row.mapping_quality) > (
                    previous.alignment_score, previous.primary, previous.mapping_quality):
                selected[key] = row
        group_alignments = sorted(selected.values(), key=lambda row: (
            row.molecule_id, row.mate or 0, row.candidate_id, row.reference_start, row.query_start,
            row.reverse, row.cigar, row.query_sequence))
        targets = list(group_contexts.values())
        fragments = recruit_mapping_fragments(group_alignments, targets, unmapped_mates=unmapped_mates)
        evidence = extract_short_read_evidence(group_alignments, targets, [group_locus],
            recruited_fragments=fragments, round_tolerance=round_tolerance,
            insert_mean=insert_mean, insert_sd=insert_sd)
        calls, _ = infer_alleles(evidence, [group_locus], targets, sample_id, "illumina", thresholds,
            require_measurement=True, sample_mode="metagenome", round_tolerance=round_tolerance)
        call = calls[0]
        status = {"called": "SHARED_ALLELE", "low_coverage": "LOW_DEPTH",
                  "mixed": "MULTIPLE_ALLELES", "estimated": "ESTIMATED",
                  "ambiguous": "AMBIGUOUS_ALLELE"}.get(call["status"], "UNRESOLVED")
        if call["repeat_count"] in ("", None) and call.get("repeat_lower_bound"):
            status = "LOWER_BOUND"
        allele_support = call.get("molecule_allele_counts", {})
        rows.append(dict(sample=sample_id, group_id=group_id, member_loci=";".join(sorted(names)),
            locus_assignment="AMBIGUOUS_LOCUS", group_status=status,
            **{field: call.get(field, "") for field in (
                "repeat_count", "repeat_count_raw", "secondary_repeat", "allele_fractions",
                "candidate_distribution", "presence_support", "full_span_support",
                "repeat_lower_bound", "inference_method")},
            independent_support=call.get("estimate_support", call["molecule_support"]),
            primary_molecule_support=allele_support.get(call["repeat_count"], 0),
            secondary_molecule_support=allele_support.get(call.get("secondary_repeat"), 0),
            model_support=call["best_probability"],
            reason="Shared sequence evidence cannot assign an allele or presence to individual member loci. "
                   "Group fractions combine loci and do not establish strain fractions. Model support is uncalibrated."))
    return rows
