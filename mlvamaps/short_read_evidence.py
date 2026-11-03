"""Illumina-specific extraction into the shared candidate evidence model."""

from __future__ import annotations

import math
import re
import statistics
from collections import defaultdict
from functools import lru_cache

import regex

from .alignment_evidence import CandidateAlignment, CandidateEvidence
from .candidate_contexts import CandidateContext
from .calling import assembly_equivalent_product_allele
from .locus_measurement import mapped_locus_window, measure_locus_product, reference_nonrepeat_lengths, reference_product_allele
from .models import Locus, ReadPair, ReadRecord
from .sequence import revcomp
from .short_reads import _merge_overlap_is_repeat_only, merge_read_pair


_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def estimate_insert_distribution(alignments, contexts, minimum_fragments=100):
    """Calibrate only uniquely mapped inward pairs wholly outside the VNTR."""
    known = {row.candidate_id: row for row in contexts}
    grouped = defaultdict(list)
    for row in alignments:
        grouped[row.molecule_id].append(row)
    lengths = []
    margins = []
    for rows in grouped.values():
        eligible = [row for row in rows if row.primary and not row.supplementary and row.mate in (1, 2)
                    and row.mapping_quality >= 20 and row.alignment_identity >= .99 and row.query_coverage >= .99]
        if len(eligible) != 2 or {row.mate for row in eligible} != {1, 2}:
            continue
        first, second = sorted(eligible, key=lambda row: row.reference_start)
        if first.candidate_id != second.candidate_id or first.reverse or not second.reverse:
            continue
        context = known[first.candidate_id]
        if not (second.reference_end <= context.repeat_start or first.reference_start >= context.repeat_end):
            continue
        lengths.append(second.reference_end - first.reference_start)
        margins.append(min(first.reference_start, len(context.sequence) - second.reference_end))
    if len(lengths) < minimum_fragments:
        return None
    median = statistics.median(lengths)
    mad = statistics.median(abs(length - median) for length in lengths)
    retained = [length for length in lengths if abs(length - median) <= max(1, 4.5 * mad)]
    if len(retained) < minimum_fragments:
        return None
    deviation = statistics.stdev(retained)
    # Catalog ends censor longer fragments. Do not fit a narrow library from
    # short targets, or infer a distribution from a single observed size.
    if deviation <= 0 or min(margins) < 3 * deviation:
        return None
    return statistics.fmean(retained), deviation, len(retained)


def partial_repeat_lower_bound(sequence, locus, reference_nonrepeat_bp=None):
    """Exact uniquely anchored terminal motif runs provide a conservative bound."""
    motif, flank = locus.repeat_motif, locus.left_flank_sequence
    if not motif or len(flank) < 12 or set(motif + flank) - set("ACGT"):
        return None
    flank = flank[-21:]
    position = sequence.find(flank)
    if position < 0 or sequence.find(flank, position + 1) >= 0:
        return None
    suffix = sequence[position + len(flank):]
    if len(suffix) < len(motif) or any(base != motif[i % len(motif)] for i, base in enumerate(suffix)):
        return None
    unit = locus.repeat_unit_length_bp or len(motif)
    nonrepeat = (len(locus.forward_primer) + len(locus.left_flank_sequence)
                 + len(locus.right_flank_sequence) + len(locus.reverse_primer))
    if reference_nonrepeat_bp is not None:
        if len(reference_nonrepeat_bp) != 1:
            return None
        nonrepeat = reference_nonrepeat_bp[0]
    if not (locus.expected_product_size_bp and locus.nominal_repeat_units):
        return assembly_equivalent_product_allele(locus, nonrepeat + len(suffix), 0)[0]
    offset = locus.nominal_repeat_units - (locus.expected_product_size_bp - nonrepeat) / unit
    return max(0, len(suffix) / unit + offset)


def _repeat_indel_delta(alignment: CandidateAlignment, context: CandidateContext) -> float:
    reference = alignment.reference_start
    delta = 0
    unit = max(context.repeat_unit_length, 1)
    for length_text, operation in _CIGAR.findall(alignment.cigar):
        length = int(length_text)
        if operation == "D":
            delta -= max(0, min(reference + length, context.repeat_end) - max(reference, context.repeat_start))
        elif operation == "I" and context.repeat_start <= reference <= context.repeat_end:
            delta += length
        if operation in "MDN=X":
            reference += length
    return delta / unit


def extract_short_read_evidence(
    alignments: list[CandidateAlignment],
    contexts: list[CandidateContext],
    loci: list[Locus],
    *, round_tolerance: float = 0.25,
    recruited_fragments: dict | None = None,
    insert_mean: float | None = None,
    insert_sd: float | None = None,
) -> list[CandidateEvidence]:
    if (insert_mean is None) != (insert_sd is None) or any(
        value is not None and (not math.isfinite(value) or value <= 0)
        for value in (insert_mean, insert_sd)
    ):
        raise ValueError("Insert-size mean and standard deviation must both be finite and positive")
    context_by_id = {context.candidate_id: context for context in contexts}
    locus_by_id = {locus.locus_id: locus for locus in loci}
    nonrepeat_lengths = reference_nonrepeat_lengths(loci, contexts)
    span_patterns = {}
    if recruited_fragments is not None:
        for locus in loci:
            calibrated = len(nonrepeat_lengths.get(locus.locus_id, ())) == 1
            anchors = (locus.forward_primer, revcomp(locus.reverse_primer),
                       locus.left_flank_sequence[-21:] if calibrated else locus.left_flank_sequence,
                       locus.right_flank_sequence[:21] if calibrated else locus.right_flank_sequence)
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
    def measure_sequence(locus_id: str, sequence: str, quality: str | None, oriented=False):
        # Cache only the compact outcome, including failed measurements.
        # The cache is bounded and local to this sample/rounding setting;
        # quality and locus remain part of the key, and molecules stay distinct.
        measurements = []
        orientations = [(sequence, quality)] if oriented else [(sequence, quality), (revcomp(sequence), quality[::-1] if quality else None)]
        for oriented, qualities in orientations:
            patterns = span_patterns.get(locus_id)
            if patterns and all(pattern is not None for pattern in patterns) and not (
                (patterns[0].search(oriented) and patterns[1].search(oriented))
                or (patterns[2].search(oriented) and patterns[3].search(oriented))
            ):
                continue
            measurements.append(measure_locus_product(
                oriented, locus_by_id[locus_id], qualities, source="illumina_molecule", round_tolerance=round_tolerance,
                reference_nonrepeat_bp=nonrepeat_lengths.get(locus_id),
            ))
        if not measurements:
            return None, "PRESENCE_ONLY", 0, None, None, "", {}
        best = max(measurements, key=lambda item: (item.status == "FULL_PRODUCT", item.status == "REPEAT_INFORMATIVE", item.confidence or 0))
        return (best.called_allele, best.status, best.confidence or 0,
                getattr(best, "raw_repeat_count", None),
                (getattr(best, "metadata", None) or {}).get("repeat_sigma"), best.product_sequence,
                {**(getattr(best, "metadata", None) or {}),
                 "product_quality": getattr(best, "product_quality", None) or ""})

    # Equal repeat labels can have different paired geometry or CIGAR evidence.
    # Retain each reference background; inference still groups by molecule ID.
    grouped: dict[tuple[str, str, str], list[CandidateAlignment]] = defaultdict(list)
    best_mates = {}
    for alignment in alignments:
        grouped[(alignment.molecule_id, alignment.locus_id, alignment.candidate_id)].append(alignment)
        key = (alignment.molecule_id, alignment.locus_id, alignment.mate or 1)
        if key not in best_mates or alignment.alignment_score > best_mates[key].alignment_score:
            best_mates[key] = alignment
    evidence = []
    direct_by_molecule: dict[tuple[str, str], tuple[int | float, str]] = {}
    conflicting_by_molecule = defaultdict(dict)

    def select_measurement(molecule, locus_id, measurements):
        usable = [item for item in measurements if item[0] is not None
                  and item[1] in {"FULL_PRODUCT", "REPEAT_INFORMATIVE"}]
        # A complete sample product determines its own calibration; partial
        # spans must not veto it based on unobserved nonrepeat sequence.
        usable = [item for item in usable if item[1] == "FULL_PRODUCT"] or usable
        if not usable:
            return None
        if len({(item[0], item[3]) for item in usable}) > 1:
            for item in usable:
                conflicting_by_molecule[(molecule, locus_id)][(item[0], item[3], item[4], item[5])] = item
            return None
        best = max(usable, key=lambda item: item[2])
        if len({item[5] for item in usable if item[1] == "FULL_PRODUCT"}) > 1:
            best = (*best[:5], "", best[6])
        return best

    guided_segments = {}
    fragment_methods = {}
    for locus_id, fragments in (recruited_fragments or {}).items():
        for molecule, fragment in fragments.items():
            pair = fragment.pair
            reads = [pair.read1, pair.read2]
            measured = []
            segments = []
            for mate, read in enumerate(reads, 1):
                if read is None:
                    continue
                mapped = best_mates.get((molecule, locus_id, mate))
                if mapped and mapped.query_sequence == read.sequence:
                    context = context_by_id[mapped.candidate_id]
                    # Partial mates need a local measurement window, not full
                    # product boundary realignment; the fragment graph joins them.
                    can_span = (mapped.query_start >= mapped.reference_start - 8
                                and len(read.sequence) - mapped.query_end >= len(context.sequence) - mapped.reference_end - 8)
                    sequence, quality, _, _ = mapped_locus_window(mapped, context, refine=False)
                    measurement = measure_sequence(locus_id, sequence, quality, True)
                    # Complete local primer products already verify the endpoints.
                    # Refine only unresolved products before consensus recovery.
                    if can_span and measurement[1] != "FULL_PRODUCT":
                        sequence, quality, _, spanning = mapped_locus_window(mapped, context)
                        if spanning:
                            segments.append((sequence, quality))
                        measurement = measure_sequence(locus_id, sequence, quality, True)
                    measured.append(measurement)
                else:
                    measured.append(measure_sequence(locus_id, read.sequence, read.quality))
            complete = [item for item in measured if item[0] is not None and item[1] in {"FULL_PRODUCT", "REPEAT_INFORMATIVE"}]
            if segments and len({item[0] for item in complete}) <= 1:
                sequence, quality = max(segments, key=lambda item: (
                    sum(map(ord, item[1] or "")) / max(1, len(item[1] or "")), len(item[0])))
                guided_segments[(molecule, locus_id)] = {"locus_sequence": sequence, "locus_quality": quality}
            method = "direct_spanning"
            if not complete:
                merged = merge_read_pair(pair, require_unique=True)
                if merged is not None and not _merge_overlap_is_repeat_only(pair, merged, locus_by_id[locus_id]):
                    item = measure_sequence(locus_id, merged.sequence, merged.quality)
                    if item[0] is not None and item[1] in {"FULL_PRODUCT", "REPEAT_INFORMATIVE"}:
                        complete.append(item)
                        method = "merged_fragment"
            best = select_measurement(molecule, locus_id, complete)
            if best is not None:
                direct_by_molecule[(molecule, locus_id)] = best
                fragment_methods[(molecule, locus_id)] = method
    observed_mates = defaultdict(list)
    for (molecule, locus_id, _mate), alignment in best_mates.items():
        observed_mates[(molecule, locus_id)].append(alignment)
    for (molecule, locus_id), rows in observed_mates.items():
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
        measurements = [measure_sequence(locus_id, read.sequence, read.quality) for read in sequences]
        measured = select_measurement(molecule, locus_id, measurements)
        if measured and measured[0] is not None and measured[1] in {"FULL_PRODUCT", "REPEAT_INFORMATIVE"}:
            direct_by_molecule[(molecule, locus_id)] = measured

    for (molecule, locus_id, _candidate), rows in grouped.items():
        best = max(rows, key=lambda row: (row.alignment_score, row.alignment_identity))
        context = context_by_id[best.candidate_id]
        repeat = context.repeat_count
        from .mapped_reconstruction import aligned_flanks
        anchor_sizes = [(row, aligned_flanks(row, context_by_id[row.candidate_id])) for row in rows]
        minimum_left = max(1, min(21, context.repeat_start))
        minimum_right = max(1, min(21, len(context.sequence) - context.repeat_end))
        left = any(sizes[0] >= minimum_left and row.reference_end > context.repeat_start for row, sizes in anchor_sizes)
        right = any(sizes[1] >= minimum_right and row.reference_start < context.repeat_end for row, sizes in anchor_sizes)
        full = any(sizes[0] >= minimum_left and sizes[1] >= minimum_right for row, sizes in anchor_sizes)
        indel_delta = _repeat_indel_delta(best, context)
        spanning = [row for row, sizes in anchor_sizes if row.candidate_id == best.candidate_id
                    and sizes[0] >= minimum_left and sizes[1] >= minimum_right]
        spanning_best = max(spanning, key=lambda row: row.alignment_score) if spanning else None
        measured = direct_by_molecule.get((molecule, locus_id))
        # Paired geometry is candidate-relative: deviation from the aligned
        # candidate product is represented in repeat units without using MAPQ.
        pair_rows = [row for row in rows if row.template_length]
        geometry = 0.0
        if pair_rows and any(row.reference_start < context.repeat_start for row in rows) and any(row.reference_end > context.repeat_end for row in rows):
            geometry = 1.0
        pair_estimate = {}
        if insert_mean is not None and measured is None and context.repeat_unit_length > 0:
            mates = {}
            for row in rows:
                if row.candidate_id == best.candidate_id and row.mate in (1, 2):
                    if row.mate not in mates or row.alignment_score > mates[row.mate].alignment_score:
                        mates[row.mate] = row
            if len(mates) == 2:
                first, second = sorted(mates.values(), key=lambda row: row.reference_start)
                # Only complete, inward-facing reads on opposite unique flanks
                # constrain the unsequenced interval. TLEN alone is insufficient.
                if (not first.reverse and second.reverse
                        and first.reference_end <= context.repeat_start
                        and second.reference_start >= context.repeat_end
                        and all(row.query_coverage >= .95 and row.alignment_identity >= .95
                                and row.reference_end - row.reference_start >= 21 for row in (first, second))):
                    reference_raw = reference_product_allele(context.sequence, locus_by_id[locus_id], 0)[0]
                    if reference_raw is not None:
                        pair_estimate = {
                            "pair_repeat_count": reference_raw + (insert_mean - (
                                second.reference_end - first.reference_start)) / context.repeat_unit_length,
                            "pair_repeat_sigma": insert_sd / context.repeat_unit_length,
                            "pair_alignment_score": first.alignment_score + second.alignment_score,
                        }
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
                "spanning_repeat_indel_delta": _repeat_indel_delta(spanning_best, context) if spanning_best else None,
                "aligned_flank_bases": max(
                    0, min(best.reference_end, context.repeat_start) - best.reference_start,
                    best.reference_end - max(best.reference_start, context.repeat_end),
                ),
                "product_sequence": measured[5] if measured and measured[1] == "FULL_PRODUCT" else "",
                **(measured[6] if measured else {}),
                **guided_segments.get((molecule, locus_id), {}),
                **pair_estimate,
            },
            raw_repeat_count=measured[3] if measured else None,
            repeat_sigma=measured[4] if measured else None,
            orientation=";".join(sorted({"reverse" if row.reverse else "forward" for row in rows})),
            quality_flags=tuple(flag for flag, present in (
                ("supplementary", any(row.supplementary for row in rows)),
                ("boundary_ambiguous", measured and measured[6].get("boundary_ambiguous")),
            ) if present),
        ))
    # Candidate-independent evidence also covers mates absent from every
    # competitive alignment. Shared inference groups by the original molecule
    # ID, so seed, mate, and mapping paths never multiply fragment support.
    for locus_id, fragments in (recruited_fragments or {}).items():
        for molecule, fragment in fragments.items():
            measured = direct_by_molecule.get((molecule, locus_id))
            bounds = [bound for read in (fragment.pair.read1, fragment.pair.read2) if read is not None
                      for sequence in (read.sequence, revcomp(read.sequence))
                      if (bound := partial_repeat_lower_bound(sequence, locus_by_id[locus_id],
                          nonrepeat_lengths.get(locus_id))) is not None]
            evidence.append(CandidateEvidence(
                locus_id, measured[0] if measured else None, molecule,
                locus_confidence=1.0,
                direct_product_measurement=bool(measured and measured[1] == "FULL_PRODUCT"),
                full_repeat_span=measured is not None,
                # Flank recruitment alone establishes presence, not repeat size.
                technology="illumina", measured_repeat_count=measured[0] if measured else None,
                evidence_tier=fragment_methods.get((molecule, locus_id), "flank_recruitment"),
                metadata={"left_anchor": fragment.left, "right_anchor": fragment.right,
                          "mate_rescue": fragment.rescued,
                          "repeat_lower_bound": max(bounds) if bounds and measured is None else None,
                          **(measured[6] if measured else {}),
                          **guided_segments.get((molecule, locus_id), {}),
                          "product_sequence": measured[5] if measured and measured[1] == "FULL_PRODUCT" else ""},
                raw_repeat_count=measured[3] if measured else None,
                repeat_sigma=measured[4] if measured else None,
                quality_flags=("boundary_ambiguous",) if measured and measured[6].get("boundary_ambiguous") else (),
            ))
    for (molecule, locus_id), alternatives in conflicting_by_molecule.items():
        for measured in alternatives.values():
            evidence.append(CandidateEvidence(
                locus_id, measured[0], molecule, technology="illumina", locus_confidence=1.0,
                measured_repeat_count=measured[0], raw_repeat_count=measured[3], repeat_sigma=measured[4],
                direct_product_measurement=measured[1] == "FULL_PRODUCT", full_repeat_span=True,
                evidence_tier="conflicting_measurements", quality_flags=("conflicting_measurements",),
                metadata={**measured[6], "product_sequence": measured[5]},
            ))
    return evidence
