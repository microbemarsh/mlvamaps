"""Paired-read flank recruitment and bounded, evidence-only local assembly."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, replace
import re

import regex

from .candidate_contexts import CandidateContext
from .models import Locus, ReadPair
from .sequence import revcomp


@dataclass(frozen=True)
class RecruitedFragment:
    pair: ReadPair
    left: bool
    right: bool
    rescued: bool


def _seed_matcher(seeds):
    # Share literal prefixes explicitly: a flat regex alternation scans every
    # seed at every base. The trie still executes entirely in the native engine.
    if not seeds:
        return None
    trie = {}
    for seed in seeds:
        node = trie
        for base in seed:
            node = node.setdefault(base, {})

    def pattern(node):
        choices = [base + pattern(child) for base, child in sorted(node.items())]
        return "(?:" + "|".join(choices) + ")" if len(choices) > 1 else "".join(choices)

    return re.compile("(?=(" + pattern(trie) + "))")


class FlankRecruiter:
    """Scan with the installed native regex engine; retain only relevant pairs."""

    def __init__(self, contexts: list[CandidateContext], seed_length: int = 21):
        owners = defaultdict(set)
        flanks = set()
        repeats = set()
        for context in contexts:
            flanks.add((context.locus_id, "left", context.sequence[:context.repeat_start]))
            flanks.add((context.locus_id, "right", context.sequence[context.repeat_end:]))
            repeats.add(context.sequence[context.repeat_start:context.repeat_end])
            if context.repeat_motif and set(context.repeat_motif) <= set("ACGT"):
                repeats.add(context.repeat_motif * (2 + seed_length // len(context.repeat_motif)))
        for locus, side, sequence in flanks:
            for start in range(len(sequence) - seed_length + 1):
                seed = sequence[start:start + seed_length]
                if set(seed) <= set("ACGT") and len(set(seed)) >= 3:
                    owners[min(seed, revcomp(seed))].add((locus, side))
        self.seeds = {}
        for seed, locations in owners.items():
            if len({locus for locus, _side in locations}) == 1:
                self.seeds[seed] = self.seeds[revcomp(seed)] = locations
        # A flank-like seed found inside any configured repeat is not unique.
        matcher = _seed_matcher(self.seeds)
        if matcher:
            for sequence in repeats:
                for match in matcher.finditer(sequence):
                    seed = match.group(1)
                    self.seeds.pop(seed, None)
                    self.seeds.pop(revcomp(seed), None)
        self.matcher = _seed_matcher(self.seeds)
        self.fragments: dict[str, dict[str, RecruitedFragment]] = defaultdict(dict)
        self.ambiguous_fragments = 0

    def add(self, pair: ReadPair) -> None:
        if self.matcher is None or pair.read2 is None:
            return
        hits = [
            {location for match in self.matcher.finditer(read.sequence)
             for location in self.seeds[match.group(1)]}
            for read in (pair.read1, pair.read2)
        ]
        locations = hits[0] | hits[1]
        loci = {locus for locus, _side in locations}
        if len(loci) != 1:
            self.ambiguous_fragments += bool(loci)
            return
        locus = next(iter(loci))
        self.fragments[locus].setdefault(pair.molecule_id, RecruitedFragment(
            pair, (locus, "left") in locations, (locus, "right") in locations,
            not hits[0] or not hits[1],
        ))


def assemble_fragments(
    fragments: dict[str, RecruitedFragment], contexts: list[CandidateContext],
    minimum_fragments: int = 3,
) -> tuple[str, set[str], str]:
    """Return only a unique, covered, acyclic flank-to-flank path.

    Perfect tandem-repeat cycles cannot identify copy number. Neither candidate
    length nor graph depth is used to choose how often to traverse such a cycle.
    """
    if len(fragments) < max(3, minimum_fragments):
        return "", set(), "insufficient_coverage"
    # ponytail: cap this Python graph at 2,000 fragments / 50,000 edges. If
    # larger pools become useful, move graph construction to a native backend.
    if len(fragments) > 2000:
        return "", set(), "fragment_limit"
    anchors = {
        (context.sequence[:context.repeat_start], context.sequence[context.repeat_end:])
        for context in contexts
    }
    if not anchors or len(anchors) > 8:
        return "", set(), "ambiguous_flanks"
    sequences = {
        name: (fragment.pair.read1.sequence, fragment.pair.read2.sequence)
        for name, fragment in fragments.items() if fragment.pair.read2 is not None
    }
    if len(sequences) < max(3, minimum_fragments):
        return "", set(), "insufficient_coverage"
    shortest = min(len(read) for reads in sequences.values() for read in reads)
    k = min(41, shortest - 10, *(len(anchor) + 1 for pair in anchors for anchor in pair))
    if k < 21:
        return "", set(), "short_anchors"
    # Native overlapping extraction; a fragment contributes at most once per
    # edge, including when both mates cover it or contain repeated copies.
    kmers = regex.compile(f"(?=([ACGT]{{{k}}}))")
    counts = Counter()
    for reads in sequences.values():
        observed = {word for read in reads for word in kmers.findall(read)}
        counts.update(observed | {revcomp(word) for word in observed})
        if len(counts) > 50000:
            return "", set(), "graph_limit"
    edges = defaultdict(set)
    for word, count in counts.items():
        if count >= 2:
            edges[word[:-1]].add(word[1:])
    products = set()
    failure = "no_path"
    for left, right in sorted(anchors):
        node, end = left[:k - 1], right[-(k - 1):]
        path, visited = [node], set()
        while node != end:
            if node in visited:
                failure = "repeat_cycle"
                break
            visited.add(node)
            successors = edges.get(node, set())
            if len(successors) != 1:
                if successors:
                    failure = "competing_paths"
                break
            node = next(iter(successors))
            path.append(node[-1])
        else:
            sequence = "".join(path)
            if sequence.startswith(left) and sequence.endswith(right):
                products.add(sequence)
    if len(products) != 1:
        return "", set(), "competing_paths" if products else failure
    sequence = products.pop()
    support = {
        name for name, reads in sequences.items()
        if any(read in sequence or revcomp(read) in sequence for read in reads)
    }
    if len(support) < max(3, minimum_fragments):
        return "", set(), "insufficient_coverage"
    return sequence, support, "assembled"


def resolve_fragment_calls(calls, evidence, molecule_calls, alignments, contexts,
                           loci, fragments, sample_id, thresholds, round_tolerance):
    """Try the graph only without a well-supported direct sequence measurement."""
    from .alignment_evidence import CandidateEvidence
    from .allele_inference import infer_alleles
    from .locus_measurement import measure_locus_product

    contexts_by_locus = defaultdict(list)
    rows_by_locus = defaultdict(list)
    alignments_by_locus = defaultdict(list)
    for context in contexts:
        contexts_by_locus[context.locus_id].append(context)
    for row in evidence:
        rows_by_locus[row.locus_id].append(row)
    for alignment in alignments:
        alignments_by_locus[alignment.locus_id].append(alignment)
    calls_by_locus = {row["locus"]: row for row in calls}
    diagnostics, products = [], []
    for locus in loci:
        locus_id = locus.locus_id
        call = calls_by_locus[locus_id]
        pool = fragments.get(locus_id, {})
        rows = rows_by_locus[locus_id]
        accepted = set(pool) | {row.molecule_id for row in rows}
        mapped = [row for row in alignments_by_locus[locus_id] if row.molecule_id in accepted]
        measured = {row.molecule_id: row.measured_repeat_count for row in rows if row.measured_repeat_count is not None}
        distribution = Counter(measured.values())
        dominant = distribution.most_common(1)
        strong = bool(dominant and dominant[0][1] >= thresholds.minimum_molecules
                      and dominant[0][1] / len(measured) >= thresholds.minimum_probability
                      and len(distribution) == 1)
        method = "direct_spanning" if measured else "candidate_likelihood" if call["repeat_count"] != "" else "existing_fallback"
        fragment_methods = {row.evidence_tier for row in rows if row.candidate_id == "" and row.measured_repeat_count is not None}
        if fragment_methods == {"merged_fragment"}:
            method = "merged_fragment"
        softclipped = {row.molecule_id for row in mapped if "S" in row.cigar}
        mapped_mates = defaultdict(set)
        for row in mapped:
            mapped_mates[row.molecule_id].add(row.mate)
        unmapped = {name for name in pool if not {1, 2} <= mapped_mates[name]}
        diagnostic = {
            "sample_id": sample_id, "locus_id": locus_id,
            "recruited_fragment_count": len(set(pool) | {row.molecule_id for row in mapped}),
            "flank_recruited_fragment_count": len(pool),
            "mate_rescued_fragment_count": sum(item.rescued for item in pool.values()),
            "existing_mapper_only_count": len({row.molecule_id for row in mapped} - set(pool)),
            "spanning_fragment_count": len(measured),
            "left_anchor_count": sum(item.left for item in pool.values()),
            "right_anchor_count": sum(item.right for item in pool.values()),
            "softclipped_fragment_count": len(softclipped),
            "unmapped_mate_count": len(unmapped),
            "dominant_product_fragment_count": dominant[0][1] if dominant else 0,
            "direct_repeat_count": dominant[0][0] if dominant else "",
            "candidate_repeat_count": call["repeat_count"],
            "microassembly_attempted": "no", "microassembly_success": "no",
            "microassembly_sequence_length": "", "microassembly_reason": "strong_direct" if strong else "no_flank_fragments",
        }
        # Once direct measurements agree, candidate-relative partial mappings
        # cannot outvote them. Keep weaker fragments as presence evidence.
        if strong:
            rows = [row if row.molecule_id in measured else replace(
                row, full_repeat_span=False, left_boundary_span=False,
                right_boundary_span=False, repeat_indel_support=0, pair_geometry_support=0,
            ) for row in rows]
        elif pool:
            diagnostic["microassembly_attempted"] = "yes"
            sequence, support, reason = assemble_fragments(pool, contexts_by_locus[locus_id], thresholds.minimum_molecules)
            diagnostic["microassembly_reason"] = reason
            if sequence:
                measurement = measure_locus_product(sequence, locus, source="short_read_microassembly", round_tolerance=round_tolerance)
                if measurement.called_allele is not None and measurement.status == "FULL_PRODUCT":
                    # Preserve each input fragment's ID; assembled support is
                    # not reported as independent complete spanning reads.
                    rows = [row for row in rows if row.molecule_id not in support]
                    rows.extend(CandidateEvidence(
                        locus_id, measurement.called_allele, name,
                        locus_confidence=1.0, left_boundary_span=pool[name].left,
                        right_boundary_span=pool[name].right, technology="illumina",
                        measured_repeat_count=measurement.called_allele,
                        evidence_tier="microassembly",
                    ) for name in sorted(support))
                    rows = [row if row.measured_repeat_count is not None else replace(
                        row, full_repeat_span=False, left_boundary_span=False,
                        right_boundary_span=False, repeat_indel_support=0, pair_geometry_support=0,
                    ) for row in rows]
                    diagnostic.update(microassembly_success="yes", microassembly_sequence_length=len(sequence))
                    products.append((locus_id, sequence))
                    method = "microassembly"
                else:
                    diagnostic["microassembly_reason"] = "product_not_measurable"
            if diagnostic["microassembly_success"] == "no" and not measured:
                # Candidate-relative boundary coordinates do not resolve how
                # many copies occur in an unbridged repeat. Keep presence and
                # the original candidate estimate in diagnostics, not a capped
                # allele fabricated from this inconclusive graph.
                rows = [replace(row, full_repeat_span=False, left_boundary_span=False,
                                right_boundary_span=False, repeat_indel_support=0,
                                pair_geometry_support=0) for row in rows]
                method = "existing_fallback"
        rows_by_locus[locus_id] = rows
        diagnostic["call_method"] = method
        diagnostics.append(diagnostic)
    evidence = [row for locus in loci for row in rows_by_locus[locus.locus_id]]
    calls, molecule_calls = infer_alleles(evidence, loci, contexts, sample_id, "illumina", thresholds)
    for diagnostic, call in zip(diagnostics, calls):
        diagnostic["final_repeat_count"] = call["repeat_count"]
    return calls, evidence, molecule_calls, diagnostics, products
