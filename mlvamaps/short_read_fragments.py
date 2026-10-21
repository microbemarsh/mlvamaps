"""Short-read flank recruitment and bounded, evidence-only local assembly."""

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
    source: str = "flank_seed"


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
    """Scan with the installed native regex engine; retain relevant molecules."""

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
        if self.matcher is None:
            return
        hits = [
            {location for match in self.matcher.finditer(read.sequence)
             for location in self.seeds[match.group(1)]}
            for read in (pair.read1, pair.read2) if read is not None
        ]
        locations = set().union(*hits)
        loci = {locus for locus, _side in locations}
        if len(loci) != 1:
            self.ambiguous_fragments += bool(loci)
            return
        locus = next(iter(loci))
        self.fragments[locus].setdefault(pair.molecule_id, RecruitedFragment(
            pair, (locus, "left") in locations, (locus, "right") in locations,
            len(hits) == 2 and (not hits[0] or not hits[1]),
        ))


def assemble_fragments(
    fragments: dict[str, RecruitedFragment], contexts: list[CandidateContext],
    minimum_fragments: int = 3,
    mapped_alignments=None,
) -> tuple[str, set[str], str]:
    """Return only a unique, covered, acyclic flank-to-flank path.

    Increase word length when shorter overlaps collapse a repeat. Neither
    reference length nor depth chooses how often to traverse a repeat cycle.
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
    if not anchors or any(len(anchor) < 21 for pair in anchors for anchor in pair):
        return "", set(), "short_anchors"
    # Only the terminal seeds locate a path; the interior flanks are recovered
    # from reads, so sample SNPs/indels are not replaced by reference sequence.
    anchors = {(left[:21], right[-21:]) for left, right in anchors}
    # Mappings reaching the reference ends can locate sample-derived terminal
    # seeds, including variants absent from every database flank.
    lengths = {context.candidate_id: len(context.sequence) for context in contexts}
    observed_left, observed_right = defaultdict(set), defaultdict(set)
    for row in mapped_alignments or ():
        if row.molecule_id not in fragments or row.candidate_id not in lengths or row.alignment_identity < .9:
            continue
        sequence = revcomp(row.query_sequence) if row.reverse else row.query_sequence
        if row.query_end - row.query_start < 21:
            continue
        if row.reference_start == 0:
            observed_left[sequence[row.query_start:row.query_start + 21]].add(row.molecule_id)
        if row.reference_end == lengths[row.candidate_id]:
            observed_right[sequence[row.query_end - 21:row.query_end]].add(row.molecule_id)
    left_seeds = {seed for seed, names in observed_left.items() if len(seed) == 21 and len(names) >= minimum_fragments}
    right_seeds = {seed for seed, names in observed_right.items() if len(seed) == 21 and len(names) >= minimum_fragments}
    if left_seeds and right_seeds:
        if len(left_seeds) * len(right_seeds) > 8:
            return "", set(), "ambiguous_flanks"
        anchors = {(left, right) for left in left_seeds for right in right_seeds}
    if len(anchors) > 8:
        return "", set(), "ambiguous_flanks"
    sequences = {
        name: tuple(read.sequence for read in (fragment.pair.read1, fragment.pair.read2)
                    if read is not None)
        for name, fragment in fragments.items()
    }
    longest = max(len(read) for reads in sequences.values() for read in reads)
    # ponytail: at most seven graph sizes, capped at 151 bases. Longer-read
    # assembly belongs in a native assembler if this bounded fallback grows.
    maximum_k = min(151, longest - 10)
    if maximum_k < 22:
        return "", set(), "short_anchors"
    failure = "no_path"
    for k in sorted(set(range(min(41, maximum_k), maximum_k, 20)) | {maximum_k}):
        # Native overlapping extraction; a molecule contributes once per edge,
        # including when both mates cover it or contain repeated copies.
        kmers = regex.compile(f"(?=([ACGT]{{{k}}}))")
        counts = Counter()
        for reads in sequences.values():
            observed = {word for read in reads for word in kmers.findall(read)}
            counts.update(observed | {revcomp(word) for word in observed})
            if len(counts) > 50000:
                return "", set(), "graph_limit"
        edges = defaultdict(set)
        starts = defaultdict(set)
        left_seeds = {left for left, _right in anchors}
        for word, count in counts.items():
            if count >= 2:
                edges[word[:-1]].add(word[1:])
                if word[:21] in left_seeds:
                    starts[word[:21]].add(word[:-1])
        products = set()
        for left, right in sorted(anchors):
            for start in sorted(starts[left]):
                node, path, visited = start, [start], set()
                while not node.endswith(right):
                    if node in visited:
                        if failure != "competing_paths":
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
                    products.add("".join(path))
        if len(products) > 1:
            return "", set(), "competing_paths"
        if products:
            sequence = products.pop()
            support = {
                name for name, reads in sequences.items()
                if any(read in sequence or revcomp(read) in sequence for read in reads)
            }
            if len(support) >= max(3, minimum_fragments):
                return sequence, support, "assembled"
            failure = "insufficient_coverage"
    return "", set(), failure


def resolve_fragment_calls(calls, evidence, molecule_calls, alignments, contexts,
                           loci, fragments, sample_id, thresholds, round_tolerance,
                           *, require_measurement=False, sample_mode="metagenome"):
    """Try the graph only without a well-supported direct sequence measurement."""
    from .alignment_evidence import CandidateEvidence
    from .allele_inference import infer_alleles
    from .mapped_reconstruction import call_query_products

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
    diagnostics, products, pending = [], [], {}
    changed = not require_measurement
    for locus in loci:
        locus_id = locus.locus_id
        call = calls_by_locus[locus_id]
        pool = fragments.get(locus_id, {})
        rows = rows_by_locus[locus_id]
        accepted = set(pool) | {row.molecule_id for row in rows}
        mapped = [row for row in alignments_by_locus[locus_id] if row.molecule_id in accepted]
        measured = {row.molecule_id: row.measured_repeat_count for row in rows if row.measured_repeat_count is not None}
        complete = {row.molecule_id for row in rows
                    if row.direct_product_measurement and row.metadata.get("product_sequence")}
        distribution = Counter(measured.values())
        dominant = distribution.most_common(1)
        strong = bool(dominant and len(complete) >= thresholds.minimum_molecules
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
        unmapped = {name for name, fragment in pool.items()
                    if not mapped_mates[name] or (
                        fragment.pair.read2 is not None and not {1, 2} <= mapped_mates[name])}
        diagnostic = {
            "sample_id": sample_id, "locus_id": locus_id,
            "recruited_fragment_count": len(set(pool) | {row.molecule_id for row in mapped}),
            "flank_recruited_fragment_count": sum(item.source == "flank_seed" for item in pool.values()),
            "mapping_recruited_fragment_count": sum(item.source == "mapping" for item in pool.values()),
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
            sequence, support, reason = assemble_fragments(
                pool, contexts_by_locus[locus_id], thresholds.minimum_molecules, mapped,
            )
            diagnostic["microassembly_reason"] = reason
            if sequence:
                pending[locus_id] = (locus, sequence, support)
            # An inconclusive graph does not erase the existing candidate
            # estimate. Retain its likelihoods and label the approximation below.
        rows_by_locus[locus_id] = rows
        diagnostic["call_method"] = method
        diagnostics.append(diagnostic)
    # One native PCR batch for all reconstructed loci, rather than a subprocess
    # pair per locus. The same assembly rules select the measured product.
    pcr_results = call_query_products([(locus, sequence) for locus, sequence, support in pending.values()],
                                      sample_id, round_tolerance)
    for diagnostic in diagnostics:
        locus_id = diagnostic["locus_id"]
        if locus_id not in pcr_results:
            continue
        assembly_call, product = pcr_results[locus_id]
        if assembly_call["status"] != "PASS":
            diagnostic["microassembly_reason"] = "product_not_measurable"
            continue
        changed = True
        locus, sequence, support = pending[locus_id]
        pool = fragments[locus_id]
        allele = assembly_call["repeat_count"]
        # Preserve original molecule IDs: assembled support is not a new set
        # of independent complete spanning reads.
        rows = [row for row in rows_by_locus[locus_id] if row.molecule_id not in support]
        rows.extend(CandidateEvidence(
            locus_id, allele, name, locus_confidence=1.0,
            left_boundary_span=pool[name].left, right_boundary_span=pool[name].right,
            technology="illumina", measured_repeat_count=allele,
            evidence_tier="microassembly", raw_repeat_count=float(assembly_call["repeat_count_raw"]),
            repeat_sigma=0.08,
        ) for name in sorted(support))
        rows_by_locus[locus_id] = [row if row.measured_repeat_count is not None else replace(
            row, full_repeat_span=False, left_boundary_span=False, right_boundary_span=False,
            repeat_indel_support=0, pair_geometry_support=0,
        ) for row in rows]
        diagnostic.update(microassembly_success="yes", microassembly_sequence_length=len(sequence),
                          call_method="microassembly")
        products.append((locus_id, product))
    evidence = [row for locus in loci for row in rows_by_locus[locus.locus_id]]
    if changed:
        calls, molecule_calls = infer_alleles(
            evidence, loci, contexts, sample_id, "illumina", thresholds,
            require_measurement=require_measurement,
            sample_mode=sample_mode, round_tolerance=round_tolerance,
        )
    for diagnostic, call in zip(diagnostics, calls):
        if not require_measurement and call["best_candidate_repeat"] != "" and (
            call["status"] == "ambiguous"
            or (not diagnostic["spanning_fragment_count"] and diagnostic["microassembly_success"] == "no")
        ):
            reason = (
                "candidate likelihoods do not separate sufficiently"
                if call["status"] == "ambiguous" else "no independently resolved product"
            )
            if diagnostic["microassembly_attempted"] == "yes":
                reason += "; local assembly: " + str(diagnostic["microassembly_reason"]).replace("_", " ")
            candidates = [context.repeat_count for context in contexts_by_locus[call["locus"]]
                          if context.repeat_count is not None]
            if candidates and call["best_candidate_repeat"] in {min(candidates), max(candidates)}:
                reason += "; estimate is at the tested candidate range boundary"
            call.update(repeat_count=call["best_candidate_repeat"], status="estimated",
                        confidence="low", estimation_reason=reason)
            diagnostic["call_method"] = "candidate_likelihood"
        diagnostic["final_repeat_count"] = call["repeat_count"]
    return calls, evidence, molecule_calls, diagnostics, products
