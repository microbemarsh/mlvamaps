"""Bounded parallel recruitment, using the sample's existing CPU allocation.

The seed filter and classify_pair thresholds are unchanged. Competing loci are
resolved only when one has substantially stronger non-repeat anchor support.
Process workers avoid serializing Python motif/anchor bookkeeping behind the interpreter lock;
only retained evidence returns to the parent. No FASTQs are independently
rescanned by each worker, and recruitment finishes before locus fitting begins.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from itertools import chain, islice
import csv
import io
import json
import multiprocessing

from .concurrency import bounded_ordered_map, resolve_threads
from .repeat_likelihood import (
    MoleculeEvidence, RepeatTemplate, classify_pair, flank_insert_length,
    pair_anchor_score, read_anchor_bound,
)
from .sequence import revcomp


def anchor_score(item: MoleculeEvidence) -> int:
    """Sum the strongest non-repeat anchor per mate, avoiding flank overlap.

    Repeat sequence and evidence classes do not increase assignment strength.
    Taking only one flank per mate avoids counting the same query bases twice
    when local alignments to the two flanks overlap.
    """
    return sum(max((hit["score"] for hit in hits.values() if hit), default=0)
               for hits in item.alignment.values())


def resolve_locus_matches(matches: list[MoleculeEvidence]) -> MoleculeEvidence | None:
    """Retain a unique hit or a clear winner; near ties remain ambiguous.

    Require both a 12-point margin (two mismatch penalties) and a margin of
    at least 20% of the winning anchor score. These conservative assignment
    guards are not calibrated probabilities.
    """
    if len(matches) == 1:
        return matches[0]
    if not matches:
        return None
    ranked = sorted(((anchor_score(item), index) for index, item in enumerate(matches)),
                    reverse=True)
    best, index = ranked[0]
    if best - ranked[1][0] >= max(12, 0.2 * best):
        return matches[index]
    return None


@dataclass
class RecruitmentChunk:
    examined: int = 0
    locus_tests: int = 0
    evidence: list[MoleculeEvidence] = field(default_factory=list)
    ambiguous: list[dict] = field(default_factory=list)
    insert_lengths: list[float] = field(default_factory=list)
    audit_tsv: str = ""
    ambiguity_witnesses: list[tuple[str, str, str, str]] = field(default_factory=list)
    ambiguous_pairs: int = 0
    unmatched_pairs: int = 0
    skipped_locus_tests: int = 0


class ShortReadRecruiter:
    def __init__(self, templates: dict[str, RepeatTemplate | None], sample_id: str, audit_fields=None,
                 audit_mode="full"):
        if audit_mode not in {"compact", "full"}:
            raise ValueError("recruitment audit mode must be compact or full")
        self.templates = [template for _, template in sorted(templates.items()) if template is not None]
        # Duplicate contexts share cached anchors already; a bound pass cannot
        # separate them and only adds overhead to full ambiguity audits.
        self.use_bounds = len({(t.left, t.right) for t in self.templates}) > 3
        self.sample_id = sample_id
        self.audit_fields = audit_fields
        self.audit_mode = audit_mode
        self.all_candidates = (1 << len(self.templates)) - 1
        self.seeds: dict[str, int] = {}
        # Bit masks replace repeated set unions without lengthening seeds or
        # losing the short/error-bearing anchors accepted by the original code.
        for index, template in enumerate(self.templates):
            for flank in (template.left, template.right):
                for position in range(len(flank) - 3):
                    seed = flank[position:position + 4]
                    self.seeds[seed] = self.seeds.get(seed, 0) | (1 << index)

    def recruit(self, pairs) -> RecruitmentChunk:
        result = RecruitmentChunk()
        audit_buffer = io.StringIO(newline="") if self.audit_fields else None
        audit_writer = csv.DictWriter(audit_buffer, fieldnames=self.audit_fields,
                                     delimiter="\t", extrasaction="ignore") if audit_buffer is not None else None
        for pair in pairs:
            result.examined += 1
            candidates = 0
            reads = (pair.read1, pair.read2) if pair.read2 else (pair.read1,)
            for read in reads:
                for sequence in (read.sequence.upper(), revcomp(read.sequence)):
                    for position in range(len(sequence) - 3):
                        candidates |= self.seeds.get(sequence[position:position + 4], 0)
                        if candidates == self.all_candidates:
                            break
                    if candidates == self.all_candidates:
                        break
                if candidates == self.all_candidates:
                    break
            # Score-only SW bounds include both strands and do not reject
            # repeat-only or low-identity hits. Thus no accepted anchor can
            # exceed this bound, including classify_pair's strand tie rules.
            candidates_with_bounds = []
            use_bounds = self.use_bounds and candidates.bit_count() > 3
            while candidates:
                bit = candidates & -candidates
                candidates ^= bit
                index = bit.bit_length() - 1
                template = self.templates[index]
                bound = sum(read_anchor_bound(read.sequence, template.left, template.right)
                            for read in reads) if use_bounds else 0
                candidates_with_bounds.append((bound, index))
            if use_bounds:
                candidates_with_bounds.sort(key=lambda item: (-item[0], item[1]))
            matches = []
            best = runner_up = 0
            for position, (bound, index) in enumerate(candidates_with_bounds):
                if use_bounds and matches:
                    # The remaining bounds are descending. Skip only if none
                    # can change a clear winner. Full ambiguity audits still
                    # classify all candidates to retain every matching row.
                    clear_winner = best - max(runner_up, bound) >= max(12, 0.2 * best)
                    settled_ambiguity = (self.audit_mode == "compact" and runner_up > bound
                                         and best - runner_up < max(12, 0.2 * best))
                    if clear_winner or settled_ambiguity:
                        result.skipped_locus_tests += len(candidates_with_bounds) - position
                        break
                template = self.templates[index]
                result.locus_tests += 1
                score = pair_anchor_score(pair, template)
                if score:
                    matches.append((score, index))
                    if score >= best:
                        best, runner_up = score, best
                    elif score > runner_up:
                        runner_up = score
            # Audit ties and rows retain panel order regardless of bound order.
            matches.sort(key=lambda item: (-item[0], item[1]))
            if matches and (len(matches) == 1 or best - runner_up >= max(12, 0.2 * best)):
                template = self.templates[matches[0][1]]
                selected = classify_pair(pair, template)
                result.evidence.append(selected)
                length = flank_insert_length(pair, template)
                if length:
                    result.insert_lengths.append(length)
            elif matches:
                result.ambiguous_pairs += 1
                if self.audit_mode == "compact":
                    result.ambiguity_witnesses.append((pair.molecule_id,
                        self.templates[matches[0][1]].locus.locus_id,
                        self.templates[matches[1][1]].locus.locus_id, "yes"))
                    continue
                for _, index in sorted(matches, key=lambda item: item[1]):
                    item = classify_pair(pair, self.templates[index])
                    row = {
                        "sample_id": self.sample_id, "locus_id": item.locus_id,
                        "molecule_id": item.molecule_id, "classes": "ambiguous_locus",
                        "alignment": json.dumps(item.alignment),
                    }
                    if audit_writer is None:
                        result.ambiguous.append(row)
                    else:
                        audit_writer.writerow(row)
            else:
                result.unmatched_pairs += 1
        if audit_buffer is not None:
            result.audit_tsv = audit_buffer.getvalue()
        return result


_WORKER_RECRUITER = None


def _initialize_worker(templates, sample_id, audit_fields, audit_mode):
    global _WORKER_RECRUITER
    _WORKER_RECRUITER = ShortReadRecruiter(templates, sample_id, audit_fields, audit_mode)


def _recruit_chunk(pairs):
    return _WORKER_RECRUITER.recruit(pairs)


def recruit_short_reads(pairs, templates, sample_id, threads=1, chunk_size=256,
                        statistics=None, audit_fields=None, audit_mode="full"):
    """Yield ordered batches; at most two chunks per allocated CPU are queued."""
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    workers = resolve_threads(threads)
    iterator = iter(pairs)

    def chunks():
        while chunk := list(islice(iterator, chunk_size)):
            yield chunk

    batches = chunks()
    # Small inputs avoid process startup, and cannot launch idle workers.
    initial = list(islice(batches, workers))
    workers = min(workers, len(initial))
    if statistics is not None:
        statistics["workers"] = workers
        statistics["chunk_size"] = chunk_size
        statistics["audit_mode"] = audit_mode
    if workers <= 1:
        recruiter = ShortReadRecruiter(templates, sample_id, audit_fields, audit_mode)
        for batch in chain(initial, batches):
            yield recruiter.recruit(batch)
        return
    # Batch samples already run in threads. Explicit spawn is safe there and
    # does not inherit BLAS/Sassy state from a multithreaded parent via fork.
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn"),
                             initializer=_initialize_worker, initargs=(templates, sample_id, audit_fields, audit_mode)) as executor:
        yield from bounded_ordered_map(executor, _recruit_chunk, chain(initial, batches),
                                       max_pending=2 * workers)
