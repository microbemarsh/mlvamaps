"""Bounded parallel recruitment, using the sample's existing CPU allocation.

The seed filter and classify_pair thresholds are unchanged. Process workers
avoid serializing Python motif/anchor bookkeeping behind the interpreter lock;
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
from .repeat_likelihood import MoleculeEvidence, RepeatTemplate, classify_pair, flank_insert_length, pair_has_anchor
from .sequence import revcomp


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
            matches = []
            matched_template = None
            while candidates:
                bit = candidates & -candidates
                candidates ^= bit
                template = self.templates[bit.bit_length() - 1]
                result.locus_tests += 1
                if self.audit_mode == "compact":
                    item = pair_has_anchor(pair, template)
                else:
                    item = classify_pair(pair, template)
                if item:
                    matches.append(template.locus.locus_id if self.audit_mode == "compact" else item)
                    matched_template = template
                    if self.audit_mode == "compact" and len(matches) == 2:
                        result.ambiguity_witnesses.append((pair.molecule_id, matches[0], matches[1],
                                                           "no" if candidates else "yes"))
                        result.skipped_locus_tests += candidates.bit_count()
                        break
            if len(matches) == 1:
                result.evidence.append(classify_pair(pair, matched_template)
                                       if self.audit_mode == "compact" else matches[0])
                length = flank_insert_length(pair, matched_template)
                if length:
                    result.insert_lengths.append(length)
            elif matches:
                result.ambiguous_pairs += 1
                if self.audit_mode == "compact":
                    continue
                for item in matches:
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
