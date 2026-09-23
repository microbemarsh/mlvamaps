"""Independent microbial E/S/F/FRR repeat-length likelihood model.

Conceptual reference: Mousavi et al. (2019), doi:10.1093/nar/gkz501.
No GangSTR code or diploid genotype assumptions are used here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from collections import Counter
import math

import numpy as np
import parasail

from .models import Locus, ReadPair
from .sequence import revcomp
from .repeat_calibration import repeat_unit_length
from .mapping_classification import fit_reference_mixture

_MATRIX = parasail.matrix_create("ACGTN", 2, -4)


@dataclass(frozen=True)
class RepeatTemplate:
    locus: Locus
    left: str
    right: str
    motif: str
    unit: int

    def sequence(self, count: float) -> str:
        length = round(count * self.unit)
        return self.left + (self.motif * math.ceil(length / len(self.motif)))[:length] + self.right


def panel_template(locus: Locus) -> RepeatTemplate | None:
    motif = locus.repeat_motif.upper()
    if not (locus.left_flank_sequence and locus.right_flank_sequence and motif
            and not set(motif) - set("ACGT") and repeat_unit_length(locus)):
        return None
    return RepeatTemplate(locus, locus.forward_primer + locus.left_flank_sequence,
                          locus.right_flank_sequence + revcomp(locus.reverse_primer),
                          motif, repeat_unit_length(locus))


@dataclass(frozen=True)
class FlankHit:
    query_start: int
    query_end: int
    flank_start: int
    flank_end: int
    score: int
    edits: int


def flank_hit(sequence: str, flank: str, minimum: int = 10) -> FlankHit | None:
    if not sequence or not flank:
        return None
    result = parasail.sw_trace_striped_32(sequence, flank, 5, 1, _MATRIX)
    q, r = result.traceback.query, result.traceback.ref
    aligned = sum(a != "-" and b != "-" for a, b in zip(q, r))
    edits = sum(a != b for a, b in zip(q, r))
    if aligned < min(minimum, len(flank)) or edits / max(aligned, 1) > 0.15:
        return None
    return FlankHit(result.end_query + 1 - len(q.replace("-", "")), result.end_query + 1,
                    result.end_ref + 1 - len(r.replace("-", "")), result.end_ref + 1,
                    result.score, edits)


def repetitive(sequence: str, motif: str) -> bool:
    if not sequence or not motif:
        return False
    return max(sum(base == motif[(i + offset) % len(motif)] for i, base in enumerate(sequence))
               for offset in range(len(motif))) / len(sequence) >= 0.9


@dataclass
class MoleculeEvidence:
    molecule_id: str
    locus_id: str
    classes: tuple[str, ...]
    sequences: tuple[str, ...]
    orientations: tuple[str, ...]
    observed_repeat: float | None = None
    lower_bound: float = 0.0
    fragment_offset: float | None = None
    product_sequence: str = ""
    alignment: dict = field(default_factory=dict)
    qualities: tuple[str | None, ...] = ()


def classify_pair(pair: ReadPair, template: RepeatTemplate) -> MoleculeEvidence | None:
    """Anchor to non-repeat sequence; count each pair once, with multiple tags."""
    reads = [pair.read1] + ([pair.read2] if pair.read2 else [])
    oriented = []
    hits = []
    strands = []
    repeat_only = []
    qualities = []
    for read in reads:
        choices = []
        for strand, seq in (("+", read.sequence.upper()), ("-", revcomp(read.sequence))):
            is_repeat = repetitive(seq, template.motif)
            left = None if is_repeat else flank_hit(seq, template.left)
            right = None if is_repeat else flank_hit(seq, template.right)
            choices.append((sum(h.score for h in (left, right) if h), strand, seq, left, right, is_repeat))
        _, strand, seq, left, right, is_repeat = max(choices, key=lambda v: (v[0], v[5]))
        oriented.append(seq)
        hits.append((left, right))
        strands.append(strand)
        qualities.append(read.quality[::-1] if strand == "-" and read.quality else read.quality)
        repeat_only.append(is_repeat)
    if not any(left or right for left, right in hits):
        return None  # Unanchored FRRs cannot identify a microbial locus.
    classes = set()
    observed = []
    lower = 0.0
    product = ""
    offsets = []
    for seq, quality, (left, right) in zip(oriented, qualities, hits):
        start = left.query_end + len(template.left) - left.flank_end if left else None
        end = right.query_start - right.flank_start if right else None
        if start is not None and end is not None and 0 <= start <= end <= len(seq):
            classes.add("E")
            observed.append((end - start) / template.unit)
            # Only a complete observed primer-bounded product is a phased haplotype.
            begin = start - len(template.left)
            finish = end + len(template.right)
            if begin >= 0 and finish <= len(seq) and (quality is None or all(ord(q)-33 >= 20 for q in quality[begin:finish])):
                product = seq[begin:finish]
        else:
            segment = seq[start:] if start is not None and 0 <= start < len(seq) else (
                seq[:end] if end is not None and 0 < end <= len(seq) else "")
            if segment and repetitive(segment, template.motif):
                classes.add("F")
                lower = max(lower, len(segment) / template.unit)
        offsets.append((start, end))
    fragment_offset = None
    if len(reads) == 2:
        # Physical inward orientation is required, including swapped mate order.
        for i, j in ((0, 1), (1, 0)):
            left, right = hits[i][0], hits[j][1]
            if left and right:
                if strands[i] != "+" or strands[j] != "-":
                    continue
                fragment_offset = offsets[i][0] + len(oriented[j]) - offsets[j][1]
                if fragment_offset >= 0:
                    classes.add("S")
                    break
                fragment_offset = None
        if any(repeat_only):
            classes.add("FRR")
            lower = max(lower, max((len(seq) / template.unit for seq, rep in zip(oriented, repeat_only) if rep), default=0))
        if not classes and any(hits[0]) and any(hits[1]) and strands[0] == strands[1]:
            classes.add("discordant")
    if observed and max(observed) - min(observed) > 0.5:
        classes = {"discordant"}
        observed = []
        product = ""
    return MoleculeEvidence(pair.molecule_id, template.locus.locus_id, tuple(sorted(classes)),
                            tuple(oriented), tuple(strands),
                            float(np.median(observed)) if observed else None, lower,
                            fragment_offset, product,
                            {str(i): {side: vars(hit) if hit else None for side, hit in zip(("left", "right"), pair_hits)}
                             for i, pair_hits in enumerate(hits)}, tuple(qualities))


@dataclass(frozen=True)
class InsertDistribution:
    mean: float
    sd: float
    pairs_used: int
    source: str


def estimate_insert_distribution(lengths, mean=None, sd=None) -> InsertDistribution | None:
    if (mean is None) != (sd is None):
        raise ValueError("insert mean and SD must be supplied together")
    if mean is not None:
        if not all(math.isfinite(v) and v > 0 for v in (mean, sd)):
            raise ValueError("insert mean and SD must be finite and positive")
        return InsertDistribution(mean, sd, 0, "override")
    values = np.asarray([v for v in lengths if math.isfinite(v) and v > 0], dtype=float)
    if len(values) < 10:
        return None
    median = np.median(values)
    mad = float(np.median(abs(values - median)))
    values = values[abs(values - median) <= max(5.0, 4.5 * 1.4826 * mad)]
    if len(values) < 10:
        return None
    return InsertDistribution(float(np.mean(values)), max(1.0, float(np.std(values))), len(values), "robust_flank_pairs")


def flank_insert_length(pair: ReadPair, template: RepeatTemplate) -> float | None:
    """Estimate inserts only from two uniquely placed reads on the SAME flank.

    VNTR-spanning pairs must not estimate their own unknown repeat length.
    """
    if pair.read2 is None:
        return None
    for first, second in ((pair.read1.sequence, revcomp(pair.read2.sequence)),
                          (pair.read2.sequence, revcomp(pair.read1.sequence))):
        for flank in (template.left, template.right):
            a, b = flank.find(first), flank.find(second)
            if a >= 0 and b >= a and flank.find(first, a + 1) < 0 and flank.find(second, b + 1) < 0:
                return b + len(second) - a
    return None


def _scores(evidence: list[MoleculeEvidence], template: RepeatTemplate, states: np.ndarray,
            insert: InsertDistribution | None) -> np.ndarray:
    sequences = [template.sequence(float(r)) for r in states]
    matrix = np.zeros((len(evidence), len(states)))
    # Repeated reads are common in amplicon data; cache independent realignments.
    cache = {}
    for i, item in enumerate(evidence):
        if "discordant" in item.classes:
            continue
        if "E" in item.classes:
            matrix[i] -= 0.5 * ((states - item.observed_repeat) / 0.20) ** 2
        elif "S" in item.classes and insert:
            fragment = item.fragment_offset + states * template.unit
            matrix[i] -= 0.5 * ((fragment - insert.mean) / insert.sd) ** 2
        if item.lower_bound:
            matrix[i] -= 0.5 * (np.maximum(item.lower_bound - states, 0) / 0.25) ** 2
        if "E" in item.classes or "F" in item.classes:
            for read in item.sequences:
                if repetitive(read, template.motif):
                    continue
                if read not in cache:
                    scores = np.asarray([parasail.sw_striped_32(read, seq, 5, 1, _MATRIX).score for seq in sequences], dtype=float)
                    cache[read] = (scores - scores.max()) / 8.0
                matrix[i] += cache[read]
    return matrix


@dataclass
class RepeatInference:
    states: np.ndarray
    log_likelihoods: np.ndarray
    posterior: np.ndarray
    molecule_likelihoods: np.ndarray
    fractions: dict[float, float]
    best: float | None
    second: float | None
    interval: tuple[float, float] | None
    confidence: float
    identifiable: bool
    limit_reached: bool
    counts: dict[str, int]
    em_history: list[float]


def infer_repeat(evidence: list[MoleculeEvidence], template: RepeatTemplate,
                 insert: InsertDistribution | None = None, maximum: int = 100,
                 minimum_probability: float = 0.8) -> RepeatInference:
    if maximum * template.unit + len(template.left) + len(template.right) > 1000000:
        raise ValueError("candidate locus exceeds the one-million-base sequence safety limit")
    if maximum < 1:
        raise ValueError("repeat safety limit must be positive")
    observed = [e.observed_repeat for e in evidence if e.observed_repeat is not None]
    bounds = [e.lower_bound for e in evidence]
    upper = min(maximum, max(10, template.locus.expected_max_repeats + 2,
                            math.ceil(max(observed + bounds + [0])) + 3))
    if insert:
        upper = min(maximum, max(upper, math.ceil((insert.mean + 4 * insert.sd) / template.unit)))
    step = 0.5 if template.unit >= 2 else 1.0
    while True:
        candidate_count = int(upper / step) + 1
        candidate_bases = candidate_count * (len(template.left) + len(template.right) + upper * template.unit / 2)
        if candidate_bases > 50_000_000:
            raise ValueError("candidate sequences exceed the 50-million-base memory safety limit; reduce the candidate ceiling")
        states = np.arange(0, upper + step / 2, step)
        matrix = _scores(evidence, template, states, insert)
        joint = matrix.sum(axis=0)
        posterior = np.exp(np.maximum(joint - joint.max(), -700))
        posterior /= posterior.sum()
        if upper >= maximum or posterior[-min(4, len(states)):].sum() < 0.01:
            break
        upper = min(maximum, upper * 2)
    ranking = np.argsort(-posterior, kind="stable")
    winner = int(ranking[0])
    peaks = {winner}
    for value in observed:
        peaks.add(int(np.argmin(abs(states - value))))
    # Local maxima of the molecule likelihood density expose separated
    # spanning-pair components without restricting them to known alleles.
    if insert and evidence:
        density = np.exp(np.maximum(matrix - matrix.max(axis=1, keepdims=True), -700)).sum(axis=0)
        for j in range(1, len(states) - 1):
            if density[j] > density[j-1] and density[j] >= density[j+1] and density[j] >= max(2, .01*len(evidence)):
                peaks.add(j)
    columns = sorted(peaks)
    fractions, history = fit_reference_mixture(matrix[:, columns], np.ones(len(evidence))) if evidence else (np.ones(len(columns)), [])
    components = {float(states[j]): float(f) for j, f in zip(columns, fractions)}
    best = max(components, key=lambda r: (components[r], posterior[int(np.argmin(abs(states-r)))])) if evidence else None
    # Confidence for mixtures uses component-conditioned molecules, not the
    # (incorrect) probability that every molecule has the dominant allele.
    confidence = float(posterior[winner])
    if len(columns) > 1 and evidence:
        selected = matrix[:, columns]
        assignment = selected.argmax(axis=1)
        dominant = columns.index(int(np.argmin(abs(states - best))))
        support = matrix[assignment == dominant]
        if len(support):
            scores = support.sum(axis=0)
            probs = np.exp(np.maximum(scores - scores.max(), -700))
            confidence = float(probs[int(np.argmin(abs(states-best)))] / probs.sum())
    selected = []
    mass = 0.0
    for index in ranking:
        selected.append(float(states[index]))
        mass += posterior[index]
        if mass >= 0.95:
            break
    identifiable = bool(observed or (insert and any("S" in e.classes for e in evidence)))
    limit = upper == maximum and (posterior[-min(2, len(states)):].sum() >= 0.01 or max(observed + bounds + [0]) > maximum)
    return RepeatInference(states, joint, posterior, matrix, components, best,
                           float(states[ranking[1]]) if len(ranking) > 1 else None,
                           (min(selected), max(selected)), confidence,
                           identifiable and confidence >= minimum_probability and not limit,
                           limit, dict(Counter(c for e in evidence for c in e.classes)), history)
