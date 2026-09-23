"""Independent microbial E/S/F/FRR repeat-length likelihood model.

Conceptual reference: Mousavi et al. (2019), doi:10.1093/nar/gkz501.
No GangSTR code or diploid genotype assumptions are used here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from collections import Counter
from functools import lru_cache
import math

import numpy as np
import parasail

from .models import Locus, ReadPair
from .sequence import revcomp
from .repeat_calibration import repeat_unit_length
from .mapping_classification import fit_reference_mixture

_MATRIX = parasail.matrix_create("ACGTN", 2, -4)
_KNOWN_BASES = str.maketrans('', '', 'ACGTN')


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


@lru_cache(maxsize=256)
def _read_profile(sequence: str):
    return parasail.profile_create_32(sequence, _MATRIX)


@lru_cache(maxsize=256)
def _score_profile(sequence: str):
    return parasail.profile_create_16(sequence, _MATRIX)


@lru_cache(maxsize=256)
def _flank_profile(sequence: str):
    return parasail.profile_create_8(sequence, _MATRIX)


@lru_cache(maxsize=128)
def _minimum_anchor_score(required: int) -> int:
    return min(2*n - 6*(3*n//20) for n in range(required, required + 7))


@lru_cache(maxsize=4096)
def flank_hit(sequence: str, flank: str, minimum: int = 10) -> FlankHit | None:
    if not sequence or not flank:
        return None
    # Short-read anchors normally fit in SIMD byte lanes. Long/high-scoring
    # anchors explicitly fall back on saturation; never trust clipped results.
    result = parasail.sw_trace_striped_profile_8(_flank_profile(sequence), flank, 5, 1)
    if result.saturated:
        result = parasail.sw_trace_striped_profile_32(_read_profile(sequence), flank, 5, 1)
    required = min(minimum, len(flank))
    # An accepted alignment with n paired bases and e <= floor(.15*n) edits
    # scores at least 2*n - 6*e (a gap costs no more than a mismatch here).
    # Its minimum occurs in the first seven n values: stepping n by seven
    # raises this bound by at least two. Avoid decoding impossible traces.
    if result.score < _minimum_anchor_score(required):
        return None
    traceback = result.traceback
    q, r = traceback.query, traceback.ref
    aligned = sum(a != "-" and b != "-" for a, b in zip(q, r))
    edits = sum(a != b for a, b in zip(q, r))
    if aligned < min(minimum, len(flank)) or edits / max(aligned, 1) > 0.15:
        return None
    return FlankHit(result.end_query + 1 - len(q.replace("-", "")), result.end_query + 1,
                    result.end_ref + 1 - len(r.replace("-", "")), result.end_ref + 1,
                    result.score, edits)


@lru_cache(maxsize=128)
def _primitive_motif(motif: str) -> str:
    return motif[:(motif + motif).find(motif, 1)]


@lru_cache(maxsize=128)
def _motif_phases(motif: str, length: int) -> np.ndarray:
    """Cache byte templates, not Python per-base comparisons, for each phase."""
    # A database template can itself contain repeated units. Duplicate phases
    # are identical evidence and need not be tested more than once.
    motif = _primitive_motif(motif)
    bases = np.frombuffer(motif.encode("ascii"), dtype=np.uint8)
    positions = (np.arange(len(motif))[:, None] + np.arange(length)) % len(motif)
    phases = bases[positions]
    phases.flags.writeable = False
    return phases


@lru_cache(maxsize=128)
def _cyclic_search_text(motif: str, length: int):
    size = len(motif) + length - 1
    text = (motif * math.ceil(size / len(motif)))[:size]
    return text, np.frombuffer(text.encode('ascii'), dtype=np.uint8)


def _seeded_repetitive(sequence: str, motif: str) -> bool | None:
    """Exact cyclic Hamming test using necessary, native substring searches.

    At most floor(L/10) substitutions are allowed. Splitting the query into
    that many plus one disjoint pieces guarantees at least one exact piece
    in every accepted phase. Only phases proposed by exact pieces need the
    full native Hamming comparison. None requests the dense fallback for
    low-complexity inputs with too many seed hits.
    """
    length, period = len(sequence), len(motif)
    errors = length // 10
    text, cyclic = _cyclic_search_text(motif, length)
    bases = np.frombuffer(sequence.encode('ascii'), dtype=np.uint8)
    checked = set()
    seed_hits = 0
    for piece in range(errors + 1):
        start = piece * length // (errors + 1)
        end = (piece + 1) * length // (errors + 1)
        seed = sequence[start:end]
        stop = period + len(seed) - 1
        position = text.find(seed, 0, stop)
        while position >= 0:
            phase = (position - start) % period
            if phase not in checked:
                if np.count_nonzero(bases == cyclic[phase:phase + length]) >= length - errors:
                    return True
                checked.add(phase)
            seed_hits += 1
            if seed_hits >= 128 or len(checked) >= 64:
                return None
            position = text.find(seed, position + 1, stop)
    return False


@lru_cache(maxsize=4096)
def repetitive(sequence: str, motif: str) -> bool:
    if not sequence or not motif:
        return False
    # Same cyclic Hamming test and 90% cutoff as before; NumPy performs the
    # phase/base loops in native code. Bound cached arrays for long contexts.
    if sequence.isascii() and motif.isascii():
        motif = _primitive_motif(motif)
        if len(motif) >= 128:
            accepted = _seeded_repetitive(sequence, motif)
            if accepted is not None:
                return accepted
        bases = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
        if len(sequence) * len(motif) <= 250_000:
            matches = np.count_nonzero(_motif_phases(motif, len(sequence)) == bases, axis=1)
            return bool(matches.max() / len(sequence) >= 0.9)
        motif_bases = np.frombuffer(motif.encode("ascii"), dtype=np.uint8)
        positions = np.arange(len(sequence))
        for offset in range(0, len(motif), 64):
            phases = np.arange(offset, min(offset + 64, len(motif)))[:, None]
            expected = motif_bases[(positions + phases) % len(motif)]
            if np.count_nonzero(expected == bases, axis=1).max() / len(sequence) >= 0.9:
                return True
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


def pair_has_anchor(pair: ReadPair, template: RepeatTemplate) -> bool:
    """Whether classify_pair would retain this pair, without building evidence.

    Every accepted flank hit has a positive score, so an orientation with an
    anchor always outranks one without anchors. One hit therefore proves that
    classify_pair returns evidence, regardless of its eventual evidence tags.
    Repeat-only orientations remain excluded exactly as in classify_pair.
    """
    for read in (pair.read1, pair.read2):
        if read is None:
            continue
        for sequence in (read.sequence.upper(), revcomp(read.sequence)):
            if not repetitive(sequence, template.motif) and (
                    flank_hit(sequence, template.left) or flank_hit(sequence, template.right)):
                return True
    return False


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


def _repeat_alignment_scores(read, template, states, previous):
    """Exact local scores with redundant long periodic targets collapsed.

    With match=2 and gap extension=1, a positive local alignment of an L-base
    read spans fewer than 3L reference bases: at most L paired bases and fewer
    than 2L reference-only gaps. Once a repeat exceeds 3L plus one motif period,
    removing whole periods preserves all possible scoring windows, including
    both flank boundaries and every internal repeat phase. This optimization
    is score-only; it must not be used for tracebacks or reference coordinates.
    """
    motif = _primitive_motif(template.motif)
    period = len(motif)
    bound = 3 * len(read) + period
    lengths = [round(float(state) * template.unit) for state in states]
    equivalent = [n if n <= bound else bound + (n - bound) % period for n in lengths]
    known = dict(zip(equivalent, previous))
    scores = list(previous)
    profile = None
    exact_match_scoring = not read.translate(_KNOWN_BASES)
    for length in equivalent[len(previous):]:
        if length not in known:
            target = template.left + (motif * math.ceil(length / period))[:length] + template.right
            if read and exact_match_scoring and read in target:
                score = 2 * len(read)  # The maximum attainable local score.
            elif not read:
                score = 0
            else:
                if profile is None:
                    profile = _score_profile(read)
                result = parasail.sw_striped_profile_16(profile, target, 5, 1)
                if result.saturated:
                    result = parasail.sw_striped_profile_32(_read_profile(read), target, 5, 1)
                score = result.score
            known[length] = score
        scores.append(known[length])
    return np.asarray(scores, dtype=float)


def _scores(evidence: list[MoleculeEvidence], template: RepeatTemplate, states: np.ndarray,
            insert: InsertDistribution | None, alignment_cache: dict | None = None) -> np.ndarray:
    matrix = np.zeros((len(evidence), len(states)))
    # Repeated reads are common in amplicon data; cache independent realignments.
    cache = {}
    raw_cache = {} if alignment_cache is None else alignment_cache
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
                    previous = raw_cache.get(read, np.empty(0))
                    scores = _repeat_alignment_scores(read, template, states, previous)
                    raw_cache[read] = scores
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
    alignment_cache = {}
    while True:
        candidate_count = int(upper / step) + 1
        candidate_bases = candidate_count * (len(template.left) + len(template.right) + upper * template.unit / 2)
        if candidate_bases > 50_000_000:
            raise ValueError("candidate sequences exceed the 50-million-base memory safety limit; reduce the candidate ceiling")
        states = np.arange(0, upper + step / 2, step)
        matrix = _scores(evidence, template, states, insert, alignment_cache)
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
