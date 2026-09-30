"""Flank-anchored molecule evidence for targeted paired-end reconstruction."""
from __future__ import annotations

from dataclasses import dataclass, field
from collections import Counter
from functools import lru_cache
import math

import numpy as np
import parasail
import regex

from .models import Locus, ReadPair
from .sequence import revcomp
from .repeat_calibration import repeat_unit_length

_MATRIX = parasail.matrix_create("ACGTN", 2, -4)
_KNOWN_BASES = str.maketrans('', '', 'ACGTN')


@dataclass(frozen=True)
class RepeatTemplate:
    locus: Locus
    left: str
    right: str
    motif: str
    unit: int
    primer_only: bool = False

    def sequence(self, count: float) -> str:
        if self.primer_only:
            raise ValueError('Primer-only templates cannot synthesize unobserved locus sequence')
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


@lru_cache(maxsize=4096)
def primer_bounds(sequence: str, primer: str, start: int = 0):
    """Use native flank alignment for concrete primers; keep IUPAC support."""
    if not primer:
        return None
    position = sequence.find(primer, start)
    if position >= 0:
        return position, position+len(primer)
    if set(primer)-set('ACGT'):
        hit = _primer_hit(sequence[start:], primer)
        return (start + hit.query_start, start + hit.query_end) if hit else None
    maximum_edits = min(2, max(1, len(primer)//5))
    hit = flank_hit(sequence[start:], primer, max(1, len(primer)-maximum_edits))
    if not hit:
        return None
    missing_left, missing_right = hit.flank_start, len(primer)-hit.flank_end
    begin = start+hit.query_start-missing_left
    end = start+hit.query_end+missing_right
    if hit.edits+missing_left+missing_right <= maximum_edits and 0 <= begin < end <= len(sequence):
        return begin, end
    return None


@lru_cache(maxsize=4096)
def primer_product(sequence, template):
    """Recover physical primer bounds, independently of a reference allele."""
    forward = primer_bounds(sequence, template.locus.forward_primer)
    reverse = primer_bounds(sequence, revcomp(template.locus.reverse_primer),
                            forward[1] if forward else 0)
    if forward and reverse and reverse[0] >= forward[1]:
        return sequence[forward[0]:reverse[1]]
    return ''


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
def _flank_profile(sequence: str):
    return parasail.profile_create_8(sequence, _MATRIX)


@lru_cache(maxsize=128)
def _minimum_anchor_score(required: int) -> int:
    return min(2*n - 6*(3*n//20) for n in range(required, required + 7))


@lru_cache(maxsize=4096)
def flank_score_bound(sequence: str, flank: str) -> int:
    """Upper bound on an accepted flank score without decoding a traceback.

    Ignoring anchor length, identity and repeat-only rejection can only raise
    the score. Saturated byte lanes fall back to an unclipped wide alignment.
    """
    if not sequence or not flank:
        return 0
    result = parasail.sw_striped_profile_8(_flank_profile(sequence), flank, 5, 1)
    if result.saturated:
        result = parasail.sw_striped_profile_32(_read_profile(sequence), flank, 5, 1)
    return result.score


@lru_cache(maxsize=4096)
def read_anchor_bound(sequence: str, left_flank: str, right_flank: str) -> int:
    """Bound a mate's strongest anchor; share repeated mates across pairs."""
    return max(flank_score_bound(oriented, flank)
               for oriented in (sequence.upper(), revcomp(sequence))
               for flank in (left_flank, right_flank))


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
def cyclic_match(sequence: str, motif: str) -> bool:
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


@lru_cache(maxsize=4096)
@lru_cache(maxsize=4096)
def repeat_fraction(sequence: str, motif: str) -> float:
    """Matched bases in the best cyclic, gapped repeat alignment / read length.

    SIMD local alignment tolerates substitutions and small indels; both strands
    and all cyclic phases occur in the extended motif target. Only recruited
    molecules are tested, never all input reads against all panel motifs.
    """
    if not sequence or not motif:
        return 0.0
    sequence, motif = sequence.upper(), _primitive_motif(motif.upper())
    target = motif * math.ceil((len(sequence) * 1.2 + len(motif)) / len(motif))
    if sequence in target or revcomp(sequence) in target:
        return 1.0
    best = 0.0
    for oriented in (sequence, revcomp(sequence)):
        result = parasail.sw_trace_striped_profile_32(_read_profile(oriented), target, 5, 1)
        q, r = result.traceback.query, result.traceback.ref
        matches = sum(a == b and a in "ACGT" for a, b in zip(q, r))
        best = max(best, matches / len(sequence))
    return best


def repetitive(sequence: str, motif: str) -> bool:
    # Strict gate for excluding repeat-only anchors. Rich-read tagging below
    # uses a configurable lower fraction without suppressing genuine flanks.
    return cyclic_match(sequence, motif)


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
    template: RepeatTemplate | None = None


def pair_has_anchor(pair: ReadPair, template: RepeatTemplate) -> bool:
    """Whether classify_pair would retain this pair, without building evidence.

    Every accepted flank hit has a positive score, so an orientation with an
    anchor always outranks one without anchors. One hit therefore proves that
    classify_pair returns evidence, regardless of its eventual evidence tags.
    Repeat-only orientations remain excluded exactly as in classify_pair.
    """
    return pair_anchor_score(pair, template) > 0


@lru_cache(maxsize=4096)
def _primer_hit(sequence: str, primer: str):
    if not sequence or not primer:
        return None
    if not set(primer)-set('ACGT'):
        exact = sequence.find(primer)
        if exact >= 0:
            return FlankHit(exact, exact+len(primer), 0, len(primer), 2*len(primer), 0)
        # With k edits, one of k+1 disjoint primer pieces must be exact,
        # including indels. Native substring searches reject most background.
        pieces = min(2, max(1, len(primer)//5)) + 1
        if not any(primer[i*len(primer)//pieces:(i+1)*len(primer)//pieces] in sequence
                   for i in range(pieces)):
            return None
    pattern = _primer_pattern(primer)
    best = None
    # Separate concrete runs so an unknown read base cannot become a fuzzy
    # wildcard match. Matching stays in-process, including degenerate primers.
    for segment in regex.finditer('[ACGT]+', sequence):
        hit = pattern.search(segment.group())
        if hit:
            candidate = (sum(hit.fuzzy_counts), segment.start()+hit.start(), segment.start()+hit.end())
            if best is None or candidate < best:
                best = candidate
    if best is not None:
        edits, start, end = best
        return FlankHit(start, end, 0, len(primer), max(1, 2*len(primer)-6*edits), edits)


@lru_cache(maxsize=256)
def _primer_pattern(primer):
    from .locus_measurement import _IUPAC
    pattern = ''.join('[' + ''.join(sorted(_IUPAC.get(base, {base}))) + ']' for base in primer)
    return regex.compile('(?:' + pattern + '){e<=' + str(min(2, max(1, len(primer)//5))) + '}', regex.BESTMATCH)


@lru_cache(maxsize=4096)
def _oriented_anchors(sequence: str, motif: str, left_flank: str, right_flank: str, primer_only=False):
    """Cache quality-independent anchors for scoring and final classification."""
    choices = []
    for strand, seq in (("+", sequence.upper()), ("-", revcomp(sequence))):
        if primer_only:
            left, right = _primer_hit(seq, left_flank), _primer_hit(seq, right_flank)
            choices.append((sum(h.score for h in (left, right) if h), strand, seq, left, right, False))
            continue
        is_repeat = repetitive(seq, motif)
        left = None if is_repeat else flank_hit(seq, left_flank)
        right = None if is_repeat else flank_hit(seq, right_flank)
        # Short local matches occur by chance in whole-genome background.
        # Rescue 10..19-base anchors only when adjacent motif sequence supplies
        # an independently recognizable repeat boundary.
        if left and left.query_end-left.query_start < 20 and len(left_flank) >= 20:
            tail = seq[left.query_end:left.query_end+20]
            if len(left_flank)-left.flank_end > 3 or len(tail) < 10 or not repetitive(tail, motif):
                left = None
        if right and right.query_end-right.query_start < 20 and len(right_flank) >= 20:
            head = seq[max(0, right.query_start-20):right.query_start]
            if right.flank_start > 3 or len(head) < 10 or not repetitive(head, motif):
                right = None
        choices.append((sum(h.score for h in (left, right) if h), strand, seq, left, right, is_repeat))
    return max(choices, key=lambda v: (v[0], v[5]))


def pair_anchor_score(pair: ReadPair, template: RepeatTemplate) -> int:
    """Exact recruitment score without repeat geometry or evidence objects."""
    score = 0
    for read in (pair.read1, pair.read2):
        if read is not None:
            _, _, _, left, right, _ = _oriented_anchors(
                read.sequence, template.motif, template.left, template.right, template.primer_only)
            score += max(left.score if left else 0, right.score if right else 0)
    return score


def classify_pair(pair: ReadPair, template: RepeatTemplate, repeat_threshold: float = .7) -> MoleculeEvidence | None:
    """Anchor to non-repeat sequence; count each pair once, with multiple tags."""
    reads = [pair.read1] + ([pair.read2] if pair.read2 else [])
    oriented = []
    hits = []
    strands = []
    repeat_only = []
    qualities = []
    for read in reads:
        _, strand, seq, left, right, is_repeat = _oriented_anchors(
            read.sequence, template.motif, template.left, template.right, template.primer_only)
        oriented.append(seq)
        hits.append((left, right))
        strands.append(strand)
        qualities.append(read.quality[::-1] if strand == "-" and read.quality else read.quality)
        repeat_only.append(is_repeat)
    if not any(left or right for left, right in hits):
        return None  # Repeat-only molecules cannot identify a microbial locus.
    if template.primer_only:
        if len(reads) == 2:
            for i, j in ((0, 1), (1, 0)):
                if not any(hits[i]) and any(hits[j]):
                    strands[i] = '-' if strands[j] == '+' else '+'
                    oriented[i] = revcomp(reads[i].sequence) if strands[i] == '-' else reads[i].sequence.upper()
                    qualities[i] = reads[i].quality[::-1] if strands[i] == '-' and reads[i].quality else reads[i].quality
        products = set()
        for seq, quality, (left, right) in zip(oriented, qualities, hits):
            if left and right and left.query_end <= right.query_start:
                start, end = left.query_start, right.query_end
                if quality is None or min(quality[start:end], default='I') >= '5':
                    products.add(seq[start:end])
        classes = ('discordant',) if len(products) > 1 else ('FULL_SPAN',) if products else ('UNINFORMATIVE',)
        return MoleculeEvidence(pair.molecule_id, template.locus.locus_id, classes,
            tuple(oriented), tuple(strands), product_sequence=next(iter(products)) if len(products) == 1 else '',
            alignment={str(i): {side: vars(hit) if hit else None for side, hit in zip(('left', 'right'), pair_hits)}
                       for i, pair_hits in enumerate(hits)}, qualities=tuple(qualities), template=template)
    classes = set()
    observed = []
    lower = 0.0
    product = ""
    offsets = []
    for seq, quality, (left, right) in zip(oriented, qualities, hits):
        start = left.query_end + len(template.left) - left.flank_end if left else None
        end = right.query_start - right.flank_start if right else None
        # Internal flank matches locate mates for fragment geometry, but do not
        # observe a repeat boundary. Allow only terminal clipping of up to
        # three bases for boundary evidence, not extrapolation across a flank.
        left_boundary = left is not None and len(template.left) - left.flank_end <= 3
        right_boundary = right is not None and right.flank_start <= 3
        if left_boundary and right_boundary and 0 <= start <= end <= len(seq):
            classes.add("FULL_SPAN")
            observed.append((end - start) / template.unit)
            # Primer positions, rather than reference flank lengths, bound
            # physical products even when an observed flank contains an indel.
            candidate = primer_product(seq, template)
            begin = seq.find(candidate) if candidate else -1
            finish = begin + len(candidate)
            if candidate and (quality is None or min(quality[begin:finish], default="I") >= "5"):
                product = candidate
        else:
            segment = seq[start:] if left_boundary and 0 <= start < len(seq) else (
                seq[:end] if right_boundary and 0 < end <= len(seq) else "")
            if segment and repetitive(segment, template.motif):
                classes.add("LEFT_BOUNDARY" if left_boundary else "RIGHT_BOUNDARY")
                lower = max(lower, len(segment) / template.unit)
        if left_boundary and start is not None and start < len(seq) and not right_boundary:
            classes.add("SOFTCLIP_LEFT")
        if right_boundary and end is not None and end > 0 and not left_boundary:
            classes.add("SOFTCLIP_RIGHT")
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
                    classes.add("FLANK_PAIR")
                    break
                fragment_offset = None
        if any(repeat_fraction(seq, template.motif) >= repeat_threshold for seq in oriented):
            classes.add("ANCHORED_REPEAT")
            lower = max(lower, max((len(seq) / template.unit for seq, rep in zip(oriented, repeat_only) if rep), default=0))
        if not classes and any(hits[0]) and any(hits[1]) and strands[0] == strands[1]:
            classes.add("discordant")
    if any(repeat_fraction(seq, template.motif) >= repeat_threshold for seq in oriented):
        classes.add("REPEAT_RICH")
    if not classes:
        classes.add("UNINFORMATIVE")
    if observed and max(observed) - min(observed) > 0.5:
        classes = {"discordant"}
        observed = []
        product = ""
    return MoleculeEvidence(pair.molecule_id, template.locus.locus_id, tuple(sorted(classes)),
                            tuple(oriented), tuple(strands),
                            float(np.median(observed)) if observed else None, lower,
                            fragment_offset, product,
                            {str(i): {side: vars(hit) if hit else None for side, hit in zip(("left", "right"), pair_hits)}
                             for i, pair_hits in enumerate(hits)}, tuple(qualities), template)


@dataclass(frozen=True)
class InsertDistribution:
    mean: float
    sd: float
    pairs_used: int
    source: str
    median: float | None = None
    mad: float | None = None
    q05: float | None = None
    q95: float | None = None


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
    return InsertDistribution(float(np.mean(values)), max(1.0, float(np.std(values))), len(values), "robust_flank_pairs", float(median), mad,
                              float(np.quantile(values, .05)), float(np.quantile(values, .95)))


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



import re
from .alignment_evidence import CandidateAlignment
from .candidate_contexts import CandidateContext

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def _repeat_indel_delta(alignment: CandidateAlignment, context: CandidateContext) -> float:
    reference = alignment.reference_start
    delta = 0
    unit = max(context.repeat_unit_length, 1)
    for length_text, operation in _CIGAR.findall(alignment.cigar):
        length = int(length_text)
        if operation in "M=X":
            reference += length
        elif operation == "D":
            if context.repeat_start - unit <= reference <= context.repeat_end + unit and length % unit == 0:
                delta -= length / unit
            reference += length
        elif operation == "I":
            if context.repeat_start - unit <= reference <= context.repeat_end + unit and length % unit == 0:
                delta += length / unit
    return delta
