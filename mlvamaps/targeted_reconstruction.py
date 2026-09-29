"""Direct observation, bounded overlap reconstruction, then haploid inference.

No reference sequence is emitted as an observed product. Inferred intervals are
represented by N bases so repeat number cannot fabricate intra-repeat SNPs.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from functools import lru_cache
import math

import numpy as np
import parasail
import regex

from .locus_products import LocusProduct, product_repeat_allele
from .short_read_evidence import MoleculeEvidence, RepeatTemplate, InsertDistribution, _read_profile, repetitive, primer_product, primer_bounds

EVIDENCE_CLASSES = ('FULL_SPAN', 'LEFT_BOUNDARY', 'RIGHT_BOUNDARY', 'FLANK_PAIR',
                    'REPEAT_RICH', 'ANCHORED_REPEAT', 'SOFTCLIP_LEFT', 'SOFTCLIP_RIGHT', 'UNINFORMATIVE')


@dataclass
class LocusRecovery:
    products: list[LocusProduct] = field(default_factory=list)
    method: str = 'NO_CALL'
    confidence: float = 0.0
    states: np.ndarray = field(default_factory=lambda: np.array([]))
    log_likelihoods: np.ndarray = field(default_factory=lambda: np.array([]))
    posterior: np.ndarray = field(default_factory=lambda: np.array([]))
    best: float | None = None
    second: float | None = None
    interval: tuple[float, float] | None = None
    identifiable: bool = False
    limit_reached: bool = False
    counts: dict[str, int] = field(default_factory=dict)
    reason: str = ''


@lru_cache(maxsize=4096)
def _overlaps(left, right, motif, minimum=20, unit=0, mismatch_fraction=0):
    # Every valid offset is considered. A repeat-only overlap is incapable of
    # determining the number of traversals and must never stitch an allele.
    # Without a supplied motif, require enough observed sequence to test two
    # repeat units. Shorter overlaps cannot rule out a repeat-only join.
    candidates = range(minimum, min(len(left), len(right)) + 1)
    if mismatch_fraction:
        candidates = set()
        # At 2% substitutions, every accepted >=20 bp overlap contains an
        # intact 15-base block. Native substring searches bound verification.
        for start in range(0, len(right)-14, 15):
            seed = right[start:start+15]
            position = left.find(seed)
            while position >= 0:
                size = len(left)-position+start
                if minimum <= size <= min(len(left), len(right)):
                    candidates.add(size)
                position = left.find(seed, position+1)
    return tuple(size for size in sorted(candidates)
                 if (left[-size:] == right[:size] or mismatch_fraction and
                     sum(a != b for a,b in zip(left[-size:], right[:size])) <= int(size*mismatch_fraction))
                 and (bool(motif) or unit > 0 and size >= 2*unit)
                 and not repetitive(left[-size:], motif or right[:unit])
                 and not repetitive(right[:size], motif or right[:unit]))


def merge_molecule(item, template):
    if len(item.sequences) != 2 or set(item.orientations) != {'+', '-'}:
        return ''
    i = item.orientations.index('+')
    left, right = item.sequences[i], item.sequences[1-i]
    overlaps = _overlaps(left, right, template.motif, unit=template.unit)
    if len(overlaps) != 1:
        return ''
    merged = left + right[overlaps[0]:]
    return primer_product(merged, template)


def _reliable(item):
    return all(q is None or min(q, default='I') >= '5' for q in item.qualities)


def direct_products(items, template, sample_id, min_fraction, min_secondary_reads):
    groups = defaultdict(list)
    for item in items:
        if not item.product_sequence and not _reliable(item):
            continue
        sequence = item.product_sequence or merge_molecule(item, item.template or template)
        if sequence:
            item.classes = tuple(sorted(set(item.classes) | {"FULL_SPAN"}))
        if not sequence and item.observed_repeat is not None:
            # A flank-bounded span is direct length evidence. Keep observed bases
            # and mask the unobserved primer/flank tails, never copy reference SNPs.
            for i, seq in enumerate(item.sequences):
                hits = item.alignment.get(str(i), {})
                left, right = hits.get('left'), hits.get('right')
                context = item.template or template
                if left and right:
                    start = left['query_start'] - left['flank_start']
                    end = right['query_end'] + len(context.right) - right['flank_end']
                    sequence = 'N'*max(0, -start) + seq[max(0, start):min(len(seq), end)] + 'N'*max(0, end-len(seq))
                    break
        if sequence:
            groups[sequence].append(item)
    # Unlinked singleton sequencing errors should not choose an arbitrary SNP
    # haplotype. A same-length consensus is allowed only without a supported
    # exact haplotype; minority bases remain unknown below 70% agreement.
    by_length = defaultdict(list)
    for seq, members in groups.items():
        by_length[len(seq)].append((seq, members))
    for length, rows in by_length.items():
        depth = sum(len(members) for _, members in rows)
        if depth >= 3 and max(len(members) for _, members in rows) < min_secondary_reads:
            consensus = []
            for position in range(length):
                votes = Counter()
                for seq, members in rows:
                    votes[seq[position]] += len(members)
                base, support = votes.most_common(1)[0]
                consensus.append(base if support/depth >= .7 else 'N')
            members = [item for _, group in rows for item in group]
            for seq, _ in rows:
                del groups[seq]
            groups[''.join(consensus)] = members
    total = sum(map(len, groups.values()))
    products = []
    for seq, members in sorted(groups.items(), key=lambda v: (-len(v[1]), len(v[0]), v[0])):
        fraction = len(members) / total
        products.append(LocusProduct(sample_id, template.locus.locus_id, 'short_read',
            f'{template.locus.locus_id}|v{len(products)+1}', seq, support_count=len(members),
            effective_depth=len(members), estimated_fraction=fraction,
            reconstruction_confidence=1 - .01**len(members),
            evidence={'call_method': 'DIRECT', 'molecule_ids': [m.molecule_id for m in members],
                      'meaningful': 'yes' if len(members) >= min_secondary_reads and fraction >= min_fraction else 'no',
                      'uncovered_bases': seq.count('N')}))
    # Conflicting singleton lengths are not a confidently observed haploid call.
    support_by_repeat = Counter()
    for product in products:
        support_by_repeat[product_repeat_allele(product, template.locus)[1]] += product.support_count
    substantial = [n for n in support_by_repeat.values() if n >= min_secondary_reads and n/total >= min_fraction]
    if products and len(substantial) < 2:
        from dataclasses import replace
        products = [replace(p, reconstruction_confidence=p.reconstruction_confidence *
                    support_by_repeat[product_repeat_allele(p, template.locus)[1]]/total) for p in products]
    return products


def _pileup_consensus(rows, length):
    votes = [dict() for _ in range(length)]
    for sequence, offset, members in rows:
        for position, base in enumerate(sequence, offset):
            for molecule in members:
                previous = votes[position].get(molecule, base)
                votes[position][molecule] = base if previous == base else 'N'
    consensus = []
    for column in votes:
        base, support = Counter(column.values()).most_common(1)[0] if column else ('N', 0)
        consensus.append(base if column and support/len(column) >= .7 else 'N')
    return ''.join(consensus)


def _pileup_sequences(sequences):
    """Consolidate substitution errors at uniquely anchored read coordinates.

    Same-length, ungapped placements preserve observed lengths. An indel or an
    ambiguous repeat offset remains a separate graph node, never a forced join.
    Molecules vote once per position, including overlapping mates.
    """
    index = defaultdict(list)
    representatives, groups = [], []
    for sequence in sorted(sequences, key=lambda s: (-len(s), -len(sequences[s]), s)):
        placements = set()
        for start in range(0, len(sequence)-14, 15):
            for node, position in index.get(sequence[start:start+15], ()):
                offset = position-start
                if offset == 0 and len(sequence) == len(representatives[node]):
                    placements.add((node, offset))
        matches = []
        for node, offset in placements:
            target = representatives[node][offset:offset+len(sequence)]
            errors = sum(a != b for a, b in zip(sequence, target))
            if errors <= int(.02*len(sequence)):
                matches.append((errors, node, offset))
        matches.sort()
        if matches and (len(matches) == 1 or matches[0][0] < matches[1][0]):
            _, node, offset = matches[0]
            groups[node].append((sequence, offset))
        else:
            node = len(representatives)
            representatives.append(sequence)
            groups.append([(sequence, 0)])
            positions = defaultdict(list)
            for start in range(len(sequence)-14):
                positions[sequence[start:start+15]].append(start)
            for seed, starts in positions.items():
                if len(starts) == 1 and 'N' not in seed:
                    index[seed].append((node, starts[0]))
    result = defaultdict(set)
    for representative, group in zip(representatives, groups):
        consensus = _pileup_consensus(((sequence, offset, sequences[sequence])
                                      for sequence, offset in group), len(representative))
        result[consensus].update(set().union(*(sequences[sequence] for sequence, _ in group)))
    return result


def _compatible_pairs(items, sequence, insert):
    for item in items:
        if len(item.sequences) != 2 or set(item.orientations) != {'+', '-'}:
            continue
        first = item.orientations.index('+')
        a, b = item.sequences[first], item.sequences[1-first]
        x, y = sequence.find(a), sequence.find(b)
        if x >= 0 and y >= 0 and sequence.find(a, x+1) < 0 and sequence.find(b, y+1) < 0:
            fragment = y + len(b) - x
            if y < x or (insert and abs(fragment-insert.mean) > 4*insert.sd):
                return False
    return True


def microassemble(items, template, insert=None, max_nodes=2048, max_paths=64):
    """Seeded overlap graph with unique offsets and bounded paths.

    High-depth reads first form coordinate pileups; limits apply to compacted
    sequences. Seed-indexed overlaps and shared coordinate layouts keep redundant
    coverage from consuming the graph/path budget.
    ponytail: at most 2048 compacted nodes and 64 alternative path expansions; complex
    unresolved graphs still defer to inference, not arbitrary repeat joins.
    Repeat-only edges and multiple overlap offsets are rejected, not resolved
    using reference length. All accepted contigs must contain both primers.
    """
    sequences = defaultdict(set)
    for item in items:
        for index, seq in enumerate(item.sequences):
            quality = item.qualities[index] if item.qualities else None
            if quality is None:
                sequences[seq].add(item.molecule_id)
            else:
                for segment in regex.finditer(r'[5-~]{20,}', quality):
                    sequences[seq[segment.start():segment.end()]].add(item.molecule_id)
    if max_nodes <= 0:
        return '', [], 'assembly_node_limit'
    if len(sequences) > 256:
        sequences = _pileup_sequences(sequences)
    pileup_rows = {sequence: [(sequence, 0, set(members))] for sequence, members in sequences.items()}
    nodes = sorted(sequences, key=lambda s: (-len(s), s))
    retained = []
    for sequence in nodes:
        container = next((larger for larger in retained if sequence in larger), None)
        if container is not None:
            offset = container.find(sequence)
            if container.find(sequence, offset+1) < 0:
                # Contained reads vote only where they were observed, not
                # across the containing read's unsupported extensions.
                pileup_rows[container].extend((read, offset+start, members)
                                             for read, start, members in pileup_rows[sequence])
                sequences[container].update(sequences[sequence])
        else:
            retained.append(sequence)
    nodes = retained
    if len(nodes) > max_nodes:
        return '', [], 'assembly_node_limit'
    seed_index = defaultdict(set)
    for j, right in enumerate(nodes):
        for start in range(0, len(right)-14, 15):
            seed_index[right[start:start+15]].add(j)
    edges = defaultdict(list)
    ambiguous = False
    for i, left in enumerate(nodes):
        candidates = set()
        for start in range(len(left)-14):
            candidates.update(seed_index.get(left[start:start+15], ()))
        for j in sorted(candidates):
            if i == j:
                continue
            right = nodes[j]
            # Full coverage of the target adds no sequence. Near-duplicate
            # reads can otherwise form reciprocal edges at the same genomic
            # position, falsely appearing to be a repeat traversal cycle.
            overlaps = _overlaps(left, right, template.motif, unit=template.unit, mismatch_fraction=.02)
            if len(right) in overlaps:
                continue
            if len(overlaps) == 1:
                edges[i].append((j, overlaps[0]))
            elif len(overlaps) > 1:
                ambiguous = True
    # Ambiguous repeat-only joins are omitted. They must not veto an
    # independently anchored complete path elsewhere in a high-depth pool.
    # A consistent overlap component is already a coordinate system. Pile up
    # all its reads once instead of enumerating thousands of equivalent paths.
    neighbors = defaultdict(list)
    for i, outgoing in edges.items():
        for j, overlap in outgoing:
            offset = len(nodes[i])-overlap
            neighbors[i].append((j, offset))
            neighbors[j].append((i, -offset))
    unseen = set(range(len(nodes)))
    layouts, consistent = [], True
    while unseen:
        seed = min(unseen)
        positions, pending = {seed: 0}, [seed]
        unseen.remove(seed)
        while pending:
            i = pending.pop()
            for j, offset in neighbors[i]:
                position = positions[i]+offset
                if j in positions:
                    consistent &= positions[j] == position
                else:
                    positions[j] = position
                    unseen.discard(j)
                    pending.append(j)
        layouts.append(positions)
    if consistent:
        products = {}
        for positions in layouts:
            start = min(positions.values())
            length = max(position+len(nodes[i]) for i, position in positions.items())-start
            sequence = _pileup_consensus(((read, position-start+offset, members)
                                         for i, position in positions.items()
                                         for read, offset, members in pileup_rows[nodes[i]]), length)
            product = primer_product(sequence, template)
            if product and _compatible_pairs(items, sequence, insert):
                products.setdefault(product, set()).update(set().union(*(sequences[nodes[i]] for i in positions)))
        if len(products) == 1:
            product, members = next(iter(products.items()))
            if len(members) >= 2:
                return product, sorted(members), ''
        if products:
            return '', [], 'multiple_reconstructions'
        return '', [], 'ambiguous_overlap_offsets' if ambiguous else 'no_complete_path'
    starts = [i for i, seq in enumerate(nodes) if primer_bounds(seq, template.locus.forward_primer)]
    paths = [(i, nodes[i], ((i, 0),)) for i in starts]
    products = {}
    visited = set()
    while paths:
        i, sequence, used = paths.pop()
        key = i, sequence, tuple((j, offset) for j, offset in used if nodes[j] not in sequence)
        if key in visited:
            continue
        visited.add(key)
        if len(visited) > max_paths:
            return '', [], 'assembly_path_limit'
        product = primer_product(sequence, template)
        if product:
            if any(nodes[j] != sequence[offset:offset+len(nodes[j])] for j, offset in used):
                # Overlaps already establish physical coordinates. Realigning
                # these reads with POA can delete a repeat traversal; pile up
                # at the validated offsets and retain uncertain bases as N.
                sequence = _pileup_consensus(((read, offset+start, members)
                                             for j, offset in used
                                             for read, start, members in pileup_rows[nodes[j]]), len(sequence))
                product = primer_product(sequence, template)
                if not product:
                    continue
            if _compatible_pairs(items, sequence, insert):
                members = set().union(*(sequences[nodes[j]] for j, _ in used))
                products.setdefault(product, set()).update(members)
            continue
        for j, overlap in edges.get(i, ()):
            if any(j == previous for previous, _ in used):
                return '', [], 'cyclic_overlap_graph'
            paths.append((j, sequence + nodes[j][overlap:], used + ((j, len(sequence)-overlap),)))
    if products and len({len(product) for product in products}) == 1:
        # Base uncertainty need not discard an agreed primer-bounded length.
        # Keep primer ends identical so this cannot hide differences in PCR
        # size calibration. Mask disputed interior bases rather than selecting
        # an arbitrary path/SNP or counting shared molecules more than once.
        ends = {(p[:len(template.locus.forward_primer)],
                 p[-len(template.locus.reverse_primer):]) for p in products}
        if len(ends) > 1:
            return '', [], 'multiple_reconstructions'
        sequence = ''.join(column[0] if len(set(column)) == 1 else 'N'
                           for column in zip(*products))
        members = sorted(set().union(*products.values()))
        if len(members) >= 2:
            return sequence, members, ''
    return '', [], 'multiple_reconstructions' if products else 'ambiguous_overlap_offsets' if ambiguous else 'no_complete_path'


def candidate_likelihood(items, template, insert, maximum, minimum_probability, context_templates=None):
    """Haploid, sequence/fragment log likelihood; no genotype or component EM.

    Marginalize over distinct recruited natural contexts for each repeat length.
    Bound-only observations remain uncertain even at very high read depth.
    """
    if maximum < 1 or maximum * template.unit > 1_000_000:
        raise ValueError('invalid candidate ceiling or candidate locus exceeds sequence safety limit')
    upper = min(maximum, max(10, template.locus.expected_max_repeats + 2,
                             math.ceil(max([e.lower_bound for e in items] + [0])) + 3))
    if insert:
        upper = min(maximum, max(upper, math.ceil((insert.mean + 4*insert.sd)/template.unit)))
    step = .5 if template.unit > 1 else 1.
    context_templates = context_templates or [e.template for e in items if e.template]
    contexts = list(dict.fromkeys((t.left, t.right, t.motif) for t in context_templates))
    contexts = contexts or [(template.left, template.right, template.motif)]
    if len(contexts) > 64:
        return LocusRecovery(method='AMBIGUOUS', reason='candidate_context_limit')
    cache = {}
    while True:
        states = np.arange(0, upper + step/2, step)
        joint = np.zeros(len(states))
        if len(states) * (len(template.left)+len(template.right)+upper*template.unit) > 50_000_000:
            raise ValueError('candidate sequences exceed memory safety limit')
        for item in items:
            # Flank-only/unsupported clips establish presence, not length.
            # Including their local scores can bias length through chance hits.
            if not item.lower_bound and 'FLANK_PAIR' not in item.classes:
                continue
            scores_for_molecule = np.zeros(len(states))
            if item.lower_bound:
                scores_for_molecule[states < item.lower_bound-.25] = -50
            if insert and 'FLANK_PAIR' in item.classes and item.fragment_offset is not None:
                fragment = item.fragment_offset + states * template.unit
                # Student-t tails reduce the influence of discordant fragments.
                scores_for_molecule -= 2.5 * np.log1p(((fragment-insert.mean)/insert.sd)**2/4)
            for read in item.sequences:
                if read not in cache or len(cache[read]) != len(states):
                    scores = list(cache.get(read, []))
                    for state in states[len(scores):]:
                        backgrounds = []
                        for left, right, motif in contexts:
                            length = round(state*template.unit)
                            target = left + (motif*math.ceil(length/len(motif)))[:length] + right
                            score = 2*len(read) if read in target else parasail.sw_striped_profile_32(_read_profile(read), target, 5, 1).score
                            backgrounds.append((score-2*len(read))/8)
                        top = max(backgrounds)
                        scores.append(top + math.log(sum(math.exp(v-top) for v in backgrounds)/len(backgrounds)))
                    cache[read] = np.asarray(scores)
                scores_for_molecule += cache[read]
            joint += scores_for_molecule
        posterior = np.exp(np.maximum(joint-joint.max(), -700))
        posterior /= posterior.sum()
        if upper >= maximum or posterior[-2:].sum() < .01:
            break
        upper = min(maximum, upper*2)
    ranking = np.argsort(-posterior, kind='stable')
    best, second = map(int, ranking[:2])
    selected = ranking[:np.searchsorted(np.cumsum(posterior[ranking]), .95)+1]
    limit = upper == maximum and posterior[-2:].sum() >= .01
    identifiable = bool(insert and any('FLANK_PAIR' in e.classes for e in items))
    # Detect separated insert-implied clusters without fitting a mixture model.
    implied = [round((insert.mean-e.fragment_offset)/template.unit) for e in items
               if insert and 'FLANK_PAIR' in e.classes and e.fragment_offset is not None]
    modes = Counter(implied).most_common()
    mixed = len(modes)>1 and modes[1][1] >= max(3, .2*len(implied)) and abs(modes[0][0]-modes[1][0])*template.unit > 4*insert.sd
    confidence = float(posterior[best])
    informative = np.ptp(joint) > 1e-8
    tied = np.count_nonzero(np.isclose(joint, joint[best], rtol=0, atol=1e-8)) > 1
    return LocusRecovery(method='MIXED' if mixed else 'INFERRED' if identifiable and confidence >= minimum_probability and not limit else 'AMBIGUOUS',
        confidence=confidence if informative else 0, states=states, log_likelihoods=joint, posterior=posterior,
        best=float(states[best]) if informative else None, second=float(states[second]) if informative else None,
        interval=(float(states[selected].min()), float(states[selected].max())) if informative else None,
        identifiable=informative and identifiable and confidence >= minimum_probability and not limit and not mixed,
        limit_reached=limit,
        reason='no_repeat_length_information' if not informative else
               'minimum_likelihood_tie; repeat_length_lower_bound' if tied and not insert else
               'minimum_likelihood_tie' if tied else
               'candidate_limit_reached' if limit else 'best_likelihood_estimate')


def recover_locus(items, template, sample_id, insert=None, maximum=100, minimum_probability=.8,
                  minimum_spanning_pairs=2, min_fraction=.01, min_secondary_reads=2, context_templates=None):
    counts = dict(Counter(c for e in items for c in e.classes))
    usable = [e for e in items if 'discordant' not in e.classes]
    if not usable or template is None:
        return LocusRecovery(counts=counts)
    products = direct_products(usable, template, sample_id, min_fraction, min_secondary_reads)
    counts = dict(Counter(c for e in items for c in e.classes))
    if products and all('N' in product.sequence for product in products):
        sequence, members, _ = microassemble(usable, template, insert)
        # Existing direct length evidence must agree with the reconstructed path.
        if sequence and len(sequence) == len(products[0].sequence):
            product = LocusProduct(sample_id, template.locus.locus_id, 'short_read', template.locus.locus_id+'|v1',
                sequence, support_count=len(members), effective_depth=len(members), reconstruction_confidence=.99,
                evidence={'call_method': 'RECONSTRUCTED', 'meaningful': 'yes', 'molecule_ids': members})
            return LocusRecovery(products=[product], method='RECONSTRUCTED', confidence=.99, identifiable=True, counts=counts)
    if products:
        meaningful = [p for p in products if p.evidence['meaningful'] == 'yes']
        return LocusRecovery(products=products,
                             method='MIXED' if len(meaningful)>1 else 'DIRECT' if products[0].reconstruction_confidence >= minimum_probability else 'AMBIGUOUS',
                             confidence=products[0].reconstruction_confidence,
                             identifiable=products[0].reconstruction_confidence >= minimum_probability, counts=counts)
    sequence, members, reason = microassemble(usable, template, insert)
    if sequence:
        product = LocusProduct(sample_id, template.locus.locus_id, 'short_read', template.locus.locus_id+'|v1',
            sequence, support_count=len(members), effective_depth=len(members), reconstruction_confidence=.99,
            evidence={'call_method': 'RECONSTRUCTED', 'meaningful': 'yes', 'molecule_ids': members})
        return LocusRecovery(products=[product], method='RECONSTRUCTED', confidence=.99, identifiable=True, counts=counts)
    if template.primer_only:
        return LocusRecovery(method='AMBIGUOUS', counts=counts,
            reason='primer_only_no_complete_product:' + reason)
    result = candidate_likelihood(usable, template, insert, maximum, minimum_probability, context_templates)
    result.counts, result.reason = counts, result.reason or reason
    if counts.get('FLANK_PAIR', 0) < minimum_spanning_pairs:
        result.identifiable = False
        if result.method == 'INFERRED':
            result.method = 'AMBIGUOUS'
    if result.identifiable:
        # Partial SNPs use quality-filtered unique flank placements; uncovered
        # positions (including the unobserved repeat) stay explicitly unknown.
        from .locus_reconstruction import _project_flanks
        synthetic = template.sequence(result.best)
        votes = defaultdict(Counter)
        for item in usable:
            for pos, allele in _project_flanks(item, synthetic, template, result.best).items():
                votes[pos][allele] += 1
        seq = []
        for pos in range(len(synthetic)):
            if votes[pos]:
                (insertion, base), support = votes[pos].most_common(1)[0]
                seq.append(insertion + ('' if base == '-' else base) if support >= 3 and support/sum(votes[pos].values()) >= .9 else 'N')
            else:
                seq.append('N')
        result.products = [LocusProduct(sample_id, template.locus.locus_id, 'short_read', template.locus.locus_id+'|v1',
            ''.join(seq), support_count=len(usable), effective_depth=len(usable), reconstruction_confidence=result.confidence,
            repeat_likelihoods=dict(zip(map(float, result.states), map(float, result.posterior))),
            evidence={'call_method': 'INFERRED', 'meaningful': 'yes', 'molecule_ids': [e.molecule_id for e in usable],
                      'uncovered_bases': ''.join(seq).count('N')})]
    return result
