"""Bounded, sample-derived recruitment and repeat-graph construction.

No reference sequences or nominal repeat counts are used to learn anchors.
The graph is represented by two observed arms and a variable repeat edge;
targeted_reconstruction scores its bounded paths with native alignments.
"""
from collections import Counter, defaultdict
from functools import lru_cache
from os.path import commonprefix

import regex

from .models import ReadPair, ReadRecord
from .sequence import revcomp
from .short_read_evidence import MoleculeEvidence, RepeatTemplate, classify_pair, cyclic_match, primer_bounds
from .targeted_reconstruction import _overlaps, _pileup_consensuses, _pileup_sequences


def reliable_sequences(items):
    sequences = defaultdict(set)
    for item in items:
        if 'discordant' in item.classes:
            continue
        for i, sequence in enumerate(item.sequences):
            quality = item.qualities[i] if item.qualities else None
            for match in regex.finditer('[5-~]{21,}', quality or 'I'*len(sequence)):
                sequences[sequence[match.start():match.end()]].add(item.molecule_id)
    return _pileup_sequences(sequences) if len(sequences) > 256 else sequences


@lru_cache(maxsize=8192)
def repeat_runs(sequence, unit):
    """Discover tandem units by lag identity; tolerate isolated substitutions.

    A run needs two units and >=90% identity to its phase consensus. Windows
    with >=90% lag identity locate candidates, including imperfect long units.
    """
    if unit < 1 or len(sequence) < 2*unit or unit > 500:
        return ()
    equal = [a == b and a in 'ACGT' for a, b in zip(sequence, sequence[unit:])]
    score = sum(equal[:unit])
    good = []
    for start in range(len(equal)-unit+1):
        if start:
            score += equal[start+unit-1] - equal[start-1]
        good.append(score >= unit - unit//10)
    intervals = []
    for match in regex.finditer('1+', ''.join('1' if value else '0' for value in good)):
        start, end = match.start(), match.end()-1+2*unit
        segment = sequence[start:end]
        motif = ''.join(Counter(segment[i::unit]).most_common(1)[0][0] for i in range(unit))
        if sum(base == motif[i % unit] for i, base in enumerate(segment)) >= .9*len(segment):
            # Rotations are equivalent motifs; strand is already locus-oriented.
            canonical = min(motif[i:]+motif[:i] for i in range(unit))
            intervals.append((start, end, canonical))
    return tuple(intervals)


def _arms(sequences, template, motif):
    """Extend primer anchors only through agreed non-periodic overlaps."""
    def extend(pool, primer, repeat):
        starts = []
        for seq, members in pool.items():
            bounds = primer_bounds(seq, primer)
            if bounds:
                starts.append((seq[bounds[0]:], 0, members))
        if not starts:
            return ''
        # Short primer fragments must not truncate independently observed
        # boundaries. HTSlib counts each read only at the bases it covers.
        length = max(len(seq) for seq, _, _ in starts)
        arm = _pileup_consensuses([(starts, length)])[0].split('N', 1)[0]
        # ponytail: at most eight extensions and 256 pool sequences. This is
        # a conservative arm builder, not a whole-genome assembler.
        for _ in range(8):
            extensions = []
            for seq in pool:
                offsets = _overlaps(arm, seq, repeat, unit=template.unit)
                if len(offsets) == 1 and offsets[0] < len(seq):
                    extensions.append(arm + seq[offsets[0]:])
            if not extensions:
                break
            extended = commonprefix(extensions)
            if len(extended) <= len(arm):
                break
            arm = extended
        return arm
    left = extend(sequences, template.locus.forward_primer, motif)
    right = extend({revcomp(s): members for s, members in sequences.items()},
                   template.locus.reverse_primer, revcomp(motif))
    right = revcomp(right)
    # Reads may extend beyond the opposite primer into unrelated repeats.
    end = primer_bounds(left, revcomp(template.locus.reverse_primer))
    start = primer_bounds(right, template.locus.forward_primer)
    return left[:end[1]] if end else left, right[start[0]:] if start else right


def learn_template(items, template, sequences=None):
    """Require independently supported motif and both observed primer arms."""
    if not template.primer_only or not template.unit:
        return template, {'source': 'panel', 'motifs': []}
    sequences = reliable_sequences(items) if sequences is None else sequences
    votes = defaultdict(set)
    for seq, members in sequences.items():
        for start, end, motif in repeat_runs(seq, template.unit):
            if end-start >= max(12, 2*template.unit):
                votes[motif].update(members)
    ranked = sorted(votes, key=lambda motif: (-len(votes[motif]), motif))
    info = {'source': 'sample', 'motif_selection': 'primer_linked_arms',
            'motifs': [{'sequence': m, 'molecules': len(votes[m])} for m in ranked[:3]],
            'graph_ready': False}
    if not ranked or len(votes[ranked[0]]) < 2:
        return template, info
    # Deterministic depth-first ordering bounds noisy pools without selecting
    # a nominal allele. A capped pool may learn less sequence, never more.
    pool = {s: sequences[s] for s in sorted(sequences, key=lambda s: (-len(sequences[s]), -len(s), s))[:256]}
    left_arm, right_arm = _arms(pool, template, '')
    left_runs = repeat_runs(left_arm, template.unit)
    right_runs = repeat_runs(right_arm, template.unit)
    # Global motif abundance includes retained mates and outward-facing reads.
    # Only families observed on both primer-linked arms can identify the VNTR.
    # Permit substitutions (85% identity, as for flank anchors), never indels
    # or a change in unit length; use a common cyclic phase on both boundaries.
    candidates = []
    motifs = sorted({m for _, _, m in left_runs + right_runs},
                    key=lambda m: (-len(votes.get(m, ())), m))
    for motif in motifs:
        pattern = regex.compile(f'(?:{motif}){{s<={3*template.unit//20}}}', regex.BESTMATCH)
        if any(pattern.search(t.motif+t.motif[:-1]) for t in candidates):
            continue
        left, right = set(), set()
        for start, end, found in left_runs:
            if pattern.search(found+found[:-1]):
                hit = pattern.search(left_arm, start, min(end, start+2*template.unit-1))
                if hit and hit.start() >= len(template.locus.forward_primer):
                    left.add(hit.start())
        for start, end, found in right_runs:
            if pattern.search(found+found[:-1]):
                hit = pattern.search(right_arm, max(start, end-2*template.unit+1), end)
                if hit and len(right_arm)-hit.end() >= len(template.locus.reverse_primer):
                    right.add(hit.end())
        # Substituted repeat units can split lag-identity runs. Join their
        # boundaries only when the intervening observed tract is still periodic.
        if len(left) > 1 and cyclic_match(left_arm[min(left):max(left)+template.unit], motif):
            left = {min(left)}
        if len(right) > 1 and cyclic_match(right_arm[min(right)-template.unit:max(right)], motif):
            right = {max(right)}
        if len(left) == len(right) == 1:
            support = set().union(*(members for found, members in votes.items()
                                   if pattern.search(found+found[:-1])))
            if len(support) >= 2:
                candidates.append(RepeatTemplate(template.locus, left_arm[:left.pop()],
                                                 right_arm[right.pop():], motif, template.unit))
    if len(candidates) != 1:
        info['reason'] = 'ambiguous_sample_motif' if candidates else 'incomplete_sample_repeat_boundaries'
        return template, info
    learned = candidates[0]
    info.update(graph_ready=True, left_bp=len(learned.left), right_bp=len(learned.right), motif=learned.motif,
                left_sequence=learned.left, right_sequence=learned.right,
                structure=f'{learned.left}({learned.motif})*{learned.right}')
    return learned, info


def evidence_pair(item):
    """Restore original read orientation before applying learned flanks."""
    reads = []
    for i, (seq, strand) in enumerate(zip(item.sequences, item.orientations)):
        quality = item.qualities[i] if item.qualities else None
        reads.append(ReadRecord(item.molecule_id+f'/{i+1}', revcomp(seq) if strand == '-' else seq,
                                quality[::-1] if strand == '-' and quality else quality))
    return ReadPair(item.molecule_id, reads[0], reads[1] if len(reads) > 1 else None)


def _kmers(sequence, quality=None):
    """Rolling canonical 21-mers; ambiguous/low-quality bases break seeds."""
    forward = reverse = valid = 0
    codes = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    mask = (1 << 42)-1
    for i, base in enumerate(sequence):
        code = codes.get(base)
        if code is None or quality is not None and quality[i] < '5':
            forward = reverse = valid = 0
            continue
        forward = ((forward << 2) | code) & mask
        reverse = (reverse >> 2) | ((3-code) << 40)
        valid += 1
        if valid >= 21:
            yield i-20, min(forward, reverse), '+' if forward <= reverse else '-'


class SampleRecruiter:
    """One frozen index across all loci; require >=40 bp of anchor support.

    Seeds inside candidate tandem repeats and low-complexity seeds do not
    recruit. Both mates vote; shared-locus ties and strand conflicts abstain.
    This is exact seed mapping, so divergent reads can remain unrecruited.
    """
    def __init__(self, evidence, templates, sequence_pools=None):
        self.templates = templates
        self.index = defaultdict(dict)
        self.capped_loci = []
        for locus_id, items in sorted(evidence.items()):
            template = templates[locus_id]
            sequences = reliable_sequences(items) if sequence_pools is None else sequence_pools[locus_id]
            if len(sequences) > 256:
                self.capped_loci.append(locus_id)
            pool = sorted(sequences, key=lambda s: (-len(sequences[s]), -len(s), s))[:256]
            periodic_seeds = {(template.motif*22)[phase:phase+21] for phase in range(len(template.motif))}
            for seq in pool:
                if len(seq) < max(40, 2*template.unit):
                    continue  # Too short to exclude a whole-unit repeat anchor.
                masked = set()
                for start, end, _ in repeat_runs(seq, template.unit):
                    masked.update(range(max(0, start-20), end))
                for start, key, strand in _kmers(seq):
                    seed = seq[start:start+21]
                    if start in masked:
                        continue
                    # Also reject shorter periodic sequences and motifs whose
                    # full period is longer than the seed.
                    if any(seed == (seed[:unit]*21)[:21] for unit in range(1, 8)):
                        continue
                    if seed in periodic_seeds:
                        continue
                    entry = self.index[key]
                    previous = entry.get(locus_id, strand)
                    entry[locus_id] = strand if previous == strand else '?'

    def recruit_pair(self, pair):
        hits = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0, 0])))
        reads = [pair.read1] + ([pair.read2] if pair.read2 else [])
        for i, read in enumerate(reads):
            sequence = read.sequence.upper()
            for start, key, strand in _kmers(sequence, read.quality):
                for locus_id, target_strand in self.index.get(key, {}).items():
                    if target_strand != '?':
                        orientation = '+' if strand == target_strand else '-'
                        coverage = hits[locus_id][i][orientation]
                        coverage[0] += min(21, start+21-coverage[1])
                        coverage[1] = start+21
        ranked = []
        for locus_id, mates in hits.items():
            scores, strands = [], []
            for i in range(len(reads)):
                orientations = mates.get(i, {})
                ordered = sorted(((coverage[0], strand) for strand, coverage in orientations.items()), reverse=True)
                if len(ordered) > 1 and ordered[1][0] >= 20:
                    break
                scores.append(ordered[0][0] if ordered else 0)
                strands.append(ordered[0][1] if ordered else '')
            else:
                minimum = max(40, 2*self.templates[locus_id].unit)
                if sum(scores) >= minimum and (len(strands) == 1 or not all(strands) or strands[0] != strands[1]):
                    ranked.append((sum(scores), locus_id, strands))
        ranked.sort(key=lambda row: (-row[0], row[1]))
        if not ranked or len(ranked) > 1 and ranked[0][0]-ranked[1][0] < max(12, .2*ranked[0][0]):
            return None, [row[1] for row in ranked]
        _, locus_id, strands = ranked[0]
        if len(strands) == 2:
            for i in range(2):
                if not strands[i]:
                    strands[i] = '-' if strands[1-i] == '+' else '+'
        template = self.templates[locus_id]
        item = classify_pair(pair, template)
        if item is None:
            item = MoleculeEvidence(pair.molecule_id, locus_id, ('UNINFORMATIVE',),
                tuple(revcomp(read.sequence) if strand == '-' else read.sequence.upper()
                      for read, strand in zip(reads, strands)), tuple(strands),
                qualities=tuple(read.quality[::-1] if strand == '-' and read.quality else read.quality
                                for read, strand in zip(reads, strands)), template=template)
        return item, []


    def recruit(self, pairs):
        from .short_read_recruitment import RecruitmentChunk
        result = RecruitmentChunk()
        for pair in pairs:
            result.examined += 1
            item, ambiguous = self.recruit_pair(pair)
            if item:
                result.evidence.append(item)
            elif ambiguous:
                result.ambiguous_pairs += 1
                result.ambiguity_witnesses.append((pair.molecule_id, ambiguous[0], ambiguous[1], 'yes'))
        return result


def _reclassify_locus(task):
    locus_id, items, template = task
    updated = []
    for item in items:
        if item.template == template:
            updated.append(item)
            continue  # Rescue already classified these reads against these exact arms.
        classified = classify_pair(evidence_pair(item), template)
        # Measurements and flank coordinates belong to the template that
        # produced them. Retain unanchored reads without stale measurements.
        updated.append(classified or MoleculeEvidence(item.molecule_id, item.locus_id,
            ('UNINFORMATIVE',), item.sequences, item.orientations,
            qualities=item.qualities, template=template))
    return locus_id, updated


def recruit_sample_reads(evidence, templates, replay, progress, sample_id, max_rounds=2,
                        threads=1, audit_writer=None, locus_executor=None):
    """Freeze the index per round; skip pairs already assigned to any locus."""
    assigned = {item.molecule_id for items in evidence.values() for item in items}
    learned = dict(templates)
    metadata = {'rounds': [], 'loci': {}}
    previous_index = None
    learned_sizes, sequence_pools = {}, {}
    for round_index in range(max_rounds):
        for locus_id, template in templates.items():
            with progress.phase(f'[{sample_id}] Sample anchor learning {locus_id}',
                                f"round {round_index+1}; {len(evidence[locus_id]):,} molecules"):
                if learned_sizes.get(locus_id) != len(evidence[locus_id]):
                    sequence_pools[locus_id] = reliable_sequences(evidence[locus_id])
                    learned[locus_id], metadata['loci'][locus_id] = learn_template(
                        evidence[locus_id], template, sequence_pools[locus_id])
                    learned_sizes[locus_id] = len(evidence[locus_id])
        with progress.phase(f'[{sample_id}] Sample anchor index', f'round {round_index+1}'):
            recruiter = SampleRecruiter(evidence, learned, sequence_pools)
        if not recruiter.index:
            metadata['stop_reason'] = 'no_unique_sample_anchors'
            break
        if recruiter.index == previous_index:
            metadata['stop_reason'] = 'anchor_index_unchanged'
            break
        previous_index = recruiter.index
        stats = {'round': round_index+1, 'examined': 0, 'recruited': 0, 'ambiguous': 0,
                 'index_seeds': len(recruiter.index), 'capped_loci': recruiter.capped_loci}
        progress.step(f'[{sample_id}] Recruiting against sample-derived anchors, round {round_index+1}/{max_rounds}')
        iterator = replay()
        def unassigned_pairs():
            for pair in iterator:
                stats['examined'] += 1
                if pair.molecule_id not in assigned:
                    yield pair
        try:
            from .short_read_recruitment import recruit_short_reads
            for chunk in recruit_short_reads(unassigned_pairs(), {}, sample_id, threads,
                                             recruiter=recruiter):
                for item in chunk.evidence:
                    if item.molecule_id in assigned:
                        continue
                    evidence[item.locus_id].append(item)
                    assigned.add(item.molecule_id)
                    stats['recruited'] += 1
                    if audit_writer:
                        audit_writer.writerow([sample_id, round_index+1, item.molecule_id,
                                               'recruited', item.locus_id, ''])
                stats['ambiguous'] += chunk.ambiguous_pairs
                if audit_writer:
                    audit_writer.writerows([sample_id, round_index+1, molecule, 'ambiguous', a, b]
                                           for molecule, a, b, _ in chunk.ambiguity_witnesses)
                progress.count(f'[{sample_id}] Sample recruitment pairs scanned', stats['examined'],
                               detail=f"round {round_index+1}; {stats['recruited']:,} additional pairs")
        finally:
            if hasattr(iterator, 'close'):
                iterator.close()
        metadata['rounds'].append(stats)
        progress.step(f"[{sample_id}] Sample recruitment round {round_index+1} finished; "
                      f"{stats['examined']:,} pairs scanned; {stats['recruited']:,} additional pairs")
        if not stats['recruited']:
            metadata['stop_reason'] = 'no_additional_pairs'
            break
    else:
        metadata['stop_reason'] = 'round_limit'
    reclassify = []
    for locus_id, template in templates.items():
        with progress.phase(f'[{sample_id}] Final sample template {locus_id}',
                            f'{len(evidence[locus_id]):,} molecules'):
            if learned_sizes.get(locus_id) != len(evidence[locus_id]):
                learned[locus_id], metadata['loci'][locus_id] = learn_template(evidence[locus_id], template)
            if learned[locus_id] != template or any(item.template != learned[locus_id]
                                                     for item in evidence[locus_id]):
                reclassify.append((locus_id, evidence[locus_id], learned[locus_id]))
    # Reuse the allocated locus workers; recovering more graphs must not force
    # every newly informative molecule through serial native alignments.
    from .concurrency import bounded_ordered_map, resolve_threads
    with progress.phase(f'[{sample_id}] Classifying learned repeat boundaries'):
        updates = (bounded_ordered_map(locus_executor, _reclassify_locus, reclassify, resolve_threads(threads))
                   if locus_executor else map(_reclassify_locus, reclassify))
        evidence.update(updates)
    return learned, metadata
