"""Bounded, sample-derived recruitment and repeat-graph construction.

No reference sequences or nominal repeat counts are used to learn anchors.
The graph is represented by two observed arms and a variable repeat edge;
targeted_reconstruction scores its bounded paths with native alignments.
"""
from collections import Counter, defaultdict
from dataclasses import replace
from functools import lru_cache
from os.path import commonprefix

import regex

from .models import ReadPair, ReadRecord
from .sequence import revcomp
from .short_read_evidence import MoleculeEvidence, RepeatTemplate, classify_pair, primer_bounds
from .targeted_reconstruction import _overlaps, _pileup_sequences


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
        for seq in pool:
            bounds = primer_bounds(seq, primer)
            if bounds:
                starts.append(seq[bounds[0]:])
        if not starts:
            return ''
        arm = commonprefix(starts)
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
    right = extend([revcomp(s) for s in sequences], template.locus.reverse_primer, revcomp(motif))
    return left, revcomp(right)


def learn_template(items, template):
    """Require independently supported motif and both observed primer arms."""
    if not template.primer_only or not template.unit:
        return template, {'source': 'panel', 'motifs': []}
    sequences = reliable_sequences(items)
    votes = defaultdict(set)
    for seq, members in sequences.items():
        for start, end, motif in repeat_runs(seq, template.unit):
            if end-start >= max(12, 2*template.unit):
                votes[motif].update(members)
    ranked = sorted(votes, key=lambda motif: (-len(votes[motif]), motif))
    info = {'source': 'sample', 'motifs': [{'sequence': m, 'molecules': len(votes[m])} for m in ranked[:3]],
            'graph_ready': False}
    if not ranked or len(votes[ranked[0]]) < 2:
        return template, info
    motif = ranked[0]
    if len(ranked) > 1 and len(votes[ranked[1]]) >= .5*len(votes[motif]):
        info['reason'] = 'ambiguous_sample_motif'
        return template, info
    # Deterministic depth-first ordering bounds noisy pools without selecting
    # a nominal allele. A capped pool may learn less sequence, never more.
    pool = sorted(sequences, key=lambda s: (-len(sequences[s]), -len(s), s))[:256]
    left_arm, right_arm = _arms(pool, template, motif)
    left, right = set(), set()
    for start, end, found in repeat_runs(left_arm, template.unit):
        if found == motif:
            position = left_arm.find(motif, start, end)
            if position >= len(template.locus.forward_primer):
                left.add(left_arm[:position])
    for start, end, found in repeat_runs(right_arm, template.unit):
        if found == motif:
            position = right_arm.rfind(motif, start, end)
            if position >= 0 and len(right_arm)-position-template.unit >= len(template.locus.reverse_primer):
                right.add(right_arm[position+template.unit:])
    if len(left) != 1 or len(right) != 1:
        info['reason'] = 'incomplete_sample_repeat_boundaries'
        return replace(template, motif=motif), info
    learned = RepeatTemplate(template.locus, left.pop(), right.pop(), motif, template.unit)
    info.update(graph_ready=True, left_bp=len(learned.left), right_bp=len(learned.right), motif=motif,
                left_sequence=learned.left, right_sequence=learned.right,
                structure=f'{learned.left}({motif})*{learned.right}')
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
    def __init__(self, evidence, templates):
        self.templates = templates
        self.index = defaultdict(dict)
        self.capped_loci = []
        for locus_id, items in sorted(evidence.items()):
            template = templates[locus_id]
            sequences = reliable_sequences(items)
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


def recruit_sample_reads(evidence, templates, replay, progress, sample_id, max_rounds=2,
                        threads=1, audit_writer=None):
    """Freeze the index per round; skip pairs already assigned to any locus."""
    assigned = {item.molecule_id for items in evidence.values() for item in items}
    learned = dict(templates)
    metadata = {'rounds': [], 'loci': {}}
    previous_index = None
    for round_index in range(max_rounds):
        for locus_id, template in templates.items():
            with progress.phase(f'[{sample_id}] Sample anchor learning {locus_id}',
                                f"round {round_index+1}; {len(evidence[locus_id]):,} molecules"):
                learned[locus_id], metadata['loci'][locus_id] = learn_template(evidence[locus_id], template)
        with progress.phase(f'[{sample_id}] Sample anchor index', f'round {round_index+1}'):
            recruiter = SampleRecruiter(evidence, learned)
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
    for locus_id, template in templates.items():
        with progress.phase(f'[{sample_id}] Final sample template {locus_id}',
                            f'{len(evidence[locus_id]):,} molecules'):
            learned[locus_id], metadata['loci'][locus_id] = learn_template(evidence[locus_id], template)
            if learned[locus_id] != template:
                updated = []
                for item in evidence[locus_id]:
                    classified = classify_pair(evidence_pair(item), learned[locus_id])
                    updated.append(classified or replace(item, template=learned[locus_id]))
                evidence[locus_id] = updated
    return learned, metadata
