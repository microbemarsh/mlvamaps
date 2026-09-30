"""Read-derived lengths for loci whose overlap assembly cannot count cycles.

A weighted de Bruijn graph collapses repeated sequence into shared k-mers.
Coverage of those k-mers relative to the two single-copy primer arms estimates
their multiplicity. Count coverage in the *whole input*, not the anchor-selected
recruitment pool, which systematically excludes internal repeat reads.
"""
from collections import Counter, defaultdict, deque
from functools import lru_cache
from itertools import islice
import math
import time
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import regex

from .concurrency import bounded_ordered_map, resolve_threads
from .sequence import revcomp
from .short_read_evidence import primer_bounds
from .repeat_calibration import assembly_equivalent_product_allele, expected_nonrepeat_bp
from .targeted_reconstruction import _pileup_sequences


def _segments(sequence, quality, k):
    if quality is None:
        yield from (m.group() for m in regex.finditer('[ACGT]{%d,}' % k, sequence))
    else:
        for match in regex.finditer('[5-~]{%d,}' % k, quality):
            yield from _segments(sequence[match.start():match.end()], None, k)


def _graph(items, template, k=21):
    """Keep observed edges connecting the primers, without a depth cutoff."""
    observations = defaultdict(set)
    sequences, starts, ends = Counter(), Counter(), Counter()
    left_positions, right_positions = set(), set()
    for item in items:
        if 'discordant' in item.classes:
            continue
        for i, sequence in enumerate(item.sequences):
            quality = item.qualities[i] if item.qualities else None
            for segment in _segments(sequence, quality, k):
                observations[segment].add(item.molecule_id)
    # Reuse coordinate-based error consolidation; singleton observations remain.
    for sequence, members in _pileup_sequences(observations).items():
        for segment in _segments(sequence, None, k):
            sequences[segment] += len(members)
    for sequence, weight in sequences.items():
        left = primer_bounds(sequence, template.locus.forward_primer)
        right = primer_bounds(sequence, revcomp(template.locus.reverse_primer))
        if left and left[0]+k-1 <= len(sequence):
            starts[sequence[left[0]:left[0]+k-1]] += weight
            left_positions.add(left[0])
        if right and right[1] >= k-1:
            ends[sequence[right[1]-k+1:right[1]]] += weight
            right_positions.add(len(sequence)-right[1])
    if not starts or not ends:
        return None
    start = min(starts, key=lambda n: (-starts[n], n))
    end = min(ends, key=lambda n: (-ends[n], n))
    if start == end:
        return None
    counts = Counter()
    for sequence, weight in sequences.items():
        for i in range(len(sequence)-k+1):
            counts[sequence[i:i+k]] += weight
    outgoing, incoming = defaultdict(list), defaultdict(list)
    for edge in counts:
        outgoing[edge[:-1]].append(edge)
        incoming[edge[1:]].append(edge)

    def reachable(seed, adjacency, backward=False):
        seen, pending = {seed}, [seed]
        while pending:
            node = pending.pop()
            for edge in adjacency.get(node, ()):
                other = edge[:-1] if backward else edge[1:]
                if other not in seen:
                    seen.add(other)
                    pending.append(other)
        return seen

    forward, backward = reachable(start, outgoing), reachable(end, incoming, True)
    if end not in forward:
        return None
    edges = {e for node in forward for e in outgoing.get(node, ()) if e[1:] in backward}
    outgoing, incoming = defaultdict(list), defaultdict(list)
    for edge in sorted(edges):
        outgoing[edge[:-1]].append(edge)
        incoming[edge[1:]].append(edge)
    # Stop at the first branch/merge on either arm: repeat-cycle edges have
    # elevated multiplicity and must never define single-copy depth.
    arms = []
    for seed, adjacency, reverse_adjacency, reverse in (
            (start, outgoing, incoming, False), (end, incoming, outgoing, True)):
        node, arm, seen = seed, set(), set()
        while node not in seen and len(adjacency.get(node, ())) == 1:
            if node != seed and len(reverse_adjacency.get(node, ())) != 1:
                break
            seen.add(node)
            edge = adjacency[node][0]
            arm.add(edge)
            node = edge[:-1] if reverse else edge[1:]
        arms.append(arm)
    if not all(arms):
        return None
    distances, pending = {start: 0}, deque([start])
    while pending:
        node = pending.popleft()
        for edge in outgoing.get(node, ()):
            if edge[1:] not in distances:
                distances[edge[1:]] = distances[node]+1
                pending.append(edge[1:])
    return {'k': k, 'edges': edges, 'arms': arms, 'minimum_bp': distances[end]+k-1,
            'read_start_positions': [len(left_positions), len(right_positions)]}


def _locus_graph(task):
    name, items, template = task
    graph = _graph(items, template)
    return name, graph if graph is not None else _graph(items, template, k=15)


_DNA_CODES = bytes('ACGT'.find(chr(i)) if chr(i) in 'ACGT' else 4 for i in range(256))


def _encoded_kmers(sequence, k):
    """Exact two-bit keys for our 15/21-mers; exclude ambiguous/boundary bases."""
    bases = np.frombuffer(sequence.translate(_DNA_CODES), dtype=np.uint8)
    size = len(bases)-k+1
    if size <= 0:
        return np.array([], dtype=np.uint64)
    codes = np.zeros(size, dtype=np.uint64)
    for offset in range(k):
        np.left_shift(codes, 2, out=codes)
        np.bitwise_or(codes, bases[offset:offset+size], out=codes)
    invalid = np.concatenate(([0], np.cumsum(bases > 3)))
    return codes[invalid[k:] == invalid[:-k]]


def _depth_targets(keys, k):
    # A/C/G/T lexical order is also two-bit numeric order.
    bases = np.frombuffer(''.join(keys).encode('ascii').translate(_DNA_CODES),
                          dtype=np.uint8).reshape(-1, k)
    codes = np.zeros(len(keys), dtype=np.uint64)
    for offset in range(k):
        np.left_shift(codes, 2, out=codes)
        np.bitwise_or(codes, bases[:, offset], out=codes)
    return codes


@lru_cache(maxsize=2)
def _depth_index(paths):
    # Unique temporary paths isolate samples and rounds. Each process attaches
    # once; the OS shares read-only pages instead of copying the index per task.
    return {k: np.load(path, mmap_mode='r', allow_pickle=False) for k, path in paths}


def _count_depth_chunk(task):
    pairs, paths = task
    targets_by_k = _depth_index(paths)
    minimum = min(targets_by_k)
    sequence = b'N'.join(segment.encode('ascii') for pair in pairs
        for read in (pair.read1, pair.read2) if read is not None
        for segment in _segments(read.sequence.upper(), read.quality, minimum))
    counts = {}
    for k, targets in targets_by_k.items():
        codes = _encoded_kmers(sequence, k)
        offsets = np.searchsorted(targets, codes)
        np.minimum(offsets, len(targets)-1, out=offsets)
        matched = offsets[targets[offsets] == codes]
        # Work scales with this batch's hits, not the entire graph index.
        counts[k] = np.unique(matched, return_counts=True)
    return len(pairs), counts


def estimate_graph_lengths(recoveries, loci, templates, evidence, replay_pairs,
                           round_tolerance=.25, progress=None, sample_id='',
                           *, locus_executor=None, threads=1):
    """Populate provisional counts/lengths; never manufacture a contig or SNPs."""
    graphs = {}
    diagnostics = {}
    tasks = ((locus.locus_id, evidence[locus.locus_id], templates[locus.locus_id])
             for locus, recovery in zip(loci, recoveries)
             if not (recovery.products or recovery.identifiable or recovery.method == 'MIXED')
             and evidence[locus.locus_id] and templates[locus.locus_id].unit > 0)
    results = (bounded_ordered_map(locus_executor, _locus_graph, tasks, resolve_threads(threads))
               if locus_executor else map(_locus_graph, tasks))
    for name, graph in results:
        if graph is not None:
            graphs[name] = graph
        else:
            diagnostics[name] = {'method': 'unresolved', 'reason': 'no_primer_connected_kmer_graph'}
    if not graphs:
        return diagnostics
    # Both strand spellings share one counter. Use one combined scan for all
    # unresolved loci, and retain no background reads or background k-mers.
    wanted = defaultdict(set)
    for graph in graphs.values():
        for edge in graph['edges']:
            wanted[graph['k']].update((edge, revcomp(edge)))
    wanted = {k: tuple(sorted(keys)) for k, keys in wanted.items()}
    counts = Counter()
    if progress:
        progress.step(f'[{sample_id}] Estimating {len(graphs)} unresolved locus lengths from whole-input k-mer depth; '
                      f'{sum(map(len, wanted.values())):,} target k-mers in a shared native index')
    started = time.perf_counter()
    with TemporaryDirectory(prefix='mlvamaps-depth-') as temporary:
        paths = tuple((k, str(Path(temporary)/f'{k}.npy')) for k in wanted)
        for k, path in paths:
            np.save(path, _depth_targets(wanted[k], k), allow_pickle=False)
        iterator = iter(replay_pairs())
        def chunks():
            while chunk := list(islice(iterator, 1024)):
                yield chunk, paths
        try:
            results = (bounded_ordered_map(locus_executor, _count_depth_chunk, chunks(), resolve_threads(threads))
                       if locus_executor else map(_count_depth_chunk, chunks()))
            examined = 0
            for size, hits in results:
                examined += size
                for k, (positions, multiplicities) in hits.items():
                    counts.update({wanted[k][i]: int(n) for i, n in zip(positions, multiplicities)})
                if progress:
                    progress.count(f'[{sample_id}] Read pairs counted for repeat depth', examined,
                                   detail=f'{examined/max(time.perf_counter()-started, 1e-9):,.0f} pairs/s')
        finally:
            if hasattr(iterator, 'close'):
                iterator.close()
            _depth_index.cache_clear()
    owners = Counter(e for g in graphs.values() for e in {min(e, revcomp(e)) for e in g['edges']})
    for locus, recovery in zip(loci, recoveries):
        name = locus.locus_id
        if name not in graphs:
            continue
        graph = graphs[name]
        coverage = {e: counts[e]+(counts[revcomp(e)] if e != revcomp(e) else 0) for e in graph['edges']}
        arm_depths = [float(np.median([coverage[e] for e in arm])) for arm in graph['arms']]
        depth = float(np.mean(arm_depths))
        if depth <= 0:
            diagnostics[name] = {'method': 'unresolved', 'reason': 'no_single_copy_depth'}
            continue
        total = sum(coverage.get(e, coverage.get(revcomp(e), 0))
                    for e in {min(e, revcomp(e)) for e in graph['edges']})
        length = max(graph['minimum_bp'], round(graph['k']-1+total/depth))
        # ponytail: this is a uniform-local-coverage point estimate, not a
        # calibrated genotype posterior. Shared sequence and library bias can
        # inflate it. Keep a conservative sensitivity interval and expose the
        # assumptions; upgrade to a fitted coverage model with validated data.
        variation = max(.1, abs(arm_depths[0]-arm_depths[1])/depth, 2/math.sqrt(depth))
        bounds = (max(graph['minimum_bp'], round(length/(1+variation))),
                  round(length*(1+variation)))
        template = templates[name]
        calibrated = expected_nonrepeat_bp(locus) is not None or bool(
            locus.left_flank_sequence and locus.right_flank_sequence)
        if calibrated:
            best = assembly_equivalent_product_allele(locus, length, round_tolerance)[1]
            # The historical absolute-value calibration is nonmonotone below
            # the nonrepeat length; include zero when the interval crosses it.
            values = [assembly_equivalent_product_allele(locus, b, round_tolerance)[1] for b in bounds]
            nonrepeat = expected_nonrepeat_bp(locus)
            if nonrepeat is not None and bounds[0] <= nonrepeat <= bounds[1]:
                values.append(0)
            interval = min(values), max(values)
        elif not template.primer_only:
            best = max(0, (length-len(template.left)-len(template.right))/template.unit)
            interval = tuple(max(0, (b-len(template.left)-len(template.right))/template.unit) for b in bounds)
        else:
            best, interval = None, None
        if interval is not None and recovery.interval is not None:
            interval = min(interval[0], recovery.interval[0]), max(interval[1], recovery.interval[1])
        shared = sum(owners[min(e, revcomp(e))] > 1 for e in graph['edges'])
        diagnostics[name] = {'method': 'KMER_DEPTH', 'k': graph['k'], 'graph_edges': len(graph['edges']),
            'arm_depths': arm_depths, 'read_start_positions': graph['read_start_positions'],
            'amplicon_length': length, 'length_interval': bounds,
            'shared_edges': shared, 'interval_kind': 'coverage_sensitivity'}
        recovery.product_size_bp = length
        recovery.best, recovery.interval = best, interval
        recovery.states = np.array([])
        recovery.log_likelihoods = np.array([])
        recovery.posterior = np.array([])
        recovery.second = None
        recovery.method, recovery.confidence, recovery.identifiable = 'KMER_DEPTH', 0.0, False
        recovery.reason = '; '.join(filter(None, (recovery.reason, 'read_depth_length_estimate',
            'shared_repeat_sequence' if shared else '',
            'repeat_count_uncalibrated' if best is None else '', 'interval_is_coverage_sensitivity')))
    return diagnostics
