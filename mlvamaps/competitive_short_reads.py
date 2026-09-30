"""Reference recruitment and competition over physical repeat lengths.

Reference targets recruit loci; native flank alignment measures reads. Once a read is
anchored, a one-sided repeat observation is censored length evidence, not an
alignment preference for a shorter reference. Paired fragments compete over
lengths using the measured library distribution. Neither operation needs an
allele database or a minimum molecule count.
"""
from collections import Counter
import math
import time

import numpy as np


def _fit_reference_locus(task):
    """Fit one locus's tied backgrounds inside the existing bounded worker pool."""
    import os
    from .locus_reconstruction import _fit_locus, _common_call
    from .short_read_evidence import classify_pair, flank_insert_length, estimate_insert_distribution
    name, targets, pairs, options, round_tolerance = task
    started = time.perf_counter()
    recoveries, distributions = [], {}
    for key, template in targets:
        items, lengths = [], []
        for pair in pairs:
            item = classify_pair(pair, template)
            if item is not None:
                items.append(item)
                length = flank_insert_length(pair, template)
                if length is not None:
                    lengths.append(length)
        fitting = dict(options)
        if fitting['insert'] is None:
            fitting['insert'] = estimate_insert_distribution(lengths)
        distributions[key] = fitting['insert']
        result, _, _ = _fit_locus((template, items, fitting, round_tolerance, None))
        call = _common_call(template.locus, result.products, len(items), options['sample_id'],
                            'illumina', 0, options['minimum_probability'], result)
        recoveries.append((key, template, items, result, call))
    return name, recoveries, distributions, {'seconds': time.perf_counter()-started,
        'worker_pid': os.getpid(), 'molecules': len(pairs), 'backgrounds': len(targets)}


def call_reference_loci(fitted, loci, templates, evidence, replay, output, progress,
                        *, database_path, options, round_tolerance, threads=1,
                        minimap2_bin='minimap2', keep_alignments=False,
                        reads1=None, reads2=None, orphan_path=None, executor=None):
    """Call loci from reads competitively recruited against reference targets.

    All loci compete for recruitment. Lengths require observed spans,
    reconstruction or fragment data, rather than a stored reference allele.
    """
    import filecmp
    import json
    from functools import partial
    from collections import defaultdict
    from contextlib import closing, nullcontext
    from pathlib import Path
    from tempfile import TemporaryDirectory
    from .candidate_contexts import generate_candidate_contexts, write_candidate_contexts, _repeat_template, _database_root
    from .io import normalize_read_id, write_tsv
    from .concurrency import bounded_ordered_map, resolve_threads
    from .minimap_mapping import map_reads_to_candidates_bam
    from .short_read_evidence import RepeatTemplate, _primitive_motif

    pending = {locus.locus_id: i for i, locus in enumerate(loci)}
    if not pending:
        return {}
    if not Path(database_path).is_dir():
        raise ValueError(f'SR FASTQ reference database does not exist: {database_path}')
    by_locus = {locus.locus_id: locus for locus in loci}
    contexts = generate_candidate_contexts(loci, database_path, maximum=options['maximum'])
    # Collapse allele expansions and duplicate references into flank backgrounds.
    backgrounds, context_keys, references = {}, {}, defaultdict(set)
    for context in contexts:
        locus = by_locus[context.locus_id]
        motif = _primitive_motif(_repeat_template(context, locus))
        if not motif or not context.repeat_unit_length:
            continue
        left, right = context.sequence[:context.repeat_start], context.sequence[context.repeat_end:]
        key = (context.locus_id, left, right, motif, context.repeat_unit_length)
        backgrounds[key] = RepeatTemplate(locus, left, right, motif, context.repeat_unit_length)
        context_keys[context.candidate_id] = key
        references[key].update(filter(None, (context.reference_id or context.reference_accession).split(';')))
    work = Path(output)/'reference_calling'
    written = write_candidate_contexts(contexts, work)
    paths = {f'reference_calling_{key}': value for key, value in written.items()}
    root = _database_root(database_path)
    resource = root/'competitive_mapping' if (root/'competitive_mapping').is_dir() else root
    cached_index, cached_fasta = resource/'short.mmi', resource/'candidate_contexts.fasta'
    reuse_index = (cached_index.is_file() and cached_fasta.is_file()
                   and filecmp.cmp(written['fasta'], cached_fasta, shallow=False))
    mapping_reference = cached_index if reuse_index else written['fasta']
    performance = {'candidate_contexts': len(contexts), 'reference_backgrounds': len(backgrounds),
                   'cached_index': reuse_index, 'mapping_runs': [], 'stage_seconds': {}, 'loci': {}}
    progress.step(f"[{options['sample_id']}] Calling {len(pending)} loci against reference targets")
    # One combined replay preserves QC, orphan mates and iterator-only callers.
    with TemporaryDirectory(prefix='mlvamaps-reference-') as temporary:
        first, second, single = (Path(temporary)/name for name in ('r1.fq', 'r2.fq', 'single.fq'))
        if reads1 is not None and Path(reads1).is_file():
            inputs = [(reads1, reads2)]
            orphan = Path(orphan_path) if orphan_path else Path(output)/'filtered_orphan_reads.fastq.gz'
            if orphan.is_file():
                inputs.append((orphan, None))
        else:
            paired = unpaired = 0
            iterator = replay()
            with first.open('w') as a, second.open('w') as b, single.open('w') as c, (closing(iterator) if hasattr(iterator, 'close') else nullcontext(iterator)):
                for pair in iterator:
                    if pair.read2:
                        paired += 1
                        records = ((pair.read1, a), (pair.read2, b))
                    else:
                        unpaired += 1
                        records = ((pair.read1, c),)
                    for read, handle in records:
                        handle.write(f'@{read.read_id}\n{read.sequence}\n+\n{read.quality or "I"*len(read.sequence)}\n')
            inputs = ([(first, second)] if paired else []) + ([(single, None)] if unpaired else [])
        grouped = {}
        for index, (first, second) in enumerate(inputs):
            bam = work/f'alignments_{index}.bam'
            stats = {}
            started = time.perf_counter()
            with progress.phase(f"[{options['sample_id']}] Reference mapping input {index+1}/{len(inputs)}",
                                f"{len(contexts):,} contexts; {threads} threads; cached index={reuse_index}"):
                reduced = map_reads_to_candidates_bam(mapping_reference, first, second, contexts,
                    bam, threads, 'illumina', executable=minimap2_bin,
                    retain_all_competitors=True, max_secondary=max(100, len(contexts)),
                    score_groups=context_keys, keep_alignments=keep_alignments, statistics=stats,
                    on_progress=partial(progress.count, f"[{options['sample_id']}] Reference alignment records"))
            stats['seconds'] = time.perf_counter()-started
            performance['mapping_runs'].append(stats)
            for molecule, matches in reduced.items():
                if molecule not in grouped:
                    grouped[molecule] = matches
                    continue
                for key, mates in matches.items():
                    target = grouped[molecule].setdefault(key, {})
                    for mate, score in mates.items():
                        target[mate] = max(target.get(mate, -math.inf), score)
            del reduced
            if keep_alignments:
                paths[f'reference_calling_alignments_{index}'] = bam
        performance['stage_seconds']['mapping'] = sum(run['seconds'] for run in performance['mapping_runs'])

    started = time.perf_counter()
    assignments, support = {}, defaultdict(lambda: defaultdict(float))
    for molecule, matches in grouped.items():
        scores = defaultdict(float)
        for key, mates in matches.items():
            scores[key[0]] = max(scores[key[0]], sum(mates.values()))
        ranked = sorted(scores, key=lambda name: (-scores[name], name))
        winner = ranked[0]
        if len(ranked) > 1 and scores[winner]-scores[ranked[1]] < max(12, .2*scores[winner]):
            continue
        if winner in pending:
            assignments[molecule] = winner
            for key, mates in matches.items():
                if key[0] == winner:
                    support[winner][key] += sum(mates.values())
    del grouped
    performance['stage_seconds']['competition'] = time.perf_counter()-started
    started = time.perf_counter()
    progress.step(f"[{options['sample_id']}] Recovering {len(assignments):,} assigned molecules from QC replay")
    pools = {name: {} for name in pending}
    iterator = replay()
    with closing(iterator) if hasattr(iterator, 'close') else nullcontext(iterator):
        for scanned, pair in enumerate(iterator, 1):
            if scanned % 100000 == 0:
                progress.count(f"[{options['sample_id']}] QC replay molecules scanned", scanned)
            name = assignments.get(normalize_read_id(pair.molecule_id)[0])
            if name is not None:
                pools[name][normalize_read_id(pair.molecule_id)[0]] = pair
    performance['stage_seconds']['read_replay'] = time.perf_counter()-started
    progress.step(f"[{options['sample_id']}] QC replay finished in {performance['stage_seconds']['read_replay']:.1f}s")
    diagnostics, targets = {}, {}
    for name in pending:
        row = {'sample_id': options['sample_id'], 'locus_id': name,
               'status': 'no_reference_support', 'reference_ids': '',
               'molecules': 0, 'method': '', 'repeat_count': '', 'product_size_bp': '',
               'insert_mean': '', 'insert_sd': '', 'insert_source': ''}
        scores = support[name]
        best_score = max(scores.values(), default=-math.inf)
        keys = sorted(key for key, score in scores.items() if score == best_score)
        # ponytail: evaluate at most eight equally supported backgrounds. More
        # remain unresolved; a larger reference panel needs a batched fitter.
        if len(keys) > 8:
            row['status'] = 'reference_background_limit'
            keys = []
        diagnostics[name] = row
        if keys:
            targets[name] = [(key, backgrounds[key]) for key in keys]
    tasks = ((name, target, list(pools[name].values()), options, round_tolerance)
             for name, target in targets.items())
    workers = min(resolve_threads(threads), len(targets)) if executor is not None else 1
    results = (bounded_ordered_map(executor, _fit_reference_locus, tasks, max(1, workers))
               if executor is not None else map(_fit_reference_locus, tasks))
    started = time.perf_counter()
    progress.step(f"[{options['sample_id']}] Fitting {len(targets)} reference loci with up to {workers} workers")
    for name, recoveries, distributions, stats in results:
        performance['loci'][name] = stats
        progress.step(f"[{options['sample_id']}] Reference locus {name} finished in {stats['seconds']:.1f}s; "
                      f"{stats['molecules']:,} molecules, {stats['backgrounds']} backgrounds")
        row, index = diagnostics[name], pending[name]
        if recoveries:
            calls = [entry[-1] for entry in recoveries]
            refs = ';'.join(sorted(set().union(*(references[entry[0]] for entry in recoveries))))
            row.update(reference_ids=refs, molecules=max(len(entry[2]) for entry in recoveries))
            # A tied background must agree on an identifiable physical length.
            if (len(recoveries) == 1 or all(call['status'] == 'called' for call in calls)
                    and len({(call['repeat_count'], call['product_size_bp']) for call in calls}) == 1):
                key, template, items, result, call = recoveries[0]
                result.reason = 'reference_assisted:' + refs + ('; '+result.reason if result.reason else '')
                for product in result.products:
                    product.evidence.update(reference_assisted=True, reference_ids=refs)
                fitted[index], templates[name], evidence[name] = result, template, items
                row.update(status=call['status'], reference_ids=refs, molecules=len(items),
                    method=result.method, repeat_count=call['repeat_count'], product_size_bp=call['product_size_bp'])
                if distributions[key] is not None:
                    distribution = distributions[key]
                    row.update(insert_mean=distribution.mean, insert_sd=distribution.sd,
                               insert_source=distribution.source)
            else:
                row['status'] = 'reference_background_ambiguous'
                evidence[name] = recoveries[0][2]
                fitted[index].reason = 'reference_background_ambiguous'
    performance['stage_seconds']['locus_fitting'] = time.perf_counter()-started
    performance['workers_used'] = len({row['worker_pid'] for row in performance['loci'].values()})
    diagnostics = list(diagnostics.values())
    paths['reference_calling_performance'] = work/'performance.json'
    paths['reference_calling_performance'].write_text(json.dumps(performance, indent=2)+'\n')
    paths['reference_calling'] = work/'summary.tsv'
    write_tsv(diagnostics, paths['reference_calling'], list(diagnostics[0]))
    return paths


def candidate_likelihood(items, template, insert, maximum, minimum_probability):
    """Haploid fragment likelihood with phase-corrected, censored boundaries.

    Compute geometry once per molecule and score all lengths as NumPy vectors.
    A boundary-only read has equal likelihood for every compatible length;
    sequencing depth cannot identify its unobserved end. Identical fragment
    geometries share computation while retaining their molecule multiplicity.
    """
    from .targeted_reconstruction import LocusRecovery

    if maximum < 1 or maximum * template.unit > 1_000_000:
        raise ValueError('invalid candidate ceiling or candidate locus exceeds sequence safety limit')
    nonrepeat = len(template.left)+len(template.right)
    lower = 0.
    offsets = Counter()
    for item in items:
        if 'discordant' in item.classes:
            continue
        context = item.template or template
        phase = len(context.left)+len(context.right)-nonrepeat
        if item.lower_bound:
            lower = max(lower, item.lower_bound+phase/template.unit)
        if insert and 'FLANK_PAIR' in item.classes and item.fragment_offset is not None:
            offsets[item.fragment_offset-phase] += 1
    upper = min(maximum, max(10, template.locus.expected_max_repeats+2, math.ceil(lower)+3))
    if insert and offsets:
        upper = min(maximum, max(upper, math.ceil((insert.mean+4*insert.sd-min(offsets))/template.unit)))
    if lower and not offsets:
        upper = maximum  # A censored observation supplies no finite upper bound.
    step = .5 if template.unit > 1 else 1.
    while True:
        states = np.arange(0, upper+step/2, step)
        joint = np.zeros(len(states))
        for offset, count in offsets.items():
            # Student-t tails tolerate fragment outliers without an absolute
            # support threshold. Library uncertainty is retained at low depth.
            fragment = offset+states*template.unit
            joint -= count*2.5*np.log1p(((fragment-insert.mean)/insert.sd)**2/4)
        if lower:
            # This is a censored observation. Repeating the same boundary must
            # not concentrate probability at the shortest compatible allele.
            joint[states < lower-.25] = -np.inf
        finite = np.isfinite(joint)
        if not finite.any():
            return LocusRecovery(method='AMBIGUOUS', limit_reached=True,
                reason='candidate_limit_below_observed_lower_bound')
        posterior = np.exp(np.maximum(joint-joint[finite].max(), -700))
        posterior[~finite] = 0
        posterior /= posterior.sum()
        if upper >= maximum or posterior[-2:].sum() < .01:
            break
        upper = min(maximum, upper*2)
    ranking = np.argsort(-posterior, kind='stable')
    best, second = map(int, ranking[:2])
    selected = ranking[:np.searchsorted(np.cumsum(posterior[ranking]), .95)+1]
    limit = bool(lower and not offsets) or upper == maximum and posterior[-2:].sum() >= .01
    implied = Counter()
    for offset, count in offsets.items():
        implied[round((insert.mean-offset)/template.unit)] += count
    modes = implied.most_common()
    mixed = bool(len(modes)>1 and modes[1][1] >= .2*sum(implied.values())
                 and abs(modes[0][0]-modes[1][0])*template.unit > 4*insert.sd)
    confidence = float(posterior[best])
    informative = bool(lower or offsets)
    tied = np.count_nonzero(np.isclose(joint, joint[best], rtol=0, atol=1e-8)) > 1
    identifiable = bool(offsets) and confidence >= minimum_probability and not tied and not limit and not mixed
    return LocusRecovery(method='MIXED' if mixed else 'INFERRED' if identifiable else 'AMBIGUOUS',
        confidence=confidence if informative else 0, states=states, log_likelihoods=joint, posterior=posterior,
        best=float(states[best]) if offsets and not tied else None,
        second=float(states[second]) if offsets and not tied else None,
        interval=((float(states[finite].min()), float(maximum)) if lower and not offsets else
                  (float(states[selected].min()), float(states[selected].max()))) if informative else None,
        identifiable=identifiable, limit_reached=limit,
        reason='no_repeat_length_information' if not informative else
               'repeat_length_lower_bound' if not offsets else
               'minimum_likelihood_tie' if tied else
               'candidate_limit_reached' if limit else 'competitive_fragment_likelihood')
