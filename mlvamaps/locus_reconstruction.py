"""Input-specific recovery into canonical products; no known-allele selection."""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import replace
import csv
import json
import time
from pathlib import Path

import numpy as np
import parasail

from .concurrency import resolve_threads
from .progress import ProgressReporter
from .io import read_fastq, read_fastq_pairs, write_fasta, write_tsv
from .locus_measurement import find_anchor
from .locus_products import LocusProduct, genotype_product, write_products
from .mixture import estimate_variant_mixtures
from .models import Locus, ReadPair
from .repeat_likelihood import (
    _read_profile, RepeatTemplate, estimate_insert_distribution, infer_repeat, panel_template,
)
from .sequence import revcomp


def locus_templates(loci, database_path=None):
    """Prefer rich panels; borrow sequence context only for incomplete panels."""
    templates = {l.locus_id: panel_template(l) for l in loci}
    missing = [l for l in loci if templates[l.locus_id] is None]
    if missing and database_path:
        from .candidate_contexts import _base_contexts, _repeat_template
        contexts = _base_contexts(missing, database_path)
        by_id = {l.locus_id: l for l in loci}
        for context in contexts:
            if templates[context.locus_id] is not None:
                continue
            locus = by_id[context.locus_id]
            motif = _repeat_template(context, locus)
            if motif and context.repeat_unit_length:
                templates[locus.locus_id] = RepeatTemplate(
                    locus, context.sequence[:context.repeat_start], context.sequence[context.repeat_end:],
                    motif, context.repeat_unit_length,
                )
    return templates


def recover_long_reads(pairs, loci, sample_id, min_fraction=0.01, min_secondary_reads=2,
                       max_anchor_edits=3):
    """Retain independently observed complete molecules, including novel states."""
    clusters = defaultdict(list)
    audit = []
    recruited = Counter()
    for pair in pairs:
        sequence = pair.read1.sequence.upper()
        hits = []
        for locus in loci:
            best = None
            for strand, oriented in (("+", sequence), ("-", revcomp(sequence))):
                fwd = find_anchor(locus.forward_primer, oriented, max_anchor_edits)
                rev = find_anchor(revcomp(locus.reverse_primer), oriented, max_anchor_edits,
                                  fwd.end if fwd else 0)
                score = sum(a.identity for a in (fwd, rev) if a)
                if score and (best is None or score > best[0]):
                    best = (score, strand, oriented, fwd, rev)
            if best:
                hits.append((locus, best))
        # Full primer pairs identify products independently. Single-anchor
        # ambiguous recruitment is retained in the audit but cannot genotype.
        for locus, (_, strand, oriented, fwd, rev) in hits:
            recruited[locus.locus_id] += 1
            complete = bool(fwd and rev and rev.start >= fwd.end)
            product = oriented[fwd.start:rev.end] if complete else ""
            calibrated_size = (len(product) - (fwd.end-fwd.start) - (rev.end-rev.start)
                               + len(locus.forward_primer) + len(locus.reverse_primer)) if complete else None
            count = genotype_product(LocusProduct(sample_id, locus.locus_id, "long_read", "molecule", product, calibrated_product_size_bp=calibrated_size), locus).repeat_count if complete else None
            row = {"sample_id": sample_id, "locus_id": locus.locus_id, "molecule_id": pair.molecule_id,
                   "orientation": strand, "complete": complete, "repeat_count": count,
                   "anchor_edits": sum(a.edit_distance for a in (fwd, rev) if a),
                   "forward_start": fwd.start if fwd else "", "reverse_end": rev.end if rev else "",
                   "variant_id": "", "calibrated_product_size_bp": calibrated_size}
            audit.append(row)
            if complete:
                clusters[(locus.locus_id, product)].append(row)
    rows = []
    for index, ((locus_id, sequence), members) in enumerate(sorted(clusters.items()), 1):
        variant = f"{locus_id}|v{index}"
        for member in members:
            member["variant_id"] = variant
        rows.append({"sample_id": sample_id, "locus_id": locus_id, "variant_id": variant,
                     "representative_sequence": sequence, "representative_length_bp": len(sequence),
                     "support_reads": len(members), "repeat_count": members[0]["repeat_count"]})
    mixtures = estimate_variant_mixtures(rows, min_fraction=min_fraction, min_secondary_reads=min_secondary_reads)
    mixture_by_id = {r["variant_id"]: r for r in mixtures}
    products = []
    for row in rows:
        mix = mixture_by_id[row["variant_id"]]
        members = clusters[(row["locus_id"], row["representative_sequence"])]
        products.append(LocusProduct(sample_id, row["locus_id"], "long_read", row["variant_id"],
            row["representative_sequence"], support_count=row["support_reads"],
            calibrated_product_size_bp=members[0]["calibrated_product_size_bp"],
            effective_depth=float(mix["estimated_reads"]), estimated_fraction=float(mix["estimated_fraction"]),
            evidence={"n_recruited": recruited[row["locus_id"]], "n_complete": sum(len(v) for (l, _), v in clusters.items() if l == row["locus_id"]),
                      "molecule_ids": [m["molecule_id"] for m in members],
                      "orientations": dict(Counter(m["orientation"] for m in members)),
                      "anchor_edits": sum(m["anchor_edits"] for m in members),
                      "meaningful": mix["meaningful"], "abundance_class": mix["abundance_class"]}))
    return products, audit, recruited


def _project_flanks(item, sequence, template, count, alignment_cache=None):
    """Per-molecule base/indel votes; overlapping mates never count twice."""
    repeat_start, repeat_end = len(template.left), len(template.left) + round(count * template.unit)
    votes = defaultdict(set)
    alignment_cache = {} if alignment_cache is None else alignment_cache
    for mate, read in enumerate(item.sequences):
        quality = item.qualities[mate] if item.qualities else None
        if read not in alignment_cache:
            result = parasail.sw_trace_striped_profile_32(_read_profile(read), sequence, 5, 1)
            traceback = result.traceback
            q, r = traceback.query, traceback.ref
            if len(alignment_cache) >= 4096:
                alignment_cache.clear()
            alignment_cache[read] = (q, r, result.end_ref + 1 - len(r.replace("-", "")),
                                     result.end_query + 1 - len(q.replace("-", "")))
        q, r, position, query_position = alignment_cache[read]
        insertion = ""
        for base, ref in zip(q, r):
            reliable = quality is None or base == "-" or ord(quality[query_position])-33 >= 20
            if base != "-":
                query_position += 1
            if ref == "-":
                insertion += base if reliable else "N"
                continue
            if reliable and (position < repeat_start or position >= repeat_end):
                votes[position].add((insertion, base))
            insertion = ""
            position += 1
    return {p: next(iter(values)) for p, values in votes.items() if len(values) == 1}


def reconstruct_short_products(evidence, template, inference, sample_id, min_fraction=.01, min_secondary_reads=2):
    products = []
    if not inference.identifiable:
        return products
    n = len(evidence)
    states = list(inference.fractions)
    columns = [int(np.argmin(abs(inference.states - r))) for r in states]
    scores = inference.molecule_likelihoods[:, columns]
    probs = np.exp(np.maximum(scores - scores.max(axis=1, keepdims=True), -700))
    probs *= np.asarray([inference.fractions[r] for r in states])
    probs /= probs.sum(axis=1, keepdims=True)
    for column, repeat in enumerate(states):
        fraction = inference.fractions[repeat]
        if fraction < 1e-6:
            continue
        indexes = [i for i in range(n) if probs[i, column] >= 0.8]
        members = [evidence[i] for i in indexes]
        if not members:
            continue
        synthetic = template.sequence(repeat)
        full_members = defaultdict(list)
        for item in members:
            if item.product_sequence:
                full_members[item.product_sequence].append(item)
        full = Counter({seq: len(group) for seq, group in full_members.items()})
        # Exact observed products retain linked SNP/indel haplotypes. Singleton
        # sequences remain trace diagnostics rather than confirmed mixtures.
        phased = [(seq, support, full_members[seq])
                  for seq, support in full.most_common()]
        haplotype_weights = {}
        if full:
            haplotype_rows = [{"sample_id": sample_id, "locus_id": template.locus.locus_id,
                "variant_id": str(i), "representative_sequence": seq,
                "representative_length_bp": len(seq), "support_reads": support,
                "repeat_count": repeat} for i, (seq, support, _) in enumerate(phased)]
            estimates = estimate_variant_mixtures(haplotype_rows, min_fraction=min_fraction,
                                                  min_secondary_reads=min_secondary_reads)
            haplotype_weights = {phased[int(row["variant_id"])][0]: float(row["estimated_fraction"])
                                 for row in estimates}
        if not phased:
            # Raw alignments can be reused across quality variants; full
            # projections also require identical qualities. Reset per candidate.
            alignment_cache = {}
            projection_cache = {}
            projections = []
            for item in members:
                key = (item.sequences, item.qualities)
                if key not in projection_cache:
                    if len(projection_cache) >= 4096:
                        projection_cache.clear()
                    projection_cache[key] = _project_flanks(item, synthetic, template, repeat, alignment_cache)
                # Identical sequence/quality pairs share immutable projections,
                # but still contribute one vote for every original molecule.
                projections.append(projection_cache[key])
            pileup = defaultdict(Counter)
            for projection in projections:
                for position, allele in projection.items():
                    pileup[position][allele] += 1
            heterozygous = {p for p, votes in pileup.items()
                            if sum(v >= max(2, 0.1 * sum(votes.values())) for v in votes.values()) > 1}
            # Phase only when individual molecules cover all variable sites.
            signatures = Counter(tuple((p, projection[p]) for p in sorted(heterozygous))
                                 for projection in projections if heterozygous <= projection.keys()) if heterozygous else Counter()
            groups = [(sig, support) for sig, support in signatures.items() if support >= min_secondary_reads] or [((), len(members))]
            for signature, support in groups:
                chosen = [i for i, projection in enumerate(projections)
                          if all(projection.get(p) == allele for p, allele in signature)]
                group_votes = defaultdict(Counter)
                for i in chosen:
                    for p, allele in projections[i].items():
                        group_votes[p][allele] += 1
                sequence = []
                start, end = len(template.left), len(template.left) + round(repeat * template.unit)
                for p, base in enumerate(synthetic):
                    if start <= p < end:
                        sequence.append(base)
                    elif group_votes[p]:
                        (insertion, called), depth = group_votes[p].most_common(1)[0]
                        sequence.append(insertion + ("" if called == "-" else called) if depth / sum(group_votes[p].values()) >= 0.7 else "N")
                    else:
                        sequence.append("N")  # Never impute an uncovered reference SNP.
                phased.append(("".join(sequence), support, [members[i] for i in chosen]))
        total = sum(support for _, support, _ in phased)
        for sequence, support, group in phased:
            estimated_fraction = fraction * haplotype_weights.get(sequence, support / total)
            products.append(LocusProduct(sample_id, template.locus.locus_id, "short_read",
                f"{template.locus.locus_id}|v{len(products)+1}", sequence,
                support_count=support, effective_depth=float(n * estimated_fraction),
                estimated_fraction=estimated_fraction, reconstruction_confidence=inference.confidence,
                repeat_likelihoods={float(r): float(p) for r, p in zip(inference.states, inference.posterior)},
                evidence={"inferred_repeat": repeat, "classes": inference.counts,
                          "molecule_ids": [e.molecule_id for e in group],
                          "meaningful": "yes" if support >= min_secondary_reads and estimated_fraction >= max(min_fraction, 1/(n+1)) else "no",
                          "uncovered_bases": sequence.count("N"),
                          "phase": "observed_product" if sequence in full else "molecule_supported_flank_consensus"}))
    return products


def _common_call(locus, products, recruited, sample_id, technology, minimum_molecules, minimum_probability, inference=None):
    ranked = sorted(products, key=lambda p: (-p.estimated_fraction, -p.support_count, p.variant_id))
    best = ranked[0] if ranked else None
    confidence = best.reconstruction_confidence if best else inference.confidence if inference else 0
    meaningful = [p for p in ranked if p.evidence.get("meaningful") == "yes"]
    status = "not_found" if not recruited else "detected_unresolved"
    if best and confidence >= minimum_probability:
        status = "low_coverage" if best.effective_depth < minimum_molecules else "mixed" if len(meaningful) > 1 else "called"
    repeat = genotype_product(best, locus).repeat_count if best else ""
    count_distribution = defaultdict(float)
    for product in ranked:
        count_distribution[genotype_product(product, locus).repeat_count] += product.estimated_fraction
    counts = inference.counts if inference else {"E": sum(p.support_count for p in ranked)}
    return {"sample": sample_id, "locus": locus.locus_id, "technology": technology,
            "repeat_count": repeat, "status": status, "best_probability": confidence,
            "second_best_probability": sorted(inference.posterior, reverse=True)[1] if inference else 0,
            "molecule_support": recruited, "direct_product_support": counts.get("E", 0),
            "full_span_support": counts.get("E", 0), "junction_support": counts.get("F", 0),
            "dominant_fraction": best.estimated_fraction if best else "",
            "secondary_repeat": genotype_product(ranked[1], locus).repeat_count if len(ranked)>1 else "",
            "secondary_fraction": ranked[1].estimated_fraction if len(ranked)>1 else "",
            "candidate_distribution": ";".join(f"{r}:{p:.6f}" for r,p in sorted(count_distribution.items(), key=lambda v:-v[1])),
            "confidence": confidence, "margin": confidence,
            "best_candidate_repeat": inference.best if inference else repeat,
            "dominant_repeat": repeat,
            "num_variants": len(ranked), "num_secondary": max(0, len(meaningful)-1),
            "inference_method": "repeat_likelihood" if technology == "illumina" else "spanning_molecules",
            "n_spanning": counts.get("S", 0), "repeat_count_min": inference.interval[0] if inference else repeat,
            "repeat_count_max": inference.interval[1] if inference else repeat,
            "second_best_repeat_count": inference.second if inference else ""}


MOLECULE_AUDIT_FIELDS = [
    "sample_id", "locus_id", "molecule_id", "classes", "observed_repeat", "lower_bound",
    "fragment_offset", "orientation", "alignment", "complete", "repeat_count",
    "anchor_edits", "forward_start", "reverse_end", "variant_id",
]


def run_reconstructed_fastq_inference(*, outdir, stream_molecule_evidence=False, **kwargs):
    """Optionally stream audit rows; file contents are identical in either mode."""
    if stream_molecule_evidence and kwargs.get("technology") == "illumina":
        output = Path(outdir)
        output.mkdir(parents=True, exist_ok=True)
        with (output / "molecule_candidate_evidence.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MOLECULE_AUDIT_FIELDS, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            return _run_reconstructed_fastq_inference(outdir=outdir, audit_writer=writer, audit_handle=handle, **kwargs)
    return _run_reconstructed_fastq_inference(outdir=outdir, **kwargs)


def _run_reconstructed_fastq_inference(*, reads1, reads2, loci, database_path, outdir, sample_id,
        technology, minimum_molecules, minimum_probability, maximum_candidate_repeat_count=100,
        insert_mean=None, insert_sd=None, orphan_path=None, min_fraction=.01, min_secondary_reads=2,
        max_anchor_edits=3, threads=1, minimum_spanning_pairs=2, round_tolerance=.25,
        show_progress=False, audit_writer=None, audit_handle=None, pairs=None,
        recruitment_threads=None, **_unused):
    progress = ProgressReporter(enabled=show_progress)
    stage_seconds, recruitment_stats, locus_stats = {}, {}, {}
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    if pairs is None:
        pairs = read_fastq_pairs(reads1, reads2)
    if orphan_path:
        from itertools import chain
        pairs = chain(pairs, (ReadPair(r.read_id, r) for r in read_fastq(orphan_path)))
    all_products, calls, summaries, likelihood_rows, audit = [], [], [], [], []
    insert = None
    if technology != "illumina":
        all_products, audit, recruited = recover_long_reads(pairs, loci, sample_id, min_fraction, min_secondary_reads, max_anchor_edits)
        all_products = [replace(p, round_tolerance=round_tolerance) for p in all_products]
        for locus in loci:
            products = [p for p in all_products if p.locus_id == locus.locus_id]
            calls.append(_common_call(locus, products, recruited[locus.locus_id], sample_id, technology, minimum_molecules, minimum_probability))
    else:
        templates = locus_templates(loci, database_path)
        evidence = defaultdict(list)
        insert_lengths = []
        from .short_read_recruitment import recruit_short_reads
        recruitment_threads = threads if recruitment_threads is None else recruitment_threads
        progress.step(f"[{sample_id}] Recruiting short-read molecules with up to {resolve_threads(recruitment_threads)} worker(s)")
        started = time.perf_counter()
        for chunk in recruit_short_reads(pairs, templates, sample_id, recruitment_threads, statistics=recruitment_stats,
                                        audit_fields=MOLECULE_AUDIT_FIELDS if audit_handle is not None else None):
            recruitment_stats["pairs_examined"] = recruitment_stats.get("pairs_examined", 0) + chunk.examined
            recruitment_stats["locus_tests"] = recruitment_stats.get("locus_tests", 0) + chunk.locus_tests
            for item in chunk.evidence:
                evidence[item.locus_id].append(item)
            if audit_handle is not None:
                audit_handle.write(chunk.audit_tsv)
            elif audit_writer is None:
                audit.extend(chunk.ambiguous)
            else:
                audit_writer.writerows(chunk.ambiguous)
            insert_lengths.extend(chunk.insert_lengths)
            progress.count(f"[{sample_id}] Read pairs recruited/scanned", recruitment_stats["pairs_examined"])
        stage_seconds["recruitment"] = time.perf_counter() - started
        progress.step(f"[{sample_id}] Recruitment finished in {stage_seconds['recruitment']:.1f}s; fitting {len(loci)} loci")
        insert = estimate_insert_distribution(insert_lengths, insert_mean, insert_sd)
        def fit_locus(locus):
            fit_started = time.perf_counter()
            items = [e for e in evidence[locus.locus_id] if e.classes and "discordant" not in e.classes]
            template = templates[locus.locus_id]
            result = infer_repeat(items, template, insert, maximum_candidate_repeat_count, minimum_probability) if items and template else None
            if result and not any("E" in e.classes for e in items) and sum("S" in e.classes for e in items) < minimum_spanning_pairs:
                result.identifiable = False
            elapsed = time.perf_counter() - fit_started
            locus_stats[locus.locus_id] = {"repeat_fitting_seconds": elapsed,
                "evidence_molecules": len(items), "candidate_states": len(result.states) if result else 0}
            progress.step(f"[{sample_id}] {locus.locus_id}: fitted {len(items):,} molecules in {elapsed:.1f}s")
            return result
        from concurrent.futures import ThreadPoolExecutor
        started = time.perf_counter()
        if resolve_threads(threads) > 1 and len(loci) > 1:
            with ThreadPoolExecutor(max_workers=min(len(loci), resolve_threads(threads))) as executor:
                fitted = list(executor.map(fit_locus, loci))
        else:
            fitted = [fit_locus(locus) for locus in loci]
        stage_seconds["repeat_fitting"] = time.perf_counter() - started
        progress.step(f"[{sample_id}] Repeat fitting finished in {stage_seconds['repeat_fitting']:.1f}s; writing molecule likelihoods")
        started = time.perf_counter()
        molecule_likelihood_path = output / "molecule_repeat_likelihoods.tsv"
        likelihood_records = 0
        with molecule_likelihood_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample_id", "locus_id", "molecule_id", "repeat_count", "log_likelihood"])
            for locus, inference in zip(loci, fitted):
                if inference:
                    items = [e for e in evidence[locus.locus_id] if e.classes and "discordant" not in e.classes]
                    for item, scores in zip(items, inference.molecule_likelihoods):
                        writer.writerows((sample_id, locus.locus_id, item.molecule_id, state, score)
                                         for state, score in zip(inference.states, scores))
                        likelihood_records += len(inference.states)
                        progress.count(f"[{sample_id}] Molecule likelihood rows written", likelihood_records)
        stage_seconds["molecule_likelihood_output"] = time.perf_counter() - started
        progress.step(f"[{sample_id}] Reconstructing locus sequences")
        started = time.perf_counter()
        for locus, inference in zip(loci, fitted):
            locus_started = time.perf_counter()
            items = evidence[locus.locus_id]
            usable = [e for e in items if e.classes and "discordant" not in e.classes]
            template = templates[locus.locus_id]
            products = reconstruct_short_products(usable, template, inference, sample_id, min_fraction, min_secondary_reads) if inference else []
            locus_stats[locus.locus_id]["sequence_reconstruction_seconds"] = time.perf_counter() - locus_started
            all_products.extend(products)
            calls.append(_common_call(locus, products, len(items), sample_id, technology, minimum_molecules, minimum_probability, inference))
            summaries.append({"sample_id": sample_id, "locus_id": locus.locus_id,
                "E": inference.counts.get("E", 0) if inference else 0,
                "S": inference.counts.get("S", 0) if inference else 0,
                "F": inference.counts.get("F", 0) if inference else 0,
                "FRR": inference.counts.get("FRR", 0) if inference else 0,
                "best_repeat": inference.best if inference else "", "second_repeat": inference.second if inference else "",
                "confidence": inference.confidence if inference else 0,
                "interval_min": inference.interval[0] if inference else "", "interval_max": inference.interval[1] if inference else "",
                "effective_depth": len(usable), "limit_reached": inference.limit_reached if inference else False,
                "identifiable": inference.identifiable if inference else False,
                "components": json.dumps(inference.fractions) if inference else "{}",
                "reason": "insufficient_panel_context" if template is None else ""})
            if inference:
                likelihood_rows.extend({"sample_id": sample_id, "locus_id": locus.locus_id, "repeat_count": r,
                    "log_likelihood": ll, "posterior": p} for r,ll,p in zip(inference.states, inference.log_likelihoods, inference.posterior))
            for item in items:
                row = {"sample_id": sample_id, "locus_id": locus.locus_id, "molecule_id": item.molecule_id,
                              "classes": ";".join(item.classes), "observed_repeat": item.observed_repeat,
                              "lower_bound": item.lower_bound, "fragment_offset": item.fragment_offset,
                              "orientation": ";".join(item.orientations), "alignment": json.dumps(item.alignment)}
                if audit_writer is None:
                    audit.append(row)
                else:
                    audit_writer.writerow(row)
        stage_seconds["sequence_reconstruction"] = time.perf_counter() - started
    progress.step(f"[{sample_id}] Writing canonical genotypes and evidence")
    started = time.perf_counter()
    paths = write_products(all_products, loci, output)
    paths.update(write_compatibility_products(all_products, loci, output))
    if technology == "illumina":
        paths["molecule_repeat_likelihoods"] = molecule_likelihood_path
    from .allele_inference import COMMON_LOCUS_CALL_FIELDS
    for key, filename, rows, fields in (
        ("common_locus_calls", "common_locus_calls.tsv", calls, COMMON_LOCUS_CALL_FIELDS),
        ("molecule_evidence", "molecule_candidate_evidence.tsv", audit,
         MOLECULE_AUDIT_FIELDS),
        ("short_read_repeat_evidence", "short_read_repeat_evidence.tsv", summaries,
         ["sample_id", "locus_id", "E", "S", "F", "FRR", "best_repeat", "second_repeat", "confidence", "interval_min", "interval_max", "effective_depth", "limit_reached", "identifiable", "components", "reason"]),
        ("repeat_likelihoods", "repeat_likelihoods.tsv", likelihood_rows,
         ["sample_id", "locus_id", "repeat_count", "log_likelihood", "posterior"]),
    ):
        paths[key] = output / filename
        if key != "molecule_evidence" or audit_writer is None:
            write_tsv(rows, paths[key], fields)
    dominant = {}
    for product in sorted(all_products, key=lambda p: -p.estimated_fraction):
        dominant.setdefault(product.locus_id, product.sequence)
    paths["taxonomic_query_sequences"] = output / "taxonomic_query_sequences.fasta"
    write_fasta(dominant.items(), paths["taxonomic_query_sequences"])
    paths["reconstruction_metadata"] = output / "reconstruction_metadata.json"
    stage_seconds["genotype_and_evidence_output"] = time.perf_counter() - started
    paths["reconstruction_metadata"].write_text(json.dumps({"technology": technology, "insert_size": vars(insert) if insert else None,
        "performance": {"stage_seconds": stage_seconds, "recruitment": recruitment_stats,
                        "loci": locus_stats}}, indent=2) + "\n")
    return calls, audit, {}, paths


def write_compatibility_products(products, loci, output):
    """Keep established variant, membership, mixture and read-call tables."""
    from .pipeline import ASV_FIELDS, ASV_MEMBERSHIP_FIELDS, PREDICTION_FIELDS
    from .mixture import MIXTURE_FIELDS
    by_id = {l.locus_id: l for l in loci}
    variants, memberships, predictions, mixtures = [], [], [], []
    totals = Counter()
    dominant = {}
    for product in sorted(products, key=lambda p: -p.estimated_fraction):
        totals[product.locus_id] += product.support_count
        dominant.setdefault(product.locus_id, product.variant_id)
    for p in products:
        genotype = genotype_product(p, by_id[p.locus_id])
        variants.append({"sample_id": p.sample_id, "locus_id": p.locus_id, "variant_id": p.variant_id,
                         "repeat_count": genotype.repeat_count, "support_reads": p.support_count,
                         "unique_sequences": 1, "frequency": p.estimated_fraction,
                         "representative_read_id": next(iter(p.evidence.get("molecule_ids", [])), ""),
                         "representative_sequence": p.sequence, "representative_length_bp": len(p.sequence)})
        for molecule in p.evidence.get("molecule_ids", []):
            membership = {"sample_id": p.sample_id, "read_id": molecule, "locus_id": p.locus_id,
                          "variant_id": p.variant_id, "repeat_count": genotype.repeat_count}
            memberships.append(membership)
            predictions.append({**membership, "predicted_repeat_count": genotype.repeat_count,
                                "probability": p.reconstruction_confidence,
                                "evidence_weight": 1, "raw_repeat_count_estimate": genotype.repeat_count_raw})
        mixtures.append({"sample_id": p.sample_id, "locus_id": p.locus_id, "variant_id": p.variant_id,
                         "repeat_count": genotype.repeat_count, "observed_reads": p.support_count,
                         "estimated_reads": p.effective_depth, "estimated_fraction": p.estimated_fraction,
                         "observed_fraction": p.support_count / max(totals[p.locus_id], 1),
                         "meaningful": p.evidence.get("meaningful", "yes"),
                         "abundance_class": "DOMINANT" if dominant[p.locus_id] == p.variant_id else
                                            "SECONDARY" if p.evidence.get("meaningful") == "yes" else "TRACE",
                         "evidence_class": "DOMINANT" if dominant[p.locus_id] == p.variant_id else
                                           "CONFIRMED_SECONDARY" if p.evidence.get("meaningful") == "yes" else "TRACE"})
    paths = {}
    for key, name, rows, fields in (
        ("mapped_variant_table", "mapped_variant_table.tsv", variants, ASV_FIELDS),
        ("mapped_read_memberships", "mapped_read_memberships.tsv", memberships, ASV_MEMBERSHIP_FIELDS),
        ("read_predictions", "read_level_allele_predictions.tsv", predictions, PREDICTION_FIELDS),
        ("mixture_abundance", "vntr_mixture_abundance.tsv", mixtures, MIXTURE_FIELDS),
    ):
        paths[key] = output / name
        write_tsv(rows, paths[key], fields)
    paths["mapped_variant_representatives"] = output / "mapped_variant_representatives.fasta.gz"
    write_fasta(((p.variant_id, p.sequence) for p in products), paths["mapped_variant_representatives"])
    paths.update({"asv_table": paths["mapped_variant_table"], "asv_memberships": paths["mapped_read_memberships"],
                  "asv_representatives": paths["mapped_variant_representatives"]})
    return paths
