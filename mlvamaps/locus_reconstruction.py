"""Input-specific recovery into canonical products; no known-allele selection."""
from __future__ import annotations

from collections import Counter, defaultdict
from contextlib import ExitStack
from dataclasses import replace
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import parasail

from .concurrency import resolve_threads
from .progress import ProgressReporter
from .io import read_fastq, read_fastq_pairs, write_fasta, write_tsv
from .locus_measurement import find_anchor
from .locus_products import LocusProduct, product_repeat_allele, write_products
from .mixture import estimate_variant_mixtures
from .models import Locus, ReadPair
from .short_read_evidence import (
    _read_profile, RepeatTemplate, estimate_insert_distribution, panel_template,
)
from .sequence import revcomp


def locus_templates(loci):
    """Build Illumina anchors exclusively from the supplied locus panel."""
    from .repeat_calibration import repeat_unit_length
    templates = {}
    for locus in loci:
        template = panel_template(locus)
        if template is None:
            if not locus.forward_primer or not locus.reverse_primer:
                raise ValueError(f"Illumina primer-only reconstruction requires both primers for {locus.locus_id}")
            motif = locus.repeat_motif.upper()
            template = RepeatTemplate(locus, locus.forward_primer, revcomp(locus.reverse_primer),
                motif if motif and not set(motif)-set('ACGT') else '', repeat_unit_length(locus),
                primer_only=True)
        templates[locus.locus_id] = template
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
            count = product_repeat_allele(LocusProduct(sample_id, locus.locus_id, "long_read", "molecule", product, calibrated_product_size_bp=calibrated_size), locus)[1] if complete else None
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
            anchors = item.alignment.get(str(mate), {})
            anchored = (position < repeat_start and anchors.get("left")) or (position >= repeat_end and anchors.get("right"))
            if reliable and anchored:
                votes[position].add((insertion, base))
            insertion = ""
            position += 1
    return {p: next(iter(values)) for p, values in votes.items() if len(values) == 1}


def _common_call(locus, products, recruited, sample_id, technology, minimum_molecules, minimum_probability, inference=None):
    ranked = sorted(products, key=lambda p: (-p.estimated_fraction, -p.support_count, p.variant_id))
    best = ranked[0] if ranked else None
    confidence = best.reconstruction_confidence if best else inference.confidence if inference else 0
    meaningful = [p for p in ranked if p.evidence.get("meaningful") == "yes"]
    status = "not_found" if not recruited else "detected_unresolved"
    if best and confidence >= minimum_probability:
        status = "low_coverage" if best.effective_depth < minimum_molecules else "mixed" if len(meaningful) > 1 else "called"
    if inference and inference.method == "MIXED":
        status = "mixed"
    repeat = product_repeat_allele(best, locus)[1] if best else ""
    if repeat is None:
        repeat = ""
        status = "detected_unresolved"
    count_distribution = defaultdict(float)
    for product in ranked:
        count = product_repeat_allele(product, locus)[1]
        if count is not None:
            count_distribution[count] += product.estimated_fraction
    counts = inference.counts if inference else {"FULL_SPAN": sum(p.support_count for p in ranked)}
    return {"sample": sample_id, "locus": locus.locus_id, "technology": technology,
            "repeat_count": repeat, "status": status, "best_probability": confidence,
            "second_best_probability": sorted(inference.posterior, reverse=True)[1] if inference and len(inference.posterior)>1 else 0,
            "molecule_support": recruited, "direct_product_support": counts.get("FULL_SPAN", 0),
            "full_span_support": counts.get("FULL_SPAN", 0), "junction_support": counts.get("LEFT_BOUNDARY", 0) + counts.get("RIGHT_BOUNDARY", 0),
            "dominant_fraction": best.estimated_fraction if best else "",
            "secondary_repeat": product_repeat_allele(ranked[1], locus)[1] if len(ranked)>1 else "",
            "secondary_fraction": ranked[1].estimated_fraction if len(ranked)>1 else "",
            "candidate_distribution": ";".join(f"{r}:{p:.6f}" for r,p in sorted(count_distribution.items(), key=lambda v:-v[1])),
            "confidence": confidence, "margin": confidence,
            "best_candidate_repeat": inference.best if inference else repeat,
            "dominant_repeat": repeat,
            "num_variants": len(ranked), "num_secondary": max(0, len(meaningful)-1),
            "inference_method": inference.method if inference else "spanning_molecules",
            "reason": inference.reason if inference else "",
            "n_spanning": counts.get("FLANK_PAIR", 0), "repeat_count_min": inference.interval[0] if inference and inference.interval else repeat,
            "repeat_count_max": inference.interval[1] if inference and inference.interval else repeat,
            "second_best_repeat_count": inference.second if inference else ""}


MOLECULE_AUDIT_FIELDS = [
    "sample_id", "locus_id", "molecule_id", "classes", "observed_repeat", "lower_bound",
    "fragment_offset", "orientation", "alignment", "complete", "repeat_count",
    "anchor_edits", "forward_start", "reverse_end", "variant_id",
]


def run_reconstructed_fastq_inference(*, outdir, stream_molecule_evidence=False, **kwargs):
    """Own diagnostic streams for the duration of inference."""
    if kwargs.get("technology") == "illumina":
        output = Path(outdir)
        output.mkdir(parents=True, exist_ok=True)
        with ExitStack() as stack:
            if stream_molecule_evidence:
                handle = stack.enter_context((output / "molecule_candidate_evidence.tsv").open("w", newline=""))
                writer = csv.DictWriter(handle, fieldnames=MOLECULE_AUDIT_FIELDS, delimiter="\t", extrasaction="ignore")
                writer.writeheader()
                kwargs.update(audit_writer=writer, audit_handle=handle)
            if kwargs.get("sr_recruitment_audit", "full") == "compact":
                handle = stack.enter_context((output / "ambiguous_molecules.tsv").open("w", newline=""))
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["sample_id", "molecule_id", "witness_locus_1", "witness_locus_2", "candidate_search_complete"])
                kwargs["ambiguity_writer"] = writer
            return _run_reconstructed_fastq_inference(outdir=outdir, **kwargs)
    return _run_reconstructed_fastq_inference(outdir=outdir, **kwargs)


def _run_reconstructed_fastq_inference(*, reads1, reads2, loci, database_path, outdir, sample_id,
        technology, minimum_molecules, minimum_probability, maximum_candidate_repeat_count=100,
        insert_mean=None, insert_sd=None, orphan_path=None, min_fraction=.01, min_secondary_reads=2,
        max_anchor_edits=3, threads=1, minimum_spanning_pairs=2, round_tolerance=.25,
        show_progress=False, audit_writer=None, audit_handle=None, pairs=None,
        recruitment_threads=None, sr_recruitment_audit="full", ambiguity_writer=None, repeat_threshold=.7, **_unused):
    progress = ProgressReporter(enabled=show_progress, stream=sys.stdout)
    stage_seconds, recruitment_stats, locus_stats = {}, {}, {}
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    if pairs is None:
        pairs = read_fastq_pairs(reads1, reads2)
    if orphan_path:
        from itertools import chain
        pairs = chain(pairs, (ReadPair(r.read_id, r) for r in read_fastq(orphan_path)))
    all_products, calls, summaries, likelihood_rows, audit = [], [], [], [], []
    reconstruction_paths = {}
    insert = None
    if technology != "illumina":
        all_products, audit, recruited = recover_long_reads(pairs, loci, sample_id, min_fraction, min_secondary_reads, max_anchor_edits)
        all_products = [replace(p, round_tolerance=round_tolerance) for p in all_products]
        for locus in loci:
            products = [p for p in all_products if p.locus_id == locus.locus_id]
            calls.append(_common_call(locus, products, recruited[locus.locus_id], sample_id, technology, minimum_molecules, minimum_probability))
    else:
        progress.step(f"[{sample_id}] Loading repeat templates for {len(loci)} loci")
        started = time.perf_counter()
        templates = locus_templates(loci)
        stage_seconds["template_loading"] = time.perf_counter() - started
        recruitment_stats["template_contexts"] = len(templates)
        recruitment_stats["template_source"] = "locus_panel"
        recruitment_stats["primer_only_loci"] = [name for name, t in templates.items() if t.primer_only]
        recruitment_stats["template_motif_lengths"] = {
            name: len(template.motif) for name, template in templates.items() if template is not None
        }
        progress.step(f"[{sample_id}] Repeat templates loaded in {stage_seconds['template_loading']:.1f}s; {len(templates):,} panel templates (no reference database sequences)")
        evidence = defaultdict(list)
        insert_lengths = []
        from .short_read_recruitment import recruit_short_reads
        recruitment_threads = threads if recruitment_threads is None else recruitment_threads
        progress.step(f"[{sample_id}] Recruiting short-read molecules with up to {resolve_threads(recruitment_threads)} worker(s); {sr_recruitment_audit} ambiguity audit")
        started = time.perf_counter()
        recruitment_stats.update(pairs_examined=0, locus_tests=0, uniquely_recruited_pairs=0,
                                 ambiguous_pairs=0, unmatched_pairs=0, skipped_locus_tests=0)
        for chunk in recruit_short_reads(pairs, templates, sample_id, recruitment_threads, statistics=recruitment_stats,
                                        audit_fields=MOLECULE_AUDIT_FIELDS if audit_handle is not None else None,
                                        audit_mode=sr_recruitment_audit, repeat_threshold=repeat_threshold):
            recruitment_stats["pairs_examined"] += chunk.examined
            recruitment_stats["uniquely_recruited_pairs"] += len(chunk.evidence)
            for key in ("locus_tests", "ambiguous_pairs", "unmatched_pairs", "skipped_locus_tests"):
                recruitment_stats[key] += getattr(chunk, key)
            if ambiguity_writer is not None:
                ambiguity_writer.writerows((sample_id, *row) for row in chunk.ambiguity_witnesses)
            for item in chunk.evidence:
                evidence[item.locus_id].append(item)
            if audit_handle is not None:
                audit_handle.write(chunk.audit_tsv)
            elif audit_writer is None:
                audit.extend(chunk.ambiguous)
            else:
                audit_writer.writerows(chunk.ambiguous)
            insert_lengths.extend(chunk.insert_lengths)
            rate = recruitment_stats["pairs_examined"] / max(time.perf_counter() - started, 1e-9)
            progress.count(f"[{sample_id}] Read pairs scanned", recruitment_stats["pairs_examined"],
                           detail=f"{recruitment_stats['uniquely_recruited_pairs']:,} unique, "
                                  f"{recruitment_stats['ambiguous_pairs']:,} ambiguous; {rate:,.0f} pairs/s")
        stage_seconds["recruitment"] = time.perf_counter() - started
        recruitment_stats["pairs_per_second"] = recruitment_stats["pairs_examined"] / max(stage_seconds["recruitment"], 1e-9)
        progress.step(f"[{sample_id}] Recruitment finished in {stage_seconds['recruitment']:.1f}s; fitting {len(loci)} loci")
        insert = estimate_insert_distribution(insert_lengths, insert_mean, insert_sd)
        from .targeted_reconstruction import recover_locus, EVIDENCE_CLASSES
        def fit_locus(locus):
            started = time.perf_counter()
            items = evidence[locus.locus_id]
            result = recover_locus(items, templates[locus.locus_id], sample_id, insert,
                maximum_candidate_repeat_count, minimum_probability, minimum_spanning_pairs,
                min_fraction, min_secondary_reads)
            locus_stats[locus.locus_id] = {"recovery_seconds": time.perf_counter()-started,
                "evidence_molecules": len(items), "candidate_states": len(result.states), "call_method": result.method}
            return result
        from concurrent.futures import ThreadPoolExecutor
        started = time.perf_counter()
        if resolve_threads(threads) > 1 and len(loci) > 1:
            with ThreadPoolExecutor(max_workers=min(len(loci), resolve_threads(threads))) as executor:
                fitted = list(executor.map(fit_locus, loci))
        else:
            fitted = [fit_locus(locus) for locus in loci]
        progress.step(f"[{sample_id}] Measuring reconstructed loci with assembly Sassy PCR")
        pcr_started = time.perf_counter()
        from .local_assembly import measure_reconstructed_loci
        reconstruction_paths = measure_reconstructed_loci(fitted, loci, output/'locus_reconstruction',
            sample_id, max_anchor_edits, round_tolerance)
        stage_seconds['reconstruction_sassy_pcr'] = time.perf_counter() - pcr_started
        for locus, result in zip(loci, fitted):
            locus_stats[locus.locus_id]['call_method'] = result.method
        stage_seconds["locus_recovery"] = time.perf_counter() - started
        progress.step(f"[{sample_id}] Locus recovery finished in {stage_seconds['locus_recovery']:.1f}s")
        for locus, result in zip(loci, fitted):
            items = evidence[locus.locus_id]
            usable = [e for e in items if "discordant" not in e.classes]
            products = [replace(p, round_tolerance=round_tolerance) for p in result.products]
            all_products.extend(products)
            calls.append(_common_call(locus, products, len(items), sample_id, technology,
                                      minimum_molecules, minimum_probability, result))
            read_lengths = [len(seq) for e in items for seq in e.sequences]
            product_length = len(products[0].sequence) if products else None
            template = templates[locus.locus_id]
            summaries.append({"sample_id": sample_id, "locus_id": locus.locus_id,
                "call_method": result.method, "confidence": result.confidence,
                **{("n_flank_pairs" if name == "FLANK_PAIR" else "n_"+name.lower()): result.counts.get(name, 0) for name in EVIDENCE_CLASSES},
                "best_repeat": calls[-1]['repeat_count'], "second_repeat": result.second,
                "interval_min": result.interval[0] if result.interval else calls[-1]['repeat_count'],
                "interval_max": result.interval[1] if result.interval else calls[-1]['repeat_count'],
                "effective_depth": len(usable), "limit_reached": result.limit_reached,
                "read_length": float(np.median(read_lengths)) if read_lengths else '',
                "motif_length": template.unit if template else '', "amplicon_length": product_length or '',
                "vntr_length": product_length-len(template.left)-len(template.right) if product_length and template and not template.primer_only else '',
                "amplicon_insert_ratio": product_length/insert.mean if product_length and insert else '',
                "reason": result.reason})
            likelihood_rows.extend({"sample_id": sample_id, "locus_id": locus.locus_id,
                "repeat_count": r, "log_likelihood": ll, "posterior": probability}
                for r,ll,probability in zip(result.states, result.log_likelihoods, result.posterior))
            for item in items:
                row = {"sample_id": sample_id, "locus_id": locus.locus_id, "molecule_id": item.molecule_id,
                       "classes": ";".join(item.classes), "observed_repeat": item.observed_repeat,
                       "lower_bound": item.lower_bound, "fragment_offset": item.fragment_offset,
                       "orientation": ";".join(item.orientations), "alignment": json.dumps(item.alignment)}
                if audit_writer is None:
                    audit.append(row)
                else:
                    audit_writer.writerow(row)
    progress.step(f"[{sample_id}] Writing canonical genotypes and evidence")
    started = time.perf_counter()
    output_stats = {}
    paths = write_products(all_products, loci, output, progress=progress, statistics=output_stats)
    paths.update(reconstruction_paths)
    compatibility_started = time.perf_counter()
    progress.step(f"[{sample_id}] Writing molecule memberships and allele predictions")
    paths.update(write_compatibility_products(all_products, loci, output, progress=progress))
    output_stats['compatibility_tables_seconds'] = time.perf_counter() - compatibility_started
    progress.step(f"[{sample_id}] Writing locus calls and evidence summaries")
    if technology == "illumina":
        if ambiguity_writer is not None:
            paths["ambiguous_molecules"] = output / "ambiguous_molecules.tsv"
    from .allele_inference import COMMON_LOCUS_CALL_FIELDS
    for key, filename, rows, fields in (
        ("common_locus_calls", "common_locus_calls.tsv", calls, COMMON_LOCUS_CALL_FIELDS),
        ("molecule_evidence", "molecule_candidate_evidence.tsv", audit,
         MOLECULE_AUDIT_FIELDS),
        ("short_read_repeat_evidence", "short_read_repeat_evidence.tsv", summaries,
         list(summaries[0]) if summaries else ["sample_id", "locus_id", "call_method", "confidence"]),
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
                        "loci": locus_stats, "output": output_stats}}, indent=2) + "\n")
    progress.step(f"[{sample_id}] Canonical genotype and evidence output finished in {stage_seconds['genotype_and_evidence_output']:.1f}s")
    return calls, audit, {}, paths


def write_compatibility_products(products, loci, output, *, progress=None):
    """Keep established variant, membership, mixture and read-call tables."""
    from .pipeline import ASV_FIELDS, ASV_MEMBERSHIP_FIELDS, PREDICTION_FIELDS
    from .mixture import MIXTURE_FIELDS
    by_id = {l.locus_id: l for l in loci}
    variants, mixtures, repeats = [], [], []
    totals = Counter()
    dominant = {}
    for product in sorted(products, key=lambda p: -p.estimated_fraction):
        totals[product.locus_id] += product.support_count
        dominant.setdefault(product.locus_id, product.variant_id)
    for p in products:
        raw_repeat, repeat = product_repeat_allele(p, by_id[p.locus_id])
        repeats.append((raw_repeat, repeat))
        variants.append({"sample_id": p.sample_id, "locus_id": p.locus_id, "variant_id": p.variant_id,
                         "repeat_count": repeat, "support_reads": p.support_count,
                         "unique_sequences": 1, "frequency": p.estimated_fraction,
                         "representative_read_id": next(iter(p.evidence.get("molecule_ids", [])), ""),
                         "representative_sequence": p.sequence, "representative_length_bp": len(p.sequence)})
        mixtures.append({"sample_id": p.sample_id, "locus_id": p.locus_id, "variant_id": p.variant_id,
                         "repeat_count": repeat, "observed_reads": p.support_count,
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
        ("mixture_abundance", "vntr_mixture_abundance.tsv", mixtures, MIXTURE_FIELDS),
    ):
        paths[key] = output / name
        write_tsv(rows, paths[key], fields)
    # Native CSV writers stream repeated molecule rows. Keep one row template
    # per product, rather than two dictionaries per molecule for the sample.
    from itertools import islice
    for key, name, fields in (
        ("mapped_read_memberships", "mapped_read_memberships.tsv", ASV_MEMBERSHIP_FIELDS),
        ("read_predictions", "read_level_allele_predictions.tsv", PREDICTION_FIELDS),
    ):
        paths[key] = output / name
        with paths[key].open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(fields)
            read_column = fields.index('read_id')
            written = 0
            for p, (raw_repeat, repeat) in zip(products, repeats):
                values = {"sample_id": p.sample_id, "locus_id": p.locus_id, "variant_id": p.variant_id,
                          "repeat_count": repeat, "predicted_repeat_count": repeat,
                          "probability": p.reconstruction_confidence, "evidence_weight": 1,
                          "raw_repeat_count_estimate": raw_repeat}
                row = [values.get(field, '') for field in fields]
                molecules = iter(p.evidence.get('molecule_ids', []))
                while batch := list(islice(molecules, 8192)):
                    writer.writerows(row[:read_column] + [molecule] + row[read_column + 1:]
                                     for molecule in batch)
                    written += len(batch)
                    if progress is not None:
                        progress.count(f"[{p.sample_id}] {name} rows written", written)
    paths["mapped_variant_representatives"] = output / "mapped_variant_representatives.fasta.gz"
    write_fasta(((p.variant_id, p.sequence) for p in products), paths["mapped_variant_representatives"])
    paths.update({"asv_table": paths["mapped_variant_table"], "asv_memberships": paths["mapped_read_memberships"],
                  "asv_representatives": paths["mapped_variant_representatives"]})
    return paths
