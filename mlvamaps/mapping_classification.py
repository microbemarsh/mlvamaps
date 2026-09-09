"""Reference likelihoods from observed VNTR alignments; Emu-inspired mixture EM.

Scores are relative, approximate likelihoods, not calibrated species probabilities.
Repeat-length indels are scored once, separately from sequence errors. The
unclassified component allows evidence poorly explained by the catalog to remain
unassigned. No reference-guided query consensus is used for read classification.
"""
from __future__ import annotations

import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from .alignment_evidence import CandidateAlignment
from .candidate_contexts import CandidateContext, repeat_interval
from .io import write_fasta, write_tsv
from .minimap_mapping import map_reads_to_candidates_bam
from .models import Locus
from .reference_database import read_reference_metadata, read_sequence_database
from .repeat_calibration import assembly_equivalent_product_allele, repeat_unit_length

MISSING_LOCUS_FIELDS = [
    "missing_locus_penalty",
    "penalized_missing_loci",
    "missing_locus_gate_passed",
    "well_covered_locus_fraction",
    "missing_locus_min_depth",
    "missing_locus_min_fraction",
    "missing_locus_penalty_per_locus",
]


def _coverage_gated_missing_loci(
    requested_loci: set[str],
    query_sequences: dict[str, str],
    locus_quality: dict[str, dict[str, object]] | None,
    input_mode: str,
    min_depth: float,
    min_fraction: float,
    penalty: float,
) -> tuple[set[str], dict]:
    """Use molecule support as a coverage proxy, never as confirmed absence."""
    if not math.isfinite(min_depth) or min_depth <= 0:
        raise ValueError("Missing-locus minimum depth must be finite and positive")
    if not math.isfinite(min_fraction) or not 0 < min_fraction <= 1:
        raise ValueError("Missing-locus minimum fraction must be in (0, 1]")
    if not math.isfinite(penalty) or penalty < 0:
        raise ValueError("Missing-locus penalty must be finite and non-negative")
    covered = set()
    missing = set()
    if input_mode in {"fastq", "illumina"}:
        for locus_id in requested_loci:
            quality = (locus_quality or {}).get(locus_id, {})
            status = str(quality.get("detection_status", quality.get("status", ""))).lower()
            try:
                depth = float(quality.get("depth", ""))
            except (TypeError, ValueError):
                continue
            if not math.isfinite(depth) or depth < 0:
                continue
            if status == "not_found":
                if depth == 0 and not query_sequences.get(locus_id):
                    missing.add(locus_id)
            elif depth >= min_depth and status in {
                "called", "pass", "low_coverage", "low_depth", "ambiguous",
                "mixed", "multiple_variants", "detected_unresolved",
                "present_count_unknown", "present",
            }:
                covered.add(locus_id)
    fraction = len(covered) / len(requested_loci) if requested_loci else 0.0
    passed = fraction >= min_fraction
    return (missing if passed and penalty > 0 else set()), {
        "missing_locus_gate_passed": "yes" if passed else "no",
        "well_covered_locus_fraction": f"{fraction:.6f}",
        "missing_locus_min_depth": min_depth,
        "missing_locus_min_fraction": min_fraction,
        "missing_locus_penalty_per_locus": penalty,
    }


_CS = re.compile(r":(\d+)|=([A-Za-z]+)|\*([A-Za-z]{2})|\+([A-Za-z]+)|-([A-Za-z]+)")
_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
UNKNOWN = "__unclassified__"
MATCH_FIELDS = [
    "sample_id", "match_type", "reference_id", "equivalent_references", "taxon_id", "taxon_name",
    "rank", "distance", "log_likelihood", "model_weight", "locus_balanced_fraction",
    "compared_loci", "match_status", *MISSING_LOCUS_FIELDS,
]


def observed_reference_contexts(database_path, loci):
    root = Path(database_path)
    if (root / "database").is_dir():
        root = root / "database"
    records = read_sequence_database(root, {l.locus_id for l in loci})
    by_id = {l.locus_id: l for l in loci}
    collapsed = {}
    for locus_id, references in sorted(records.items()):
        for reference_id, sequence in references:
            if reference_id == UNKNOWN:
                raise ValueError(f"Reference ID {UNKNOWN!r} is reserved")
            collapsed.setdefault((locus_id, sequence.upper()), set()).add(reference_id)
    contexts, members = [], {}
    for index, ((locus_id, sequence), references) in enumerate(sorted(collapsed.items())):
        locus = by_id[locus_id]
        start, end = repeat_interval(sequence, locus)
        count = assembly_equivalent_product_allele(locus, len(sequence))[1]
        unit = repeat_unit_length(locus)
        if count is None and unit and locus.left_flank_sequence and locus.right_flank_sequence:
            count = (end - start) / unit
        candidate_id = f"observed{index:07d}"
        contexts.append(CandidateContext(
            candidate_id, locus_id, sequence, count, unit, start, end,
            locus.left_flank_sequence, locus.right_flank_sequence,
            observed_or_synthetic="observed", source_repeat_count=count,
        ))
        members[candidate_id] = sorted(references)
    if not contexts:
        raise ValueError("No observed reference VNTR sequences match the panel")
    return root, contexts, members


def alignment_statistics(alignment, context):
    """Count sequence errors and separately measure repeat-length evidence.

    cs retains mismatches even for secondary SAM records without SEQ. Assembly
    alignments supply extended CIGAR (= / X). Generic M without cs is rejected.
    """
    operations = []
    if alignment.cs:
        tokens = list(_CS.finditer(alignment.cs))
        if "".join(token.group(0) for token in tokens) != alignment.cs:
            raise ValueError("Unsupported or malformed cs tag in classification alignment")
        for token in tokens:
            match, identical, mismatch, insertion, deletion = token.groups()
            operations.append((
                "=" if match or identical else "X" if mismatch else "I" if insertion else "D",
                int(match) if match else len(identical or insertion or deletion or "x"),
            ))
    else:
        operations = [(op, int(n)) for n, op in _CIGAR.findall(alignment.cigar)]
        if any(op == "M" for op, _ in operations):
            raise ValueError("Classification requires cs tags or extended CIGAR")
    position = alignment.reference_start
    matched = mismatches = insertions = deletions = repeat_delta = 0
    repeat_insertions = repeat_deletions = 0
    outside_bases = 0
    for operation, length in operations:
        if operation in "SH":
            continue
        if operation not in "=XID":
            raise ValueError(f"Unsupported classification CIGAR operation: {operation}")
        overlap = (
            length if operation == "I" and (
                context.repeat_start <= position < context.repeat_end
                or position == context.repeat_start == context.repeat_end
            )
            else 0 if operation == "I"
            else max(0, min(position + length, context.repeat_end) - max(position, context.repeat_start))
        )
        outside = length - overlap
        if operation == "=":
            matched += length
            outside_bases += outside
        elif operation == "X":
            mismatches += length
            outside_bases += outside
        elif operation == "I":
            insertions += outside
            repeat_delta += overlap
            repeat_insertions += overlap
        else:
            deletions += outside
            repeat_delta -= overlap
            repeat_deletions += overlap
        if operation != "I":
            position += length
    spanning = (
        context.repeat_count is not None and context.repeat_unit_length > 0
        and context.repeat_end >= context.repeat_start
        and alignment.reference_start <= context.repeat_start - 3
        and alignment.reference_end >= context.repeat_end + 3
    )
    # Opposing indels with zero net length are sequence differences, not a
    # repeat-count change. Only the net length event belongs to the repeat term.
    balanced_indels = min(repeat_insertions, repeat_deletions)
    return {
        "matches": matched, "mismatches": mismatches, "insertions": insertions + balanced_indels,
        "deletions": deletions + balanced_indels, "outside_bases": outside_bases,
        "repeat_delta": repeat_delta / context.repeat_unit_length if spanning else None,
    }


def molecule_log_likelihoods(alignments, contexts, members, repeat_scale=1.0):
    if not math.isfinite(repeat_scale) or repeat_scale <= 0:
        raise ValueError("Classification repeat scale must be finite and positive")
    context_by_id = {c.candidate_id: c for c in contexts}
    best = {}
    for alignment in alignments:
        context = context_by_id.get(alignment.candidate_id)
        if context is None or not math.isfinite(alignment.alignment_score):
            continue
        key = (alignment.molecule_id, alignment.mate, alignment.candidate_id)
        if key not in best or alignment.alignment_score > best[key].alignment_score:
            best[key] = alignment
    stats = {key: alignment_statistics(row, context_by_id[row.candidate_id]) for key, row in best.items()}
    best = {key: row for key, row in best.items() if stats[key]["outside_bases"] >= 6}
    # Estimate errors from the best observed alignment per read, with pseudocounts.
    primary = {}
    for key, row in best.items():
        read = key[:2]
        if read not in primary or row.alignment_score > best[primary[read]].alignment_score:
            primary[read] = key
    totals = Counter({"matches": 100, "mismatches": 1, "insertions": 1, "deletions": 1})
    for key in primary.values():
        for field in ("matches", "mismatches", "insertions", "deletions"):
            totals[field] += stats[key][field]
    denominator = sum(totals.values())
    errors = {field: min(0.25, max(1e-4, totals[field] / denominator)) for field in ("mismatches", "insertions", "deletions")}
    logs = {field: math.log(value / (3 if field == "mismatches" else 1)) for field, value in errors.items()}
    maximum_span = defaultdict(int)
    locus_scores = defaultdict(dict)
    for key, row in best.items():
        maximum_span[(row.molecule_id, row.mate, row.locus_id)] = max(
            maximum_span[(row.molecule_id, row.mate, row.locus_id)], stats[key]["outside_bases"]
        )
        locus_scores[row.molecule_id][row.locus_id] = max(
            locus_scores[row.molecule_id].get(row.locus_id, -math.inf), row.alignment_score
        )
    # A molecule mapping equally well to different loci is not independent
    # evidence at each of them. Exclude locus ties rather than duplicating it.
    assigned_locus = {}
    for molecule, scores in locus_scores.items():
        ranked = sorted(scores, key=lambda locus: (-scores[locus], locus))
        if len(ranked) == 1 or scores[ranked[0]] > scores[ranked[1]]:
            assigned_locus[molecule] = ranked[0]
    by_molecule = defaultdict(lambda: defaultdict(dict))
    for key, row in best.items():
        if assigned_locus.get(row.molecule_id) != row.locus_id:
            continue
        values = stats[key]
        sequence_score = sum(values[field] * logs[field] for field in logs)
        uncovered = maximum_span[(row.molecule_id, row.mate, row.locus_id)] - values["outside_bases"]
        sequence_score += uncovered * logs["insertions"]
        repeat_score = -abs(values["repeat_delta"]) / repeat_scale if values["repeat_delta"] is not None else 0.0
        for reference in members[row.candidate_id]:
            by_molecule[(row.locus_id, row.molecule_id)][reference][row.mate] = (sequence_score, repeat_score)
    results = {}
    for key, references in sorted(by_molecule.items()):
        mates = {mate for scores in references.values() for mate in scores}
        per_mate_floor = {
            mate: min(scores[mate][0] for scores in references.values() if mate in scores) - 8.0
            for mate in mates
        }
        results[key] = {}
        for reference, scores in references.items():
            # Both mates contribute sequence evidence; their shared repeat length
            # is one molecule-level observation, not two independent counts.
            sequence_score = sum(scores[mate][0] if mate in scores else per_mate_floor[mate] for mate in mates)
            repeat_score = min((value[1] for value in scores.values()), default=0.0)
            results[key][reference] = sequence_score + repeat_score
        span = sum(maximum_span[(key[1], mate, key[0])] for mate in mates)
        results[key][UNKNOWN] = -max(4.0, 0.1 * span)
    return results, errors


def fit_reference_mixture(log_likelihoods, weights, max_iterations=500, tolerance=1e-8):
    """Weighted EM; returns fractions and a monotone objective history."""
    values = np.asarray(log_likelihoods, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if max_iterations < 1 or not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("EM iteration limit and tolerance must be positive")
    if values.ndim != 2 or not all(values.shape) or weights.shape != (values.shape[0],):
        raise ValueError("EM requires a nonempty observations-by-references matrix and observation weights")
    if not np.isfinite(values).all() or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError("EM likelihoods must be finite and weights positive")
    shifted = values.max(axis=1)
    likelihood = np.exp(np.maximum(values - shifted[:, None], -700))
    fractions = np.full(values.shape[1], 1 / values.shape[1])
    history = [float(weights @ (np.log(likelihood @ fractions) + shifted))]
    for _ in range(max_iterations):
        responsibilities = likelihood * fractions
        responsibilities /= responsibilities.sum(axis=1, keepdims=True)
        updated = weights @ responsibilities / weights.sum()
        objective = float(weights @ (np.log(likelihood @ updated) + shifted))
        if objective < history[-1] - 1e-8:
            raise RuntimeError("Reference EM objective decreased")
        history.append(objective)
        change = np.max(np.abs(updated - fractions))
        fractions = updated
        if change < tolerance:
            break
    return fractions, history


def classify_molecules(likelihoods, reference_ids, *, sample_mode="isolate", penalties=None):
    if sample_mode not in {"isolate", "metagenome"}:
        raise ValueError("Classification sample mode must be isolate or metagenome")
    if any(not math.isfinite(value) or value < 0 for value in (penalties or {}).values()):
        raise ValueError("Reference penalties must be finite and non-negative")
    references = sorted(set(reference_ids) - {UNKNOWN}) + [UNKNOWN]
    observations = sorted(likelihoods)
    if not observations:
        return [], {"converged": True, "iterations": 0, "objective_history": []}
    depths = Counter(locus for locus, _ in observations)
    weights = np.array([min(depths[locus], 20) / depths[locus] for locus, _ in observations])
    # ponytail: dense molecule x reference matrix; use sparse blocks if catalogs
    # grow beyond the memory budget of the existing in-memory alignment reader.
    matrix = np.array([
        [likelihoods[key].get(ref, min(likelihoods[key].values()) - 8.0) for ref in references]
        for key in observations
    ], dtype=float)
    if not np.isfinite(matrix).all():
        raise ValueError("Reference log likelihoods must be finite")
    for j, reference in enumerate(references):
        matrix[:, j] -= (penalties or {}).get(reference, 0.0) / weights.sum()
    # Identical evidence columns share one prior and one EM component. Duplicate
    # references must not add taxon support just because the catalog is larger.
    equivalent = {}
    for j, ref in enumerate(references):
        signature = (ref == UNKNOWN, tuple(np.round(matrix[:, j], 10)))
        equivalent.setdefault(signature, []).append(j)
    groups = list(equivalent.values())
    grouped_matrix = matrix[:, [indexes[0] for indexes in groups]]
    joint = weights @ grouped_matrix
    model_weights = np.exp(np.maximum(joint - joint.max(), -700))
    model_weights /= model_weights.sum()
    if sample_mode == "metagenome":
        fractions, history = fit_reference_mixture(grouped_matrix, weights)
    else:
        fractions, history = model_weights, []
    rows = [{
        "references": [references[j] for j in indexes],
        "log_likelihood": float(joint[i]), "model_weight": float(model_weights[i]),
        "locus_balanced_fraction": float(fractions[i]),
    } for i, indexes in enumerate(groups)]
    rows.sort(key=lambda row: (-row["locus_balanced_fraction"], -row["log_likelihood"], row["references"]))
    return rows, {
        "iterations": max(0, len(history) - 1), "objective_history": history,
        "converged": len(history) < 501, "molecules": len(observations),
        "loci": len(depths), "effective_molecules": float(weights.sum()),
        "sample_mode": sample_mode,
    }


def _assembly_alignments(queries, contexts):
    import parasail
    matrix = parasail.matrix_create("ACGTN", 2, -4)
    rows = []
    for context in contexts:
        query = queries.get(context.locus_id)
        if not query:
            continue
        # Amplicons are primer-oriented by the assembly caller.
        result = parasail.nw_trace_striped_32(query, context.sequence, 5, 1, matrix)
        rows.append(CandidateAlignment(
            context.locus_id, context.locus_id, None, context.locus_id,
            context.candidate_id, context.repeat_count, "", float(result.score),
            0, 1.0, 1.0, 1.0, 0, len(query), 0, len(context.sequence),
            result.cigar.decode.decode(), "", True, False, False, query,
        ))
    return rows


def run_mapping_classification(
    *, database_path, loci, outdir, sample_id, reads1=None, reads2=None,
    query_sequences=None, technology="hifi", threads=1, minimap2_bin="minimap2",
    sample_mode="isolate", locus_quality=None, reference_metadata_path=None,
    taxon_identification=None, minimum_loci=2, repeat_scale=1.0,
    missing_locus_min_depth=3.0, missing_locus_min_fraction=0.8,
    missing_locus_penalty=8.0, keep_alignments=False,
    query_repeat_counts=None,
):
    output = Path(outdir) / "classification"
    output.mkdir(parents=True, exist_ok=True)
    root, contexts, members = observed_reference_contexts(database_path, loci)
    if minimum_loci < 1:
        raise ValueError("Minimum classification loci must be positive")
    metadata_path = reference_metadata_path or (root / "reference_metadata.tsv" if root.is_dir() else None)
    metadata = read_reference_metadata(metadata_path) if reference_metadata_path or (metadata_path and Path(metadata_path).is_file()) else {}
    fasta = output / "observed_reference_vntrs.fasta"
    write_fasta(((c.candidate_id, c.sequence) for c in contexts), fasta)
    write_tsv(({"candidate_id": c.candidate_id, "locus_id": c.locus_id,
        "reference_ids": ";".join(members[c.candidate_id]), "repeat_count": c.repeat_count,
        "repeat_start": c.repeat_start, "repeat_end": c.repeat_end} for c in contexts),
        output / "observed_reference_metadata.tsv",
        ["candidate_id", "locus_id", "reference_ids", "repeat_count", "repeat_start", "repeat_end"])
    if reads1 is not None:
        bam = output / "reference_alignments.bam"
        alignments = map_reads_to_candidates_bam(
            fasta, reads1, reads2, contexts, bam, threads, technology,
            executable=minimap2_bin, max_secondary=max(100, len(contexts)),
            include_unknown_repeats=True, retain_all_competitors=True,
        )
        if not keep_alignments:
            bam.unlink(missing_ok=True)
    else:
        alignments = _assembly_alignments(query_sequences or {}, contexts)
        # Keep the observed products and reference-relative evidence reviewable
        # independently of the final reference ranking.
        write_fasta(sorted((query_sequences or {}).items()), output / "query_amplicons.fasta")
        context_by_id = {context.candidate_id: context for context in contexts}
        assembly_evidence = []
        for alignment in alignments:
            context = context_by_id[alignment.candidate_id]
            statistics = alignment_statistics(alignment, context)
            assembly_evidence.append({
                "locus_id": alignment.locus_id,
                "reference_ids": ";".join(members[alignment.candidate_id]),
                "candidate_id": alignment.candidate_id,
                "query_repeat_count": (query_repeat_counts or {}).get(alignment.locus_id, ""),
                "reference_repeat_count": context.repeat_count,
                "cigar": alignment.cigar,
                **statistics,
            })
        write_tsv(assembly_evidence, output / "assembly_reference_evidence.tsv", [
            "locus_id", "reference_ids", "candidate_id", "query_repeat_count",
            "reference_repeat_count", "cigar", "matches", "mismatches",
            "insertions", "deletions", "outside_bases", "repeat_delta",
        ])
    likelihoods, errors = molecule_log_likelihoods(alignments, contexts, members, repeat_scale)
    reference_loci = defaultdict(set)
    for context in contexts:
        for reference in members[context.candidate_id]:
            reference_loci[reference].add(context.locus_id)
    # Mapping itself may rescue a locus unresolved by the allele caller.
    detected = {locus: "detected" for locus, _ in likelihoods}
    missing, gate = _coverage_gated_missing_loci(
        {l.locus_id for l in loci}, detected, locus_quality,
        "fastq" if reads1 is not None else "assembly", missing_locus_min_depth,
        missing_locus_min_fraction, missing_locus_penalty,
    )
    penalties = {ref: len(missing & expected) * missing_locus_penalty for ref, expected in reference_loci.items()}
    groups, diagnostics = classify_molecules(likelihoods, reference_loci, sample_mode=sample_mode, penalties=penalties)
    likelihood_path = output / "molecule_reference_likelihoods.tsv.gz"
    write_tsv(({
        "locus_id": locus, "molecule_id": molecule, "reference_id": reference,
        "log_likelihood": score,
    } for (locus, molecule), scores in likelihoods.items() for reference, score in scores.items()),
        likelihood_path, ["locus_id", "molecule_id", "reference_id", "log_likelihood"])
    matches = []
    unknown_fraction = 0.0 if groups else 1.0
    for group in groups:
        refs = group["references"]
        if refs == [UNKNOWN]:
            unknown_fraction += group["locus_balanced_fraction"]
            continue
        taxa = {str(metadata.get(ref, {}).get("taxon_id") or metadata.get(ref, {}).get("taxid") or "") for ref in refs}
        reference = refs[0]
        matches.append({
            "sample_id": sample_id, "match_type": "mapping_reference",
            "reference_id": reference, "equivalent_references": ";".join(refs),
            "taxon_id": ";".join(sorted(taxa - {""})),
            "taxon_name": "; ".join(sorted({
                str(metadata.get(ref, {}).get("taxon_name") or metadata.get(ref, {}).get("species")
                    or metadata.get(ref, {}).get("organism_name") or "") for ref in refs
            } - {""})), "rank": len(matches) + 1,
            "distance": -group["log_likelihood"], **{key: value for key, value in group.items() if key != "references"},
            "compared_loci": len(detected),
            "match_status": "EQUIVALENT_REFERENCES" if len(refs) > 1 else "MAPPING_LIKELIHOOD",
            **gate, "missing_locus_penalty": penalties[reference],
            "penalized_missing_loci": ",".join(sorted(missing & reference_loci[reference])),
        })
    match_path = output / "mapping_reference_matches.tsv"
    write_tsv(matches, match_path, MATCH_FIELDS)
    bands = []
    if matches:
        best_ref = matches[0]["reference_id"]
        for context in contexts:
            if best_ref in members[context.candidate_id]:
                bands.append({"reference_id": best_ref, "locus_id": context.locus_id,
                    "product_size_bp": len(context.sequence), "repeat_count": context.repeat_count if context.repeat_count is not None else ""})
    bands_path = output / "closest_reference_bands.tsv"
    write_tsv(bands, bands_path, ["reference_id", "locus_id", "product_size_bp", "repeat_count"])
    paths = {"classification": output, "mapping_reference_matches": match_path,
        "closest_reference_bands": bands_path,
        "mapping_likelihoods": likelihood_path}
    if reads1 is None:
        paths["assembly_reference_evidence"] = output / "assembly_reference_evidence.tsv"
        paths["classification_query_amplicons"] = output / "query_amplicons.fasta"
    paths.update(write_profile_tree(
        contexts, members, query_repeat_counts or {}, matches, sample_id, output, metadata,
    ))
    diagnostics.update({"method": "mapping_em", "error_rates": errors,
        "repeat_scale": repeat_scale, "unclassified_fraction": unknown_fraction,
        "coverage_gate": gate, "reference_penalties": penalties, "groups": groups})
    summary_path = output / "classification.json"
    summary_path.write_text(json.dumps(diagnostics, indent=2, allow_nan=False) + "\n")
    paths["classification_details"] = summary_path
    # Clear this run's taxonomic tables even when taxonomy is disabled, so an
    # older result in the same output directory cannot become the new report.
    write_tsv([], output / "taxonomic_identification.tsv", ["sample_id", "method"])
    write_tsv([], output / "taxonomic_identification_evidence.tsv", ["taxon_id"])
    if taxon_identification is not False:
        best = matches[0] if matches else {}
        best_support = best.get("locus_balanced_fraction", 0.0)
        best_refs = best.get("equivalent_references", "")
        best_taxa = {str(metadata.get(ref, {}).get("taxon_id") or metadata.get(ref, {}).get("taxid") or "")
                     for ref in best_refs.split(";")} if best_refs else set()
        best_taxon = next(iter(best_taxa)) if len(best_taxa) == 1 and "" not in best_taxa else ""
        # Identification follows the same reference ranking as the match table.
        # Taxon annotations never pool support or change the winning reference.
        status = "INSUFFICIENT_EVIDENCE"
        if best:
            status = "CLOSEST_REFERENCE_LOW_CONFIDENCE"
            if best_support >= 0.9 and len(detected) >= minimum_loci and diagnostics["converged"] and len(reference_loci) > 1:
                status = "SUPPORTED" if ";" not in best_refs else "AMBIGUOUS_REFERENCES"
            if sample_mode == "metagenome" and sum(row["locus_balanced_fraction"] >= 0.05 for row in matches) > 1:
                status = "MIXED_REFERENCES"
        assignment = "Unresolved"
        if best_taxon:
            assignment = best.get("taxon_name") or f"Taxon ID: {best_taxon}"
        elif best:
            assignment = "Taxonomy unavailable"
            if len(best_taxa - {""}) > 1:
                closest_metadata = metadata.get(best["reference_id"], {})
                closest_taxon = (closest_metadata.get("taxon_name") or closest_metadata.get("species")
                                 or closest_metadata.get("organism_name")
                                 or closest_metadata.get("taxon_id") or closest_metadata.get("taxid")
                                 or "Taxonomy unavailable")
                assignment = f"Ambiguous call: {closest_taxon}"
            elif best_taxa - {""}:
                assignment = "Unresolved taxonomy (incomplete reference metadata)"
        summary = {"sample_id": sample_id, "method": "mapping_em", "sample_mode": sample_mode,
            "reference_id": best.get("reference_id", ""), "equivalent_references": best_refs,
            "best_taxon": best_taxon, "best_species": best.get("taxon_name", "") if best_taxon else "",
            "assignment": assignment, "assignment_status": status,
            "assignment_rank": "taxon" if status == "SUPPORTED" and best_taxon else "unresolved",
            "confidence": "MODEL_SUPPORTED" if status == "SUPPORTED" else "LOW",
            "model_support": best_support, "unclassified_fraction": unknown_fraction,
            "informative_loci": len(detected), "expected_loci": len(loci),
            "loci_recovered": len(detected), "closest_reference": best.get("reference_id", ""),
            "status_reason": "MAPPING_LIKELIHOOD_SUPPORT" if status == "SUPPORTED" else "AMBIGUOUS_OR_INSUFFICIENT_MAPPING_EVIDENCE"}
        taxon_path = output / "taxonomic_identification.tsv"
        write_tsv([summary], taxon_path, list(summary))
        evidence_path = output / "taxonomic_identification_evidence.tsv"
        write_tsv(({"sample_id": sample_id, "rank": row["rank"],
                    "reference_id": row["reference_id"], "equivalent_references": row["equivalent_references"],
                    "taxon_id": row["taxon_id"], "species": row["taxon_name"],
                    "model_support": row["locus_balanced_fraction"]} for row in matches), evidence_path,
            ["sample_id", "rank", "reference_id", "equivalent_references", "taxon_id", "species", "model_support"])
        paths.update(taxonomic_identification=taxon_path, taxonomic_identification_evidence=evidence_path)
    return paths


def write_profile_tree(contexts, members, query_counts, matches, sample_id, output, metadata):
    """A phenetic NJ tree on one shared set of observed repeat-count loci.

    No EM weights or inferred absences enter these pairwise distances. Larger
    sample collections use the existing export-myoga shared-overlap workflow.
    """
    from .profile_tree import neighbor_joining_tree_from_matrix

    called = {}
    for locus, value in query_counts.items():
        try:
            count = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(count) and count >= 0:
            called[locus] = count
    profiles = defaultdict(dict)
    for context in contexts:
        if context.repeat_count is not None and math.isfinite(context.repeat_count):
            for ref in members[context.candidate_id]:
                profiles[ref][context.locus_id] = context.repeat_count
    shared = sorted(called)
    candidates = [ref for row in matches for ref in row["equivalent_references"].split(";")]
    eligible = [ref for ref in candidates if all(locus in profiles[ref] for locus in shared)][:20]
    status_path = output / "mlva_profile_tree_status.tsv"
    reason = "WRITTEN" if len(shared) >= 2 and eligible else "INSUFFICIENT_SHARED_REPEAT_CALLS"
    write_tsv([{"status": reason, "loci": ",".join(shared), "references": len(eligible),
        "distance": "mean_absolute_repeat_count_difference", "interpretation": "MLVA profile similarity, not an evolutionary phylogeny"}],
        status_path, ["status", "loci", "references", "distance", "interpretation"])
    paths = {"mlva_profile_tree_status": status_path}
    tree_path = output / "mlva_profiles.tree"
    if reason != "WRITTEN":
        tree_path.unlink(missing_ok=True)
        for name in ("mlva_profile_distances.tsv", "mlva_profile_tree_metadata.tsv"):
            (output / name).unlink(missing_ok=True)
        return paths
    query_label = sample_id
    while query_label in profiles:
        query_label = "QUERY__" + query_label
    labels = eligible + [query_label]
    values = np.array([[profiles[ref][locus] for locus in shared] for ref in eligible] + [[called[locus] for locus in shared]])
    distances = np.abs(values[:, None, :] - values[None, :, :]).mean(axis=2)
    tree_path.write_text(neighbor_joining_tree_from_matrix(labels, distances))
    matrix_path = output / "mlva_profile_distances.tsv"
    write_tsv((dict(zip(["sample_id", *labels], [label, *distances[i]])) for i, label in enumerate(labels)), matrix_path, ["sample_id", *labels])
    metadata_path = output / "mlva_profile_tree_metadata.tsv"
    write_tsv(({"sample_id": label, "record_type": "query" if label == query_label else "reference", **{key: metadata.get(label, {}).get(key, "") for key in ("taxon_id", "latitude", "longitude", "location")}} for label in labels), metadata_path,
        ["sample_id", "record_type", "taxon_id", "latitude", "longitude", "location"])
    paths.update(mlva_profile_tree=tree_path, mlva_profile_distances=matrix_path, mlva_profile_tree_metadata=metadata_path)
    return paths
