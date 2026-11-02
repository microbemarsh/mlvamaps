"""Supplementary genotype uncertainty and observed sequence alleles."""
from collections import defaultdict
import hashlib
import json

from .io import write_tsv
from .locus_measurement import measure_locus_product
from .sequence import repeat_motif_statistics


GENOTYPE_FIELDS = ["sample", "locus", "genotype_status", "repeat_count", "repeat_lower_bound",
    "repeat_interval_low", "repeat_interval_high", "independent_support", "full_span_support",
    "left_anchor_support", "right_anchor_support", "ambiguous_locus_support", "forward_support",
    "reverse_support", "inference_tier", "candidate_distribution", "log_likelihood_gap",
    "computational_fallback", "reason", "error_model"]
SEQUENCE_FIELDS = ["sample", "locus", "sequence_id", "independent_support", "repeat_count",
    "raw_repeat_count", "repeat_sequence", "repeat_pattern", "partial_terminal_bases",
    "motif_mismatches", "product_sequence"]


def write_evidence_reports(calls, evidence, loci, outdir, thresholds):
    grouped = defaultdict(list)
    for row in evidence:
        grouped[row.locus_id].append(row)
    locus_by_id = {locus.locus_id: locus for locus in loci}
    summaries, sequences = [], []
    for call in calls:
        rows = grouped[call["locus"]]
        def support(predicate):
            return len({row.molecule_id for row in rows if predicate(row)})
        ambiguous = support(lambda row: "ambiguous_locus" in row.quality_flags)
        measured = [row for row in rows if row.measured_repeat_count is not None]
        direct = support(lambda row: row.direct_product_measurement and not row.quality_flags)
        fallback = call.get("computational_fallback", "")
        reason = call.get("estimation_reason", "")
        if call["status"] == "mixed" or call.get("possible_secondary"):
            status, reason = "POSSIBLE_MIXTURE", "Independent secondary-allele evidence; strain mixture is not confirmed."
        elif fallback:
            status, reason = "UNRESOLVED", fallback
        elif not call["molecule_support"] and ambiguous:
            status, reason = "AMBIGUOUS_LOCUS", "Competing locus assignments; these molecules do not vote for an allele."
        elif call["repeat_count"] == "" and call.get("repeat_lower_bound"):
            status, reason = "LOWER_BOUND", "Anchored terminal repeat sequence; the opposite boundary is unobserved."
        elif call["repeat_count"] == "":
            status, reason = "UNRESOLVED", reason or ("No recruited evidence." if not rows else "Insufficient identifiable allele evidence.")
        elif (call["status"] == "called" and direct >= max(3, thresholds.minimum_molecules)
              and measured and all(row.measured_repeat_count == call["repeat_count"]
                  and not row.quality_flags and (row.repeat_sigma or 1) <= .25
                  and not set(row.metadata.get("product_sequence", "")) - set("ACGT") for row in measured)
              and not call.get("derived_molecule_support")):
            status, reason = "EXACT", "Independent unambiguous complete products agree with the reported allele."
        else:
            status = "ESTIMATED"
            reason = reason or "Coverage, measurement uncertainty, or reconstruction does not meet EXACT criteria."
        call["genotype_status"] = status
        call["inference_tier"] = max([call.get("inference_tier", 0)] + [row.metadata.get("inference_tier", 0) for row in rows])
        summaries.append(dict(sample=call["sample"], locus=call["locus"], genotype_status=status,
            repeat_count=call["repeat_count"], repeat_lower_bound=call.get("repeat_lower_bound", ""),
            repeat_interval_low=call.get("repeat_interval_low", ""), repeat_interval_high=call.get("repeat_interval_high", ""),
            independent_support=call["molecule_support"], full_span_support=call["full_span_support"],
            left_anchor_support=support(lambda row: row.left_boundary_span or row.metadata.get("left_anchor")),
            right_anchor_support=support(lambda row: row.right_boundary_span or row.metadata.get("right_anchor")),
            ambiguous_locus_support=ambiguous,
            forward_support=support(lambda row: "forward" in row.orientation),
            reverse_support=support(lambda row: "reverse" in row.orientation),
            inference_tier=max([call.get("inference_tier", 0)] + [row.metadata.get("inference_tier", 0) for row in rows]),
            candidate_distribution=call["candidate_distribution"], log_likelihood_gap=call.get("log_likelihood_gap", ""),
            computational_fallback=fallback, reason=reason,
            error_model=json.dumps(next((row.metadata["error_model"] for row in rows if "error_model" in row.metadata), {}), sort_keys=True)))
        products = defaultdict(set)
        ambiguous_molecules = {row.molecule_id for row in rows if "ambiguous_locus" in row.quality_flags}
        for row in measured:
            product = row.metadata.get("product_sequence")
            if (product and not set(product) - set("ACGT") and row.direct_product_measurement
                    and not row.quality_flags and row.molecule_id not in ambiguous_molecules):
                products[row.molecule_id].add(product)
        groups = defaultdict(set)
        for molecule, products_for_molecule in products.items():
            if len(products_for_molecule) == 1:
                groups[next(iter(products_for_molecule))].add(molecule)
        locus = locus_by_id[call["locus"]]
        for product, molecules in sorted(groups.items()):
            if len(molecules) < max(3, thresholds.minimum_molecules):
                continue
            measurement = measure_locus_product(product, locus)
            if not (measurement.metadata or {}).get("flanks_resolved") or (measurement.metadata or {}).get("boundary_ambiguous"):
                continue
            repeat = measurement.repeat_sequence
            unit = locus.repeat_unit_length_bp or len(locus.repeat_motif)
            pattern, mismatches, _ = repeat_motif_statistics(repeat, locus.repeat_motif, unit) if unit and locus.repeat_motif else ([], 0, 0)
            sequences.append(dict(sample=call["sample"], locus=call["locus"],
                sequence_id="sha256:" + hashlib.sha256(product.encode("ascii")).hexdigest(),
                independent_support=len(molecules), repeat_count=measurement.called_allele,
                raw_repeat_count=measurement.raw_repeat_count, repeat_sequence=repeat,
                repeat_pattern="-".join(pattern), partial_terminal_bases=len(repeat) % unit if unit else "",
                motif_mismatches=mismatches, product_sequence=product))
    paths = {"genotype_evidence": outdir / "genotype_evidence.tsv", "sequence_alleles": outdir / "sequence_alleles.tsv"}
    write_tsv(summaries, paths["genotype_evidence"], GENOTYPE_FIELDS)
    write_tsv(sequences, paths["sequence_alleles"], SEQUENCE_FIELDS)
    return paths
