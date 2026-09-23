"""Canonical recovered products and technology-independent marker typing.

Recovery confidence is deliberately separate from genotype measurement. A
synthetic product cannot turn an uncertain length estimate into a certain call.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import time
from pathlib import Path

from .combined_marker_phylogeny import decompose_marker_sequence
from .io import write_fasta, write_tsv
from .models import Locus
from .repeat_calibration import assembly_equivalent_product_allele


@dataclass(frozen=True)
class LocusProduct:
    sample_id: str
    locus_id: str
    source_type: str
    variant_id: str
    sequence: str
    orientation: str = "+"
    support_count: int = 1
    effective_depth: float = 1.0
    estimated_fraction: float = 1.0
    reconstruction_confidence: float = 1.0
    calibrated_product_size_bp: int | None = None
    round_tolerance: float = 0.25
    repeat_likelihoods: dict = field(default_factory=dict)
    evidence: dict = field(default_factory=dict)


@dataclass(frozen=True)
class ProductGenotype:
    product: LocusProduct
    repeat_count_raw: float | None
    repeat_count: int | float | None
    snp_sequence: str
    haplotype_id: str
    combined_marker: str
    masking_method: str


def product_repeat_allele(product: LocusProduct, locus: Locus, tolerance: float | None = None):
    """Measure repeat count without computing an unused SNP haplotype."""
    return assembly_equivalent_product_allele(
        locus, product.calibrated_product_size_bp or len(product.sequence),
        product.round_tolerance if tolerance is None else tolerance,
    )


def genotype_product(product: LocusProduct, locus: Locus, tolerance: float | None = None) -> ProductGenotype:
    """Use the assembly calibration and existing repeat masker for every input."""
    raw, repeat = product_repeat_allele(product, locus, tolerance)
    components = decompose_marker_sequence(locus, product.sequence)
    haplotype = hashlib.sha256(components.snp_sequence.encode()).hexdigest()[:20]
    return ProductGenotype(product, raw, repeat, components.snp_sequence, haplotype,
                           f"{repeat}:{haplotype}", components.masking_method)


def assembly_product(product: dict, sample_id: str, round_tolerance: float = 0.25) -> LocusProduct:
    return LocusProduct(sample_id, product["locus_id"], "assembly", product["product_id"],
                        product.get("sequence", ""),
                        calibrated_product_size_bp=int(product["product_size_bp"]),
                        round_tolerance=round_tolerance,
                        evidence={"original_orientation": product.get("orientation", "+")})


PRODUCT_FIELDS = ["sample_id", "locus_id", "source_type", "variant_id", "orientation",
                  "support_count", "effective_depth", "estimated_fraction",
                  "reconstruction_confidence", "round_tolerance", "repeat_count", "repeat_count_raw",
                  "haplotype_id", "combined_marker", "masking_method", "snp_sequence",
                  "sequence", "repeat_likelihoods", "evidence"]


def write_products(products: list[LocusProduct], loci: list[Locus], outdir: str | Path,
                   *, progress=None, statistics=None) -> dict[str, Path]:
    output = Path(outdir)
    by_id = {l.locus_id: l for l in loci}
    genotypes = []
    started = time.perf_counter()
    for product in products:
        if progress is not None:
            progress.step(f"[{product.sample_id}] {product.locus_id}: masking {len(product.sequence):,} bases for canonical genotyping")
        genotype = genotype_product(product, by_id[product.locus_id])
        genotypes.append(genotype)
    if statistics is not None:
        statistics['genotyping_seconds'] = time.perf_counter() - started
    def rows():
        for genotype in genotypes:
            product = genotype.product
            row = {name: getattr(product, name, "") for name in PRODUCT_FIELDS}
            for name in ("repeat_count", "repeat_count_raw", "haplotype_id", "combined_marker", "masking_method", "snp_sequence"):
                row[name] = getattr(genotype, name)
            row["repeat_likelihoods"] = json.dumps(product.repeat_likelihoods, sort_keys=True)
            row["evidence"] = json.dumps(product.evidence, sort_keys=True)
            yield row
    if progress is not None:
        progress.step(f"Writing {len(products):,} canonical product sequences and evidence records")
    started = time.perf_counter()
    fasta = output / "reconstructed_loci.fasta"
    table = output / "reconstructed_locus_variants.tsv"
    write_fasta(((p.variant_id, p.sequence) for p in products), fasta)
    write_tsv(rows(), table, PRODUCT_FIELDS)
    if statistics is not None:
        statistics['product_tables_seconds'] = time.perf_counter() - started
    started = time.perf_counter()
    alignment_paths = write_product_alignments(genotypes, output, progress=progress)
    if statistics is not None:
        statistics['canonical_alignments_seconds'] = time.perf_counter() - started
    return {"reconstructed_loci": fasta, "reconstructed_locus_variants": table,
            **alignment_paths}


def write_product_alignments(genotypes: list[ProductGenotype], output: Path, *, progress=None) -> dict[str, Path]:
    """Align repeat-masked variants once, with explicit canonical coordinates.

    SNP depths are component effective molecule counts, not pileup read counts.
    Unknown bases do not contribute. Reference coordinates are 1-based positions
    in the dominant product's repeat-masked sequence.
    """
    from collections import defaultdict, Counter
    import parasail
    from .mapping import SNP_FIELDS, MAPPING_SUMMARY_FIELDS

    by_locus = defaultdict(list)
    for genotype in genotypes:
        by_locus[genotype.product.locus_id].append(genotype)
    snps, summaries, alignments, references = [], [], [], []
    for locus_id, variants in sorted(by_locus.items()):
        if progress is not None:
            progress.step(f"[{variants[0].product.sample_id}] {locus_id}: aligning {len(variants):,} canonical SNP sequences")
        variants.sort(key=lambda g: (-g.product.estimated_fraction, -g.product.support_count, g.product.variant_id))
        dominant = variants[0]
        reference = dominant.snp_sequence
        references.append((dominant.product.variant_id, reference))
        counts = defaultdict(Counter)
        ref_bases = {}
        matrix = parasail.matrix_create("ACGTN", 2, -4)
        for genotype in variants:
            if not reference or not genotype.snp_sequence:
                continue
            if genotype.snp_sequence == reference and not set(reference) - set('ACGTN'):
                # Identical strings attain the maximum score without any DP.
                query = target = reference
                cigar, score = f"{len(reference)}=", 2 * len(reference)
            else:
                result = parasail.nw_trace_striped_32(genotype.snp_sequence, reference, 5, 1, matrix)
                query, target = result.traceback.query, result.traceback.ref
                cigar, score = result.cigar.decode.decode(), result.score
            position = 0
            for base, ref in zip(query, target):
                if ref != "-":
                    position += 1
                    ref_bases[position] = ref
                    if base in "ACGT-" and ref in "ACGT":
                        counts[position][base] += genotype.product.effective_depth
            alignments.append({"sample_id": genotype.product.sample_id, "locus_id": locus_id,
                "variant_id": genotype.product.variant_id, "reference_variant_id": dominant.product.variant_id,
                "aligned_query": query, "aligned_reference": target,
                "cigar": cigar, "alignment_score": score})
        n_snps = 0
        for position, alleles in sorted(counts.items()):
            total = sum(alleles.values())
            for base, depth in sorted(alleles.items()):
                if base == ref_bases[position] or depth <= 0:
                    continue
                n_snps += 1
                snps.append({"sample_id": dominant.product.sample_id, "locus_id": locus_id,
                    "reference_variant_id": dominant.product.variant_id, "position": position,
                    "reference_base": ref_bases[position], "alternate_base": base,
                    "depth": total, "alternate_depth": depth, "alternate_frequency": depth/total,
                    "coordinate_system": "repeat_masked_product", "method": "canonical_variant_alignment"})
        summaries.append({"sample_id": dominant.product.sample_id, "locus_id": locus_id,
            "reference_variant_id": dominant.product.variant_id, "reference_length_bp": len(reference),
            "total_reads": sum(g.product.support_count for g in variants),
            "mean_depth": sum(sum(v.values()) for v in counts.values())/max(len(reference),1),
            "covered_bases": len(counts), "coverage_percent": 100*len(counts)/max(len(reference),1),
            "snp_count": n_snps, "method": "canonical_variant_alignment"})
    paths = {"mapping_snps": output / "locus_snps.tsv",
             "mapping_summary": output / "locus_mapping_summary.tsv",
             "mapping_references": output / "locus_mapping_references.fasta.gz",
             "canonical_alignments": output / "canonical_locus_alignments.tsv"}
    write_tsv(snps, paths["mapping_snps"], SNP_FIELDS + ["coordinate_system", "method"])
    write_tsv(summaries, paths["mapping_summary"], MAPPING_SUMMARY_FIELDS + ["method"])
    write_tsv(alignments, paths["canonical_alignments"], ["sample_id", "locus_id", "variant_id", "reference_variant_id",
                                                       "aligned_query", "aligned_reference", "cigar", "alignment_score"])
    write_fasta(references, paths["mapping_references"])
    return paths
