from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

import spoars

from .in_silico_pcr import read_pcr_results, run_in_silico_pcr_loci
from .io import write_fasta
from .models import Locus, RepeatFeature


LOCAL_ASSEMBLY_FIELDS = [
    "sample_id",
    "locus_id",
    "dominant_variant_id",
    "input_reads",
    "unique_sequences",
    "observed_min_product_bp",
    "observed_modal_product_bp",
    "observed_max_product_bp",
    "poa_consensus_bp",
    "pcr_product_size_bp",
    "raw_repeat_count",
    "called_repeat_count",
    "pcr_status",
    "measurement_source",
]


def _dominant_variants(mixture_rows: list[dict]) -> dict[str, str]:
    dominant: dict[str, str] = {}
    for row in sorted(
        mixture_rows,
        key=lambda item: (
            str(item.get("locus_id", "")),
            -float(item.get("estimated_fraction") or 0),
            str(item.get("variant_id", "")),
        ),
    ):
        dominant.setdefault(str(row["locus_id"]), str(row["variant_id"]))
    return dominant


def _poa_consensus(features: list[RepeatFeature]) -> str:
    graph = spoars.Poa(
        alignment_type="global",
        scoring=spoars.Scoring.default(),
    )
    # POA construction can be order-sensitive. Start with the highest-quality
    # observations and retain a deterministic read-ID tie break.
    for feature in sorted(
        features,
        key=lambda item: (-item.mean_qscore, item.read_id),
    ):
        graph.add(feature.amplicon_sequence.upper())
    return str(graph.consensus()).upper()


def assemble_dominant_locus_products(
    features: list[RepeatFeature],
    memberships: list[dict],
    mixture_rows: list[dict],
    loci: list[Locus],
    outdir: str | Path,
    sample_id: str,
    max_primer_mismatches: int,
    round_tolerance: float,
    threads: int,
) -> tuple[list[tuple[str, str]], dict[str, dict], list[dict], dict[str, Path]]:
    """Assemble dominant locus clusters and call them through assembly logic."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    consensus_path = outdir.parent / "local_locus_products.fasta"
    dominant_by_locus = _dominant_variants(mixture_rows)
    variant_by_read = {
        (str(row["locus_id"]), str(row["read_id"])): str(row["variant_id"])
        for row in memberships
    }
    by_locus: dict[str, list[RepeatFeature]] = defaultdict(list)
    for feature in features:
        if (
            feature.amplicon_sequence
            and feature.product_size_bp > 0
            and variant_by_read.get((feature.locus_id, feature.read_id))
            == dominant_by_locus.get(feature.locus_id)
        ):
            by_locus[feature.locus_id].append(feature)

    records: list[tuple[str, str]] = []
    audit_by_locus: dict[str, dict] = {
        locus.locus_id: {
            "sample_id": sample_id,
            "locus_id": locus.locus_id,
            "dominant_variant_id": dominant_by_locus.get(locus.locus_id, ""),
            "input_reads": 0,
            "unique_sequences": 0,
            "observed_min_product_bp": "",
            "observed_modal_product_bp": "",
            "observed_max_product_bp": "",
            "poa_consensus_bp": "",
            "pcr_product_size_bp": "",
            "raw_repeat_count": "",
            "called_repeat_count": "",
            "pcr_status": "NO_DOMINANT_COMPLETE_PRODUCT",
            "measurement_source": "per_read_fallback",
        }
        for locus in loci
    }
    contig_locus: dict[str, str] = {}
    for index, (locus_id, locus_features) in enumerate(
        sorted(by_locus.items())
    ):
        consensus = _poa_consensus(locus_features)
        if not consensus:
            continue
        contig_id = f"local_poa_{index:06d}"
        contig_locus[contig_id] = locus_id
        records.append((contig_id, consensus))
        lengths = Counter(feature.product_size_bp for feature in locus_features)
        modal_length = min(
            lengths,
            key=lambda length: (-lengths[length], length),
        )
        audit_by_locus[locus_id].update(
            {
                "dominant_variant_id": dominant_by_locus[locus_id],
                "input_reads": len(locus_features),
                "unique_sequences": len(
                    {feature.amplicon_sequence for feature in locus_features}
                ),
                "observed_min_product_bp": min(lengths),
                "observed_modal_product_bp": modal_length,
                "observed_max_product_bp": max(lengths),
                "poa_consensus_bp": len(consensus),
                "pcr_status": "POA_PCR_NOT_FOUND",
            }
        )

    write_fasta(records, consensus_path)
    measurements: dict[str, dict] = {}
    # Imports are intentionally local: assembly_call imports shared output
    # helpers from pipeline during module initialization.
    from .assembly_call import (
        legacy_amplicon_bounds,
        legacy_assembly_call_rows,
        pcr_rows_to_products,
    )

    pcr_paths = run_in_silico_pcr_loci(
        consensus_path,
        loci,
        outdir,
        max_errors=max_primer_mismatches,
        threads=threads,
        max_n_fraction=1.0,
        amplicon_bounds=legacy_amplicon_bounds(loci),
    )
    products = pcr_rows_to_products(
        read_pcr_results(pcr_paths["stats"], pcr_paths["products"]),
        loci,
        sample_id,
        enforce_locus_bounds=False,
        reference_order={
            contig_id: index for index, contig_id in enumerate(contig_locus)
        },
    )
    # A locally assembled contig is evidence only for the cluster from which
    # it was built, even if a short primer cross-matches elsewhere.
    products = [
        product
        for product in products
        if contig_locus.get(str(product["contig"]))
        == str(product["locus_id"])
    ]
    product_sizes_by_locus: dict[str, set[int]] = defaultdict(set)
    for product in products:
        product_sizes_by_locus[str(product["locus_id"])].add(
            int(product["product_size_bp"])
        )
    for locus_id, product_sizes in product_sizes_by_locus.items():
        audit_by_locus[locus_id].update(
            {
                "pcr_product_size_bp": ",".join(
                    str(size) for size in sorted(product_sizes)
                ),
                "pcr_status": "PCR_PRODUCT_UNCALLABLE",
            }
        )
    product_by_id = {
        str(product["product_id"]): product for product in products
    }
    call_rows = legacy_assembly_call_rows(
        loci,
        products,
        sample_id,
        round_tolerance=round_tolerance,
    )
    for call in call_rows:
        locus_id = str(call["locus_id"])
        if call.get("present") != "yes":
            continue
        product = product_by_id.get(str(call.get("evidence", "")))
        if product is None:
            continue
        raw_count = call.get("repeat_count_raw", "")
        called_count = call.get("repeat_count", "")
        measurements[locus_id] = {
            "product_size_bp": int(product["product_size_bp"]),
            "product_sequence": str(product["sequence"]),
            "raw_repeat_count": (
                float(raw_count) if raw_count not in ("", None) else ""
            ),
            "called_repeat_count": called_count,
            "supporting_reads": audit_by_locus[locus_id]["input_reads"],
            "source": "dominant_cluster_poa_assembly",
        }
        audit_by_locus[locus_id].update(
            {
                "pcr_product_size_bp": int(product["product_size_bp"]),
                "raw_repeat_count": raw_count,
                "called_repeat_count": called_count,
                "pcr_status": "PASS",
                "measurement_source": "dominant_cluster_poa_assembly",
            }
        )

    display_records = [
        (f"{contig_locus[contig_id]}_local_primary", sequence)
        for contig_id, sequence in records
    ]
    return (
        display_records,
        measurements,
        [
            audit_by_locus[locus_id]
            for locus_id in sorted(audit_by_locus)
        ],
        pcr_paths,
    )
