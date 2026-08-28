from __future__ import annotations

import csv
from dataclasses import asdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from .assembly_call import pcr_rows_to_products
from .concurrency import DEFAULT_THREADS, resolve_threads
from .in_silico_pcr import read_pcr_results, run_in_silico_pcr_loci
from .io import open_text
from .models import Locus
from .phylogeny import (
    REFERENCE_ASSEMBLY_FIELDS,
    build_reference_phylogenies,
    canonical_assembly_digest,
)
from .primers import read_loci_or_primers
from .progress import ProgressReporter


_ASSEMBLY_SUFFIXES = (
    ".fasta.gz",
    ".fna.gz",
    ".fa.gz",
    ".fas.gz",
    ".fasta",
    ".fna",
    ".fa",
    ".fas",
)
_ID_ALIASES = (
    "reference_id",
    "genome_id",
    "sample_id",
    "isolate_id",
    "strain",
    "accession",
    "assembly_accession",
    "name",
    "id",
)
_ASSEMBLY_ALIASES = ("assembly_file", "assembly", "filename", "file", "path")

REFERENCE_BUILD_FIELDS = [
    "reference_id",
    "assembly_file",
    "locus_id",
    "status",
    "product_count",
    "best_product_count",
    "selected_product_id",
    "product_size_bp",
    "forward_mismatches",
    "reverse_mismatches",
]

REFERENCE_BUILD_STATUS_BUILT = "BUILT"
REFERENCE_BUILD_STATUS_PARTIAL = "PARTIAL"
REFERENCE_BUILD_STATUS_NO_USABLE_LOCI = "NO_USABLE_LOCI"
REFERENCE_BUILD_STATUS_INSUFFICIENT_REFERENCES = "INSUFFICIENT_REFERENCES"

REFERENCE_LOCUS_AMPLIFIABILITY_FIELDS = [
    "locus_id",
    "genomes_examined",
    "genomes_with_valid_amplicon",
    "valid_amplicons",
    "percent_genomes_amplifiable",
    "amplifiable",
    "tree_status",
]


def _strip_assembly_suffix(name: str) -> str:
    lower = name.lower()
    for suffix in _ASSEMBLY_SUFFIXES:
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _read_metadata(path: str | Path) -> tuple[list[str], list[dict[str, str]], str, str | None]:
    metadata_path = Path(path)
    with metadata_path.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        first_line = sample.splitlines()[0] if sample.splitlines() else ""
        delimiter = "\t" if "\t" in first_line else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Metadata has no header: {metadata_path}")
        fields = [str(field).strip() for field in reader.fieldnames]
        rows = [
            {str(key).strip(): str(value or "").strip() for key, value in row.items()}
            for row in reader
        ]
    lookup = {field.lower(): field for field in fields}
    id_field = next((lookup[name] for name in _ID_ALIASES if name in lookup), None)
    if id_field is None:
        raise ValueError(
            "Metadata needs an identifier column such as reference_id, genome_id, "
            "sample_id, strain, accession, or id"
        )
    assembly_field = next(
        (lookup[name] for name in _ASSEMBLY_ALIASES if name in lookup), None
    )
    return fields, rows, id_field, assembly_field


def _discover_assemblies(directory: str | Path) -> list[Path]:
    root = Path(directory)
    if not root.is_dir():
        raise ValueError(f"Assemblies path is not a directory: {root}")
    paths = [
        path
        for path in root.rglob("*")
        if path.is_file() and path.name.lower().endswith(_ASSEMBLY_SUFFIXES)
    ]
    if not paths:
        raise ValueError(f"No FASTA assemblies found under {root}")
    return sorted(paths)


def _match_assemblies(
    assembly_dir: str | Path,
    metadata_rows: list[dict[str, str]],
    id_field: str,
    assembly_field: str | None,
) -> list[tuple[str, Path, dict[str, str]]]:
    root = Path(assembly_dir)
    paths = _discover_assemblies(root)
    by_name: dict[str, list[Path]] = {}
    by_stem: dict[str, list[Path]] = {}
    for path in paths:
        by_name.setdefault(path.name, []).append(path)
        by_stem.setdefault(_strip_assembly_suffix(path.name), []).append(path)

    matched: list[tuple[str, Path, dict[str, str]]] = []
    seen_ids: set[str] = set()
    used_paths: set[Path] = set()
    for row in metadata_rows:
        reference_id = row.get(id_field, "").strip()
        if not reference_id:
            continue
        if any(character.isspace() for character in reference_id):
            raise ValueError(
                f"Reference id {reference_id!r} contains whitespace; FASTA/tree identifiers "
                "must be single tokens"
            )
        if reference_id in seen_ids:
            raise ValueError(f"Duplicate metadata identifier {reference_id!r}")
        seen_ids.add(reference_id)
        requested = row.get(assembly_field, "").strip() if assembly_field else ""
        candidates: list[Path]
        if requested:
            direct = Path(requested)
            if not direct.is_absolute():
                direct = root / direct
            if direct.is_file():
                candidates = [direct]
            else:
                candidates = by_name.get(Path(requested).name, [])
        else:
            candidates = by_stem.get(reference_id, [])
        if not candidates:
            target = requested or reference_id
            raise ValueError(f"No assembly matched metadata row {reference_id!r} ({target!r})")
        if len(candidates) > 1:
            raise ValueError(
                f"Multiple assemblies matched metadata row {reference_id!r}: "
                + ", ".join(str(path) for path in candidates)
            )
        assembly = candidates[0].resolve()
        if assembly in used_paths:
            raise ValueError(f"Assembly {assembly} was assigned to more than one reference")
        used_paths.add(assembly)
        matched.append((reference_id, assembly, row))
    if not matched:
        raise ValueError("Metadata contains no non-empty reference identifiers")
    return matched


def _write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


def _write_tsv(rows: list[dict], path: Path, fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _locus_amplifiability_rows(
    loci: list[Locus],
    records_by_locus: dict[str, list[tuple[str, str]]],
    genomes_examined: int,
    min_references_per_tree: int,
) -> list[dict[str, Any]]:
    """Summarize products retained by the existing extraction/QC policy."""
    rows = []
    for locus in loci:
        records = records_by_locus[locus.locus_id]
        genomes_with_amplicon = len({reference_id for reference_id, _sequence in records})
        valid_amplicons = len(records)
        if valid_amplicons == 0:
            tree_status = "NO_AMPLICONS"
        elif valid_amplicons < min_references_per_tree:
            tree_status = REFERENCE_BUILD_STATUS_INSUFFICIENT_REFERENCES
        else:
            tree_status = REFERENCE_BUILD_STATUS_BUILT
        rows.append(
            {
                "locus_id": locus.locus_id,
                "genomes_examined": genomes_examined,
                "genomes_with_valid_amplicon": genomes_with_amplicon,
                "valid_amplicons": valid_amplicons,
                "percent_genomes_amplifiable": (
                    round(100.0 * genomes_with_amplicon / genomes_examined, 1)
                    if genomes_examined
                    else 0.0
                ),
                "amplifiable": "TRUE" if valid_amplicons else "FALSE",
                "tree_status": tree_status,
            }
        )
    return rows


def _product_rank(product: dict) -> tuple:
    return (
        int(product["primer_error_round"]),
        int(product["forward_mismatches"]) + int(product["reverse_mismatches"]),
        int(product["product_size_bp"]),
        str(product["product_id"]),
    )


def _primer_quality(product: dict) -> tuple[int, int]:
    return (
        int(product["primer_error_round"]),
        int(product["forward_mismatches"]) + int(product["reverse_mismatches"]),
    )


def _extract_reference_products(
    reference_id: str,
    assembly: Path,
    loci: list[Locus],
    extraction_dir: Path,
    max_primer_mismatches: int,
) -> tuple[str, list[dict]]:
    """Extract one assembly in a worker process and return normalized products."""
    paths = run_in_silico_pcr_loci(
        assembly,
        loci,
        extraction_dir / reference_id,
        max_errors=max_primer_mismatches,
        threads=1,
    )
    products = pcr_rows_to_products(
        read_pcr_results(paths["stats"], paths["products"]), loci, reference_id
    )
    return reference_id, products


def build_reference_database(
    assemblies_dir: str | Path,
    primers_path: str | Path,
    metadata_path: str | Path,
    outdir: str | Path,
    *,
    loci_path: str | Path | None = None,
    multiple_products: str = "exclude",
    max_primer_mismatches: int = 2,
    min_references_per_tree: int = 3,
    threads: int = DEFAULT_THREADS,
    amplirust_bin: str = "amplirust",
    mafft_bin: str = "mafft",
    raxml_ng_bin: str = "raxml-ng",
    raxml_model: str = "DNA",
    show_progress: bool = False,
) -> dict[str, Any]:
    """Extract reference amplicons and infer a SNP tree for every usable locus."""
    if multiple_products not in {"exclude", "best", "error"}:
        raise ValueError("multiple_products must be exclude, best, or error")
    if min_references_per_tree < 2:
        raise ValueError("min_references_per_tree must be at least 2")
    loci: list[Locus] = read_loci_or_primers(loci_path, None if loci_path else primers_path)
    if not loci:
        raise ValueError("Primer/locus panel contains no loci")
    unsafe_loci = [locus.locus_id for locus in loci if Path(locus.locus_id).name != locus.locus_id]
    if unsafe_loci:
        raise ValueError(f"Locus identifiers cannot contain path separators: {unsafe_loci}")

    metadata_fields, metadata_rows, id_field, assembly_field = _read_metadata(metadata_path)
    metadata_lookup = {field.lower(): field for field in metadata_fields}
    taxon_source = next(
        (metadata_lookup[name] for name in ("taxon_id", "taxid", "ncbi_taxid") if name in metadata_lookup),
        None,
    )
    name_source = next(
        (metadata_lookup[name] for name in ("taxon_name", "species", "organism_name", "scientific_name") if name in metadata_lookup),
        None,
    )
    if taxon_source:
        missing_taxa = [str(row.get(id_field, "")) for row in metadata_rows if not str(row.get(taxon_source, "")).strip()]
        if missing_taxa:
            raise ValueError(
                "Taxonomic metadata is incomplete; blank taxon identifiers for: "
                + ", ".join(missing_taxa[:10])
            )
    matched = _match_assemblies(assemblies_dir, metadata_rows, id_field, assembly_field)
    thread_count = resolve_threads(threads)
    progress = ProgressReporter(enabled=show_progress)
    progress.step(
        f"Starting reference build for {len(matched):,} assemblies and "
        f"{len(loci):,} loci with {thread_count} worker(s)"
    )
    output = Path(outdir)
    database_dir = output / "database"
    extraction_dir = output / "extraction"
    database_dir.mkdir(parents=True, exist_ok=True)
    records_by_locus: dict[str, list[tuple[str, str]]] = {
        locus.locus_id: [] for locus in loci
    }
    manifest_rows: list[dict] = []

    progress.step("Extracting primer products from reference assemblies")
    products_by_reference: dict[str, list[dict]] = {}
    if thread_count == 1 or len(matched) == 1:
        for completed, (reference_id, assembly, _metadata) in enumerate(matched, start=1):
            _, products = _extract_reference_products(
                reference_id,
                assembly,
                loci,
                extraction_dir,
                max_primer_mismatches,
            )
            products_by_reference[reference_id] = products
            progress.count("Extracted assemblies", completed, len(matched), force=True)
    else:
        worker_count = min(thread_count, len(matched))
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            futures = {
                executor.submit(
                    _extract_reference_products,
                    reference_id,
                    assembly,
                    loci,
                    extraction_dir,
                    max_primer_mismatches,
                ): reference_id
                for reference_id, assembly, _metadata in matched
            }
            for completed, future in enumerate(as_completed(futures), start=1):
                reference_id = futures[future]
                try:
                    _, products = future.result()
                except Exception as exc:
                    raise RuntimeError(
                        f"Reference extraction failed for {reference_id!r}"
                    ) from exc
                products_by_reference[reference_id] = products
                progress.count("Extracted assemblies", completed, len(matched), force=True)

    # Consume results in metadata order so FASTA and manifest output stay deterministic.
    for reference_id, assembly, _metadata in matched:
        products = products_by_reference[reference_id]
        products_by_locus: dict[str, list[dict]] = {}
        for product in products:
            products_by_locus.setdefault(str(product["locus_id"]), []).append(product)
        for locus in loci:
            candidates = sorted(products_by_locus.get(locus.locus_id, []), key=_product_rank)
            best_candidates = (
                [candidate for candidate in candidates if _primer_quality(candidate) == _primer_quality(candidates[0])]
                if candidates
                else []
            )
            selected: dict | None = None
            if not candidates:
                status = "NOT_FOUND"
            elif len(best_candidates) == 1:
                status = "INCLUDED"
                selected = best_candidates[0]
            elif multiple_products == "error":
                raise ValueError(
                    f"Reference {reference_id!r}, locus {locus.locus_id!r} has "
                    f"{len(best_candidates)} equally good primer products"
                )
            elif multiple_products == "best":
                status = "INCLUDED_BEST_OF_MULTIPLE"
                selected = best_candidates[0]
            else:
                status = "AMBIGUOUS_EXCLUDED"
            if selected is not None:
                records_by_locus[locus.locus_id].append(
                    (reference_id, str(selected["sequence"]).upper())
                )
            manifest_rows.append(
                {
                    "reference_id": reference_id,
                    "assembly_file": str(assembly),
                    "locus_id": locus.locus_id,
                    "status": status,
                    "product_count": len(candidates),
                    "best_product_count": len(best_candidates),
                    "selected_product_id": "" if selected is None else selected["product_id"],
                    "product_size_bp": "" if selected is None else selected["product_size_bp"],
                    "forward_mismatches": "" if selected is None else selected["forward_mismatches"],
                    "reverse_mismatches": "" if selected is None else selected["reverse_mismatches"],
                }
            )

    locus_fasta_paths: list[Path] = []
    for locus in loci:
        records = records_by_locus[locus.locus_id]
        fasta_path = database_dir / f"{locus.locus_id}.fasta.gz"
        (database_dir / f"{locus.locus_id}.fasta").unlink(missing_ok=True)
        fasta_path.unlink(missing_ok=True)
        if records:
            _write_fasta(records, fasta_path)
            locus_fasta_paths.append(fasta_path)

    panel_path = database_dir / "reference_panel.tsv"
    panel_fields = list(asdict(loci[0]))
    _write_tsv([asdict(locus) for locus in loci], panel_path, panel_fields)

    normalized_metadata = []
    myoga_metadata = []
    for reference_id, _assembly, row in matched:
        normalized = {"reference_id": reference_id}
        normalized.update({field: row.get(field, "") for field in metadata_fields if field != id_field})
        if taxon_source:
            normalized["taxon_id"] = str(row.get(taxon_source, "")).strip()
        if name_source:
            normalized["taxon_name"] = str(row.get(name_source, "")).strip()
        normalized_metadata.append(normalized)
        myoga = {"genome_id": reference_id}
        myoga.update({field: row.get(field, "") for field in metadata_fields if field != id_field})
        myoga_metadata.append(myoga)
    reference_metadata_path = database_dir / "reference_metadata.tsv"
    normalized_fields = ["reference_id", *[field for field in metadata_fields if field != id_field]]
    for field in ("taxon_id", "taxon_name"):
        if field in normalized_metadata[0] and field not in normalized_fields:
            normalized_fields.append(field)
    _write_tsv(normalized_metadata, reference_metadata_path, normalized_fields)
    reference_assemblies_path = database_dir / "reference_assemblies.tsv"
    _write_tsv(
        [
            {
                "reference_id": reference_id,
                "assembly_file": str(assembly.resolve()),
                "assembly_sha256": canonical_assembly_digest(assembly),
            }
            for reference_id, assembly, _row in matched
        ],
        reference_assemblies_path,
        REFERENCE_ASSEMBLY_FIELDS,
    )
    _write_tsv(manifest_rows, output / "reference_build_manifest.tsv", REFERENCE_BUILD_FIELDS)
    locus_summary_rows = _locus_amplifiability_rows(
        loci, records_by_locus, len(matched), min_references_per_tree
    )
    locus_summary_path = output / "reference_locus_amplifiability.tsv"
    _write_tsv(
        locus_summary_rows,
        locus_summary_path,
        REFERENCE_LOCUS_AMPLIFIABILITY_FIELDS,
    )
    with (output / "myoga_metadata.csv").open("w", newline="") as handle:
        myoga_fields = ["genome_id", *[field for field in metadata_fields if field != id_field]]
        writer = csv.DictWriter(handle, fieldnames=myoga_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(myoga_metadata)

    tree_paths: dict[str, Path | None] = {}
    if locus_fasta_paths:
        progress.step("Building per-locus reference alignments and trees")
        tree_paths = build_reference_phylogenies(
            database_dir,
            output / "phylogeny",
            loci,
            thread_count,
            min_references=min_references_per_tree,
            mafft_bin=mafft_bin,
            raxml_ng_bin=raxml_ng_bin,
            raxml_model=raxml_model,
            progress=progress,
        )
        build_status = (
            REFERENCE_BUILD_STATUS_BUILT
            if all(row["tree_status"] == REFERENCE_BUILD_STATUS_BUILT for row in locus_summary_rows)
            else REFERENCE_BUILD_STATUS_PARTIAL
        )
    else:
        build_status = REFERENCE_BUILD_STATUS_NO_USABLE_LOCI
        tree_paths = {"phylogeny": None}
        progress.step(
            "No usable reference amplicons were recovered; skipping phylogeny construction."
        )
    progress.step(f"Done. Reference database: {database_dir}")
    return {
        "status": build_status,
        "outdir": output,
        "database": database_dir,
        "metadata": reference_metadata_path,
        "myoga_metadata": output / "myoga_metadata.csv",
        "manifest": output / "reference_build_manifest.tsv",
        "locus_amplifiability": locus_summary_path,
        "locus_summary_rows": locus_summary_rows,
        "reference_assemblies": reference_assemblies_path,
        **tree_paths,
    }
