from __future__ import annotations

import csv
from pathlib import Path

from .assembly_call import amplirust_rows_to_products
from .concurrency import DEFAULT_THREADS, resolve_threads
from .in_silico_pcr import read_amplirust_results, run_amplirust_loci
from .models import Locus
from .phylogeny import build_reference_phylogenies
from .primers import read_loci_or_primers


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
    with path.open("w") as handle:
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
    raxml_model: str = "GTR+G",
) -> dict[str, Path]:
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
    matched = _match_assemblies(assemblies_dir, metadata_rows, id_field, assembly_field)
    output = Path(outdir)
    database_dir = output / "database"
    extraction_dir = output / "extraction"
    database_dir.mkdir(parents=True, exist_ok=True)
    records_by_locus: dict[str, list[tuple[str, str]]] = {
        locus.locus_id: [] for locus in loci
    }
    manifest_rows: list[dict] = []

    for reference_id, assembly, _metadata in matched:
        paths = run_amplirust_loci(
            assembly,
            loci,
            extraction_dir / reference_id,
            max_errors=max_primer_mismatches,
            threads=threads,
            executable=amplirust_bin,
        )
        products = amplirust_rows_to_products(
            read_amplirust_results(paths["stats"], paths["products"]), loci, reference_id
        )
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

    for locus in loci:
        records = records_by_locus[locus.locus_id]
        fasta_path = database_dir / f"{locus.locus_id}.fasta"
        fasta_path.unlink(missing_ok=True)
        if records:
            _write_fasta(records, fasta_path)

    normalized_metadata = []
    myoga_metadata = []
    for reference_id, _assembly, row in matched:
        normalized = {"reference_id": reference_id}
        normalized.update({field: row.get(field, "") for field in metadata_fields if field != id_field})
        normalized_metadata.append(normalized)
        myoga = {"genome_id": reference_id}
        myoga.update({field: row.get(field, "") for field in metadata_fields if field != id_field})
        myoga_metadata.append(myoga)
    reference_metadata_path = database_dir / "reference_metadata.tsv"
    normalized_fields = ["reference_id", *[field for field in metadata_fields if field != id_field]]
    _write_tsv(normalized_metadata, reference_metadata_path, normalized_fields)
    _write_tsv(manifest_rows, output / "reference_build_manifest.tsv", REFERENCE_BUILD_FIELDS)
    with (output / "myoga_metadata.csv").open("w", newline="") as handle:
        myoga_fields = ["genome_id", *[field for field in metadata_fields if field != id_field]]
        writer = csv.DictWriter(handle, fieldnames=myoga_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(myoga_metadata)

    tree_paths = build_reference_phylogenies(
        database_dir,
        output / "phylogeny",
        loci,
        resolve_threads(threads),
        min_references=min_references_per_tree,
        mafft_bin=mafft_bin,
        raxml_ng_bin=raxml_ng_bin,
        raxml_model=raxml_model,
    )
    return {
        "outdir": output,
        "database": database_dir,
        "metadata": reference_metadata_path,
        "myoga_metadata": output / "myoga_metadata.csv",
        "manifest": output / "reference_build_manifest.tsv",
        **tree_paths,
    }
