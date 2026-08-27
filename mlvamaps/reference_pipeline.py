"""Prepare NCBI taxon cohorts and build isolated mlvamaps references."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import shutil
import subprocess
import time
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from .concurrency import DEFAULT_THREADS
from .reference_builder import (
    REFERENCE_BUILD_STATUS_BUILT,
    REFERENCE_BUILD_STATUS_NO_USABLE_LOCI,
    REFERENCE_BUILD_STATUS_PARTIAL,
    build_reference_database,
)


NCBI_METADATA_FIELDS = (
    "accession,organism-name,organism-tax-id,assminfo-level,"
    "assminfo-release-date,assminfo-biosample-accession,"
    "assminfo-biosample-collection-date,assminfo-biosample-geo-loc-name,"
    "assminfo-biosample-lat-lon,assminfo-biosample-host,"
    "assminfo-biosample-isolation-source"
)
PREPARED_METADATA_FIELDS = (
    "accession",
    "assembly_file",
    "biosample_accession",
    "taxid",
    "organism_name",
    "assembly_level",
    "release_date",
    "collection_date",
    "latitude",
    "longitude",
    "country",
    "location",
    "sample_type",
    "isolation_source",
    "host",
)
TAXON_REFERENCE_SUMMARY_FIELDS = [
    "taxid",
    "taxon",
    "genomes_examined",
    "loci_total",
    "loci_amplifiable",
    "loci_not_amplifiable",
    "percent_loci_amplifiable",
    "total_valid_amplicons",
    "trees_built",
    "status",
]
TAXON_LOCUS_AMPLIFIABILITY_FIELDS = [
    "taxid",
    "taxon",
    "locus_id",
    "genomes_examined",
    "genomes_with_valid_amplicon",
    "valid_amplicons",
    "percent_genomes_amplifiable",
    "amplifiable",
    "tree_status",
]
_ACCESSION_PATTERN = re.compile(r"(GC[AF]_\d+\.\d+)", re.IGNORECASE)
_SAFE_NAME_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\Z")
_MISSING_VALUES = {
    "",
    "na",
    "n/a",
    "none",
    "null",
    "missing",
    "not collected",
    "not determined",
    "not applicable",
    "not available",
    "not provided",
    "unknown",
}


@dataclass(frozen=True)
class TaxonReference:
    taxid: str
    name: str


def _clean_taxid(value: object) -> str:
    taxid = str(value or "").strip()
    if not taxid.isascii() or not taxid.isdigit() or int(taxid) < 1:
        raise ValueError(f"invalid NCBI taxonomy identifier: {taxid!r}")
    return str(int(taxid))


def _clean_name(value: object, taxid: str) -> str:
    name = str(value or "").strip() or f"taxid_{taxid}"
    if not _SAFE_NAME_PATTERN.fullmatch(name) or name in {".", ".."}:
        raise ValueError(
            f"invalid reference name {name!r}; use only letters, numbers, dots, "
            "underscores, and hyphens"
        )
    return name


def read_taxon_references(
    *, taxid: str | int | None = None, taxids_csv: str | Path | None = None
) -> list[TaxonReference]:
    """Read one taxid or a CSV containing ``taxid`` and optional ``name`` columns."""
    if bool(taxid) == bool(taxids_csv):
        raise ValueError("provide exactly one of taxid or taxids_csv")
    if taxid is not None:
        clean_taxid = _clean_taxid(taxid)
        return [TaxonReference(clean_taxid, _clean_name("", clean_taxid))]

    csv_path = Path(str(taxids_csv))
    if not csv_path.is_file():
        raise ValueError(f"taxid CSV does not exist: {csv_path}")
    with csv_path.open(encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError(f"taxid CSV has no header: {csv_path}")
        fields = {str(field).strip().lower(): field for field in reader.fieldnames}
        taxid_field = next(
            (fields[name] for name in ("taxid", "taxon_id", "ncbi_taxid") if name in fields),
            None,
        )
        if taxid_field is None:
            raise ValueError(
                f"taxid CSV needs a taxid column (accepted: taxid, taxon_id, "
                f"ncbi_taxid): {csv_path}"
            )
        name_field = next(
            (fields[name] for name in ("name", "reference_name", "database_name") if name in fields),
            None,
        )
        references = []
        for line_number, row in enumerate(reader, start=2):
            raw_taxid = str(row.get(taxid_field) or "").strip()
            if not raw_taxid and not any(str(value or "").strip() for value in row.values()):
                continue
            try:
                clean_taxid = _clean_taxid(raw_taxid)
                clean_name = _clean_name(row.get(name_field) if name_field else "", clean_taxid)
            except ValueError as exc:
                raise ValueError(f"{csv_path}:{line_number}: {exc}") from exc
            references.append(TaxonReference(clean_taxid, clean_name))
    if not references:
        raise ValueError(f"taxid CSV contains no taxids: {csv_path}")
    duplicate_taxids = sorted(
        {item.taxid for item in references if sum(x.taxid == item.taxid for x in references) > 1}
    )
    duplicate_names = sorted(
        {item.name for item in references if sum(x.name == item.name for x in references) > 1}
    )
    if duplicate_taxids:
        raise ValueError(f"taxid CSV contains duplicate taxids: {', '.join(duplicate_taxids)}")
    if duplicate_names:
        raise ValueError(
            f"taxid CSV contains duplicate reference names: {', '.join(duplicate_names)}"
        )
    return references


def _resolve_executable(value: str, label: str) -> str:
    resolved = shutil.which(value)
    if resolved is None:
        raise RuntimeError(
            f"NCBI Datasets CLI executable `{value}` ({label}) is not available on PATH"
        )
    return resolved


def _run_command(command: list[str]) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if result.returncode:
        detail = (result.stderr or result.stdout).strip()
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}"
            + (f"\n{detail[-3000:]}" if detail else "")
        )
    return result


def _download_package(
    command: list[str], package: Path, attempts: int
) -> None:
    """Run an NCBI download with bounded retries and reject partial archives."""
    if attempts < 1:
        raise ValueError("download attempts must be at least 1")
    errors: list[str] = []
    for attempt in range(1, attempts + 1):
        package.unlink(missing_ok=True)
        try:
            _run_command(command)
            if not package.is_file():
                raise RuntimeError(
                    f"NCBI Datasets did not create the requested package: {package}"
                )
            if not zipfile.is_zipfile(package):
                raise RuntimeError(
                    f"NCBI Datasets created an incomplete or invalid ZIP archive: {package}"
                )
            return
        except RuntimeError as exc:
            errors.append(str(exc))
            package.unlink(missing_ok=True)
            if attempt < attempts:
                time.sleep(min(2 ** (attempt - 1), 10))
    raise RuntimeError(
        f"NCBI Datasets download failed after {attempts} attempt(s). "
        f"Last error:\n{errors[-1]}"
    )


def _tool_version(executable: str) -> str:
    for argument in ("--version", "version"):
        result = subprocess.run(
            [executable, argument], text=True, capture_output=True, check=False
        )
        output = (result.stdout or result.stderr).strip()
        if result.returncode == 0 and output:
            return output.splitlines()[0]
    return "unknown"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_extract(package: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with zipfile.ZipFile(package) as archive:
        for member in archive.infolist():
            target = (destination / member.filename).resolve()
            if target != root and root not in target.parents:
                raise RuntimeError(
                    f"NCBI package contains an unsafe archive path: {member.filename!r}"
                )
        archive.extractall(destination)


def _nonmissing(row: dict[str, Any], *names: str) -> str:
    for name in names:
        value = str(row.get(name) or row.get(f"#{name}") or "").strip()
        if value.casefold() not in _MISSING_VALUES:
            return value
    return ""


def _coordinates(value: str) -> tuple[str, str]:
    parts = value.replace(",", " ").split()
    try:
        if len(parts) == 2:
            latitude, longitude = float(parts[0]), float(parts[1])
        elif len(parts) >= 4:
            latitude = float(parts[0]) * (-1 if parts[1].upper() == "S" else 1)
            longitude = float(parts[2]) * (-1 if parts[3].upper() == "W" else 1)
        else:
            return "", ""
    except ValueError:
        return "", ""
    if not -90 <= latitude <= 90 or not -180 <= longitude <= 180:
        return "", ""
    return str(latitude), str(longitude)


def _assembly_accession(path: Path) -> str:
    match = _ACCESSION_PATTERN.search(path.name)
    if match is None:
        match = next(
            (
                _ACCESSION_PATTERN.search(parent.name)
                for parent in path.parents
                if _ACCESSION_PATTERN.search(parent.name)
            ),
            None,
        )
    return match.group(1).upper() if match else ""


def _assembly_paths(package_dir: Path) -> dict[str, Path]:
    paths = sorted(
        path
        for path in package_dir.rglob("*")
        if path.is_file()
        and path.name.lower().endswith((".fna", ".fna.gz", ".fasta", ".fasta.gz"))
        and "_cds_from_genomic." not in path.name.lower()
        and "_rna_from_genomic." not in path.name.lower()
    )
    by_accession: dict[str, Path] = {}
    for path in paths:
        accession = _assembly_accession(path)
        if not accession:
            continue
        if accession in by_accession:
            raise RuntimeError(f"NCBI package contains multiple assemblies for {accession}")
        by_accession[accession] = path
    if not by_accession:
        raise RuntimeError("NCBI package contained no recognizable genome assemblies")
    return by_accession


def _canonical_metadata(raw_tsv: str, assemblies: dict[str, Path], root: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    reader = csv.DictReader(raw_tsv.splitlines(), delimiter="\t")
    for raw in reader:
        accession = _nonmissing(raw, "Assembly Accession", "accession").upper()
        if not accession:
            continue
        assembly = assemblies.get(accession)
        if assembly is None:
            raise RuntimeError(f"metadata has no matching downloaded assembly for {accession}")
        if accession in seen:
            raise RuntimeError(f"NCBI metadata contains duplicate accession {accession}")
        seen.add(accession)
        location = _nonmissing(
            raw, "Assembly BioSample Geographic Location", "geo_loc_name"
        )
        isolation_source = _nonmissing(
            raw, "Assembly BioSample Isolation Source", "isolation_source"
        )
        latitude, longitude = _coordinates(
            _nonmissing(raw, "Assembly BioSample Latitude Longitude", "lat_lon")
        )
        rows.append(
            {
                "accession": accession,
                "assembly_file": str(assembly.relative_to(root)),
                "biosample_accession": _nonmissing(
                    raw, "Assembly BioSample Accession", "biosample_accession"
                ),
                "taxid": _nonmissing(raw, "Organism Taxonomic ID", "taxid"),
                "organism_name": _nonmissing(raw, "Organism Name", "organism_name"),
                "assembly_level": _nonmissing(raw, "Assembly Level", "assembly_level"),
                "release_date": _nonmissing(
                    raw, "Assembly Release Date", "release_date"
                ),
                "collection_date": _nonmissing(
                    raw, "Assembly BioSample Collection Date", "collection_date"
                ).split("/", 1)[0],
                "latitude": latitude,
                "longitude": longitude,
                "country": location.split(":", 1)[0].strip() if location else "",
                "location": location,
                "sample_type": isolation_source,
                "isolation_source": isolation_source,
                "host": _nonmissing(raw, "Assembly BioSample Host", "host"),
            }
        )
    missing_metadata = sorted(set(assemblies) - seen)
    if missing_metadata:
        raise RuntimeError(
            "downloaded assemblies are missing from NCBI metadata: "
            + ", ".join(missing_metadata[:10])
        )
    if not rows:
        raise RuntimeError("NCBI returned no assembly metadata")
    return sorted(rows, key=lambda row: row["accession"])


def prepare_taxon_reference(
    taxon: TaxonReference,
    output: str | Path,
    *,
    assembly_source: str = "refseq",
    datasets_args: list[str] | None = None,
    datasets_bin: str = "datasets",
    dataformat_bin: str = "dataformat",
    resume: bool = False,
    download_retries: int = 3,
) -> dict[str, Any]:
    """Download and normalize one NCBI taxid as a portable build input."""
    if assembly_source not in {"refseq", "genbank", "all"}:
        raise ValueError("assembly_source must be refseq, genbank, or all")
    datasets = _resolve_executable(datasets_bin, "datasets")
    dataformat = _resolve_executable(dataformat_bin, "dataformat")
    output_path = Path(output).resolve()
    package = output_path / "ncbi_dataset.zip"
    if output_path.exists() and not resume:
        existing = list(output_path.iterdir())
        interrupted_download = (
            existing == [package] and not zipfile.is_zipfile(package)
        )
        if interrupted_download:
            package.unlink()
        elif existing:
            raise ValueError(f"prepared output directory is not empty: {output_path}")
    output_path.mkdir(parents=True, exist_ok=True)
    download_command = [
        datasets,
        "download",
        "genome",
        "taxon",
        taxon.taxid,
    ]
    if assembly_source != "all":
        download_command += ["--assembly-source", assembly_source]
    download_command += [
        "--include",
        "genome",
        "--filename",
        str(package),
        *(datasets_args or []),
    ]
    reused_package = resume and package.is_file() and zipfile.is_zipfile(package)
    if resume and package.is_file() and not reused_package:
        package.unlink()
    if resume and not reused_package:
        if any(output_path.iterdir()):
            raise ValueError(
                f"cannot resume preparation because a valid package is missing: {package}"
            )
        _download_package(download_command, package, download_retries)
        reused_package = False
    elif not reused_package:
        _download_package(download_command, package, download_retries)

    package_dir = output_path / "package"
    _safe_extract(package, package_dir)
    assemblies = _assembly_paths(package_dir)
    metadata_command = [
        dataformat,
        "tsv",
        "genome",
        "--package",
        str(package),
        "--fields",
        NCBI_METADATA_FIELDS,
    ]
    rows = _canonical_metadata(
        _run_command(metadata_command).stdout, assemblies, package_dir
    )
    metadata_path = output_path / "metadata.tsv"
    with metadata_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=PREPARED_METADATA_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    manifest = {
        "schema_version": "1.0",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "reference_name": taxon.name,
        "selection": {
            "taxid": taxon.taxid,
            "assembly_source": assembly_source,
        },
        "commands": [download_command, metadata_command],
        "tool_versions": {
            "datasets": _tool_version(datasets),
            "dataformat": _tool_version(dataformat),
        },
        "download_reused": reused_package,
        "package_sha256": _sha256(package),
        "assembly_count": len(assemblies),
        "metadata_count": len(rows),
        "genomes": "package",
        "metadata": "metadata.tsv",
    }
    manifest_path = output_path / "download_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return {
        "outdir": output_path,
        "genomes": package_dir,
        "metadata": metadata_path,
        "manifest": manifest_path,
        "taxid": taxon.taxid,
        "name": taxon.name,
    }


def prepare_taxon_references(
    references: list[TaxonReference],
    outdir: str | Path,
    **kwargs: Any,
) -> list[dict[str, Any]]:
    """Prepare independent NCBI input packages for multiple taxids."""
    output = Path(outdir)
    return [
        prepare_taxon_reference(reference, output / reference.name / "prepared", **kwargs)
        for reference in references
    ]


def _write_tsv(rows: list[dict[str, Any]], path: Path, fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _print_taxon_summary(summary: dict[str, Any], non_amplifiable: list[str]) -> None:
    print(f"{summary['taxon']} [taxid {summary['taxid']}]")
    print(f"  Genomes examined: {summary['genomes_examined']:,}")
    print(
        f"  Amplifiable loci: {summary['loci_amplifiable']} / {summary['loci_total']} "
        f"({summary['percent_loci_amplifiable']:.1f}%)"
    )
    print(f"  Valid amplicons: {summary['total_valid_amplicons']:,}")
    print(f"  Trees built: {summary['trees_built']:,}")
    print(f"  Status: {summary['status']}")
    if non_amplifiable:
        print(f"  Non-amplifiable loci: {', '.join(non_amplifiable)}")


def build_taxon_references(
    references: list[TaxonReference],
    primers_path: str | Path,
    outdir: str | Path,
    *,
    loci_path: str | Path | None = None,
    assembly_source: str = "refseq",
    datasets_args: list[str] | None = None,
    datasets_bin: str = "datasets",
    dataformat_bin: str = "dataformat",
    resume: bool = False,
    download_retries: int = 3,
    multiple_products: str = "exclude",
    max_primer_mismatches: int = 2,
    min_references_per_tree: int = 3,
    threads: int = DEFAULT_THREADS,
    amplirust_bin: str = "amplirust",
    mafft_bin: str = "mafft",
    raxml_ng_bin: str = "raxml-ng",
    raxml_model: str = "DNA",
    show_progress: bool = False,
    builder: Callable[..., dict[str, Any]] = build_reference_database,
) -> dict[str, Any]:
    """Prepare and build one isolated mlvamaps database per taxid."""
    output = Path(outdir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    results = []
    taxon_summary_rows: list[dict[str, Any]] = []
    locus_summary_rows: list[dict[str, Any]] = []
    for reference in references:
        prepared = prepare_taxon_reference(
            reference,
            output / reference.name / "prepared",
            assembly_source=assembly_source,
            datasets_args=datasets_args,
            datasets_bin=datasets_bin,
            dataformat_bin=dataformat_bin,
            resume=resume,
            download_retries=download_retries,
        )
        database_outdir = output / reference.name / "reference"
        built = builder(
            assemblies_dir=prepared["genomes"],
            primers_path=primers_path,
            loci_path=loci_path,
            metadata_path=prepared["metadata"],
            outdir=database_outdir,
            multiple_products=multiple_products,
            max_primer_mismatches=max_primer_mismatches,
            min_references_per_tree=min_references_per_tree,
            threads=threads,
            amplirust_bin=amplirust_bin,
            mafft_bin=mafft_bin,
            raxml_ng_bin=raxml_ng_bin,
            raxml_model=raxml_model,
            show_progress=show_progress,
        )
        taxon_locus_rows = [
            {"taxid": reference.taxid, "taxon": reference.name, **row}
            for row in built.get("locus_summary_rows", [])
        ]
        locus_summary_rows.extend(taxon_locus_rows)
        loci_total = len(taxon_locus_rows)
        loci_amplifiable = sum(row["amplifiable"] == "TRUE" for row in taxon_locus_rows)
        trees_built = sum(row["tree_status"] == REFERENCE_BUILD_STATUS_BUILT for row in taxon_locus_rows)
        if loci_amplifiable == 0:
            status = REFERENCE_BUILD_STATUS_NO_USABLE_LOCI
        elif loci_amplifiable == loci_total and trees_built == loci_total:
            status = REFERENCE_BUILD_STATUS_BUILT
        else:
            status = REFERENCE_BUILD_STATUS_PARTIAL
        # Preserve compatibility with injected/custom builders without summary data.
        if not taxon_locus_rows:
            status = str(built.get("status", REFERENCE_BUILD_STATUS_BUILT))
        taxon_summary = {
            "taxid": reference.taxid,
            "taxon": reference.name,
            "genomes_examined": (
                taxon_locus_rows[0]["genomes_examined"] if taxon_locus_rows else 0
            ),
            "loci_total": loci_total,
            "loci_amplifiable": loci_amplifiable,
            "loci_not_amplifiable": loci_total - loci_amplifiable,
            "percent_loci_amplifiable": (
                round(100.0 * loci_amplifiable / loci_total, 1) if loci_total else 0.0
            ),
            "total_valid_amplicons": sum(
                int(row["valid_amplicons"]) for row in taxon_locus_rows
            ),
            "trees_built": trees_built,
            "status": status,
        }
        taxon_summary_rows.append(taxon_summary)
        if show_progress and taxon_locus_rows:
            _print_taxon_summary(
                taxon_summary,
                [row["locus_id"] for row in taxon_locus_rows if row["amplifiable"] == "FALSE"],
            )
        results.append(
            {
                "taxid": reference.taxid,
                "name": reference.name,
                "status": status,
                "prepared": str(prepared["outdir"]),
                "reference": str(built["outdir"]),
                "database": str(built["database"]),
                "manifest": str(built["manifest"]),
                "locus_amplifiability": str(built.get("locus_amplifiability", "")),
            }
        )
    taxon_summary_path = output / "taxon_reference_summary.tsv"
    locus_summary_path = output / "taxon_locus_amplifiability.tsv"
    _write_tsv(taxon_summary_rows, taxon_summary_path, TAXON_REFERENCE_SUMMARY_FIELDS)
    _write_tsv(locus_summary_rows, locus_summary_path, TAXON_LOCUS_AMPLIFIABILITY_FIELDS)
    pipeline_manifest = output / "reference_pipeline_manifest.json"
    pipeline_manifest.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "created_at": datetime.now(timezone.utc).isoformat(),
                "taxon_summary": taxon_summary_path.name,
                "locus_amplifiability": locus_summary_path.name,
                "references": results,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "outdir": output,
        "manifest": pipeline_manifest,
        "taxon_summary": taxon_summary_path,
        "locus_amplifiability": locus_summary_path,
        "references": results,
    }
