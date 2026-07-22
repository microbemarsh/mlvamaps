from __future__ import annotations

import csv
import re
import shutil
import subprocess
from pathlib import Path

from .concurrency import DEFAULT_THREADS
from .models import Locus
from .primers import read_loci_or_primers


AMPLIRUST_TSV_FIELDS = [
    "amplicon_id",
    "reference_id",
    "source_file",
    "primer_name",
    "product_len",
    "full_len",
    "fwd_start",
    "fwd_end",
    "fwd_mismatches",
    "fwd_identity",
    "fwd_cigar",
    "rev_start",
    "rev_end",
    "rev_mismatches",
    "rev_identity",
    "rev_cigar",
    "strand",
    "is_circular_wrap",
    "product_seq",
]


def write_amplirust_primers(loci: list[Locus], path: str | Path) -> Path:
    """Write amplirust primer-pair CSV from the MLVA loci table."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "forward", "reverse"])
        writer.writeheader()
        for locus in loci:
            writer.writerow(
                {
                    "name": locus.locus_id,
                    "forward": locus.forward_primer,
                    "reverse": locus.reverse_primer,
                }
            )
    return path


def expected_amplicon_bounds(loci: list[Locus]) -> tuple[int, int]:
    mins = [locus.expected_amplicon_min_bp for locus in loci if locus.expected_amplicon_min_bp > 0]
    maxes = [locus.expected_amplicon_max_bp for locus in loci if locus.expected_amplicon_max_bp > 0]
    return (min(mins) if mins else 50, max(maxes) if maxes else 5000)


def build_amplirust_command(
    input_path: str | Path,
    primers_path: str | Path,
    output_fasta: str | Path,
    stats_tsv: str | Path,
    min_len: int,
    max_len: int,
    max_errors: int = 2,
    threads: int = DEFAULT_THREADS,
    circular: bool = False,
    search_rc: bool = True,
    trim_primers: bool = False,
    max_n_fraction: float = 0.0,
    executable: str = "amplirust",
) -> list[str]:
    command = [
        executable,
        "--input",
        str(input_path),
        "--primers",
        str(primers_path),
        "--output",
        str(output_fasta),
        "--tsv",
        str(stats_tsv),
        "--max-errors",
        str(max_errors),
        "--min-len",
        str(min_len),
        "--max-len",
        str(max_len),
        "--threads",
        str(threads),
        "--max-n-fraction",
        str(max_n_fraction),
        "--quiet",
    ]
    if circular:
        command.append("--circular")
    if search_rc:
        command.append("--search-rc")
    if trim_primers:
        command.append("--trim-primers")
    return command


def run_amplirust(
    input_path: str,
    loci_path: str | None,
    outdir: str,
    primers_path: str | None = None,
    max_errors: int = 2,
    threads: int = DEFAULT_THREADS,
    circular: bool = False,
    search_rc: bool = True,
    trim_primers: bool = False,
    max_n_fraction: float = 0.0,
    executable: str = "amplirust",
) -> dict[str, Path]:
    loci = read_loci_or_primers(loci_path, primers_path)
    return run_amplirust_loci(
        input_path,
        loci,
        outdir,
        max_errors=max_errors,
        threads=threads,
        circular=circular,
        search_rc=search_rc,
        trim_primers=trim_primers,
        max_n_fraction=max_n_fraction,
        executable=executable,
    )


def run_amplirust_loci(
    input_path: str | Path,
    loci: list[Locus],
    outdir: str | Path,
    max_errors: int = 2,
    threads: int = DEFAULT_THREADS,
    circular: bool = False,
    search_rc: bool = True,
    trim_primers: bool = False,
    max_n_fraction: float = 0.0,
    executable: str = "amplirust",
    amplicon_bounds: tuple[int, int] | None = None,
) -> dict[str, Path]:
    """Run Amplirust for an already-loaded MLVA panel."""
    if shutil.which(executable) is None:
        raise RuntimeError(
            f"Could not find {executable!r} on PATH. Install amplirust first, for example with "
            "`conda install bioconda::amplirust`, then rerun this command."
        )

    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    primers_path = write_amplirust_primers(loci, outdir_path / "amplirust_primers.csv")
    output_fasta = outdir_path / "amplirust_products.fasta"
    stats_tsv = outdir_path / "amplirust_stats.tsv"
    # Amplirust creates output files lazily only after finding its first product.
    # Seed empty outputs so zero-hit runs and reruns cannot expose stale results.
    output_fasta.write_text("")
    stats_tsv.write_text("\t".join(AMPLIRUST_TSV_FIELDS) + "\n")
    min_len, max_len = amplicon_bounds or expected_amplicon_bounds(loci)
    command = build_amplirust_command(
        input_path=input_path,
        primers_path=primers_path,
        output_fasta=output_fasta,
        stats_tsv=stats_tsv,
        min_len=min_len,
        max_len=max_len,
        max_errors=max_errors,
        threads=threads,
        circular=circular,
        search_rc=search_rc,
        trim_primers=trim_primers,
        max_n_fraction=max_n_fraction,
        executable=executable,
    )
    subprocess.run(command, check=True)
    return {"primers": primers_path, "products": output_fasta, "stats": stats_tsv}


_AMPLIRUST_POSITION_RE = re.compile(r"\bpos=(\d+)-(\d+)\b")


def read_amplirust_results(
    stats_path: str | Path, products_path: str | Path
) -> list[dict[str, str | int]]:
    """Read Amplirust TSV rows and attach original-reference product coordinates."""
    coordinates: dict[str, tuple[int, int]] = {}
    with Path(products_path).open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            fields = line[1:].rstrip().split("\t")
            match = _AMPLIRUST_POSITION_RE.search(line)
            if match:
                coordinates[fields[0]] = (int(match.group(1)), int(match.group(2)))

    rows: list[dict[str, str | int]] = []
    with Path(stats_path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            amplicon_id = row.get("amplicon_id", "")
            if amplicon_id not in coordinates:
                raise RuntimeError(
                    f"Amplirust TSV product {amplicon_id!r} has no matching FASTA record"
                )
            original_start, original_end = coordinates[amplicon_id]
            rows.append(
                {
                    **row,
                    "original_start": original_start,
                    "original_end": original_end,
                }
            )
    return rows
