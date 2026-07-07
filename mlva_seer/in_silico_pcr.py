from __future__ import annotations

import csv
import shutil
import subprocess
from pathlib import Path

from .models import Locus
from .primers import read_loci_or_primers


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
    threads: int = 0,
    circular: bool = False,
    search_rc: bool = True,
    trim_primers: bool = False,
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
    threads: int = 0,
    circular: bool = False,
    search_rc: bool = True,
    trim_primers: bool = False,
    executable: str = "amplirust",
) -> dict[str, Path]:
    if shutil.which(executable) is None:
        raise RuntimeError(
            f"Could not find {executable!r} on PATH. Install amplirust first, for example with "
            "`conda install bioconda::amplirust`, then rerun this command."
        )

    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    loci = read_loci_or_primers(loci_path, primers_path)
    primers_path = write_amplirust_primers(loci, outdir_path / "amplirust_primers.csv")
    output_fasta = outdir_path / "amplirust_products.fasta"
    stats_tsv = outdir_path / "amplirust_stats.tsv"
    min_len, max_len = expected_amplicon_bounds(loci)
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
        executable=executable,
    )
    subprocess.run(command, check=True)
    return {"primers": primers_path, "products": output_fasta, "stats": stats_tsv}
