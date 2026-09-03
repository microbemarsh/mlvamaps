from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

from .concurrency import resolve_threads
from .progress import ProgressReporter


def check_deacon(executable: str = "deacon") -> str:
    path = shutil.which(executable)
    if path is None:
        raise RuntimeError(
            f"Deacon executable {executable!r} was not found. Install Deacon "
            "from Bioconda or pass --deacon-bin."
        )
    return path


def deacon_version(executable: str = "deacon") -> str:
    executable_path = check_deacon(executable)
    result = subprocess.run(
        [executable_path, "--version"], capture_output=True, text=True, check=False
    )
    if result.returncode:
        raise RuntimeError(f"Could not determine Deacon version from {executable_path}")
    return (result.stdout or result.stderr).strip()


def build_deacon_index(
    fasta: str | Path,
    index: str | Path,
    threads: int,
    executable: str = "deacon",
    *,
    kmer_size: int = 31,
    window_size: int = 15,
) -> Path:
    """Build a broad target-group recruitment index from real genomes."""
    executable_path = check_deacon(executable)
    index = Path(index)
    index.parent.mkdir(parents=True, exist_ok=True)
    command = [
        executable_path, "index", "build", "-k", str(kmer_size), "-w",
        str(window_size), "--threads", str(resolve_threads(threads)), "--quiet",
        "--output", str(index), str(fasta),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode:
        detail = result.stderr.strip() or result.stdout.strip() or "no diagnostic output"
        raise RuntimeError(f"Deacon indexing failed ({result.returncode}): {detail}")
    if not index.is_file():
        raise RuntimeError("Deacon completed without writing its recruitment index")
    return index


def build_deacon_filter_command(
    index_path: str | Path,
    reads_path: str | Path,
    output_path: str | Path,
    summary_path: str | Path,
    threads: int,
    absolute_threshold: int = 2,
    relative_threshold: float = 0.01,
    executable: str = "deacon",
) -> list[str]:
    if absolute_threshold < 1:
        raise ValueError("taxon-screen absolute threshold must be at least 1")
    if not 0 <= relative_threshold <= 1:
        raise ValueError("taxon-screen relative threshold must be between 0 and 1")
    return [
        executable,
        "filter",
        "--abs-threshold",
        str(absolute_threshold),
        "--rel-threshold",
        str(relative_threshold),
        "--threads",
        str(resolve_threads(threads)),
        "--output",
        str(output_path),
        "--summary",
        str(summary_path),
        "--quiet",
        str(index_path),
        str(reads_path),
    ]


def run_taxon_screen(
    reads_path: str | Path,
    index_path: str | Path,
    outdir: str | Path,
    threads: int,
    absolute_threshold: int = 1,
    relative_threshold: float = 0,
    executable: str = "deacon",
    progress: ProgressReporter | None = None,
) -> tuple[Path, Path, dict]:
    """Retain reads matching a target-taxon Deacon minimizer index."""
    reads_path = Path(reads_path)
    index_path = Path(index_path)
    if not reads_path.is_file():
        raise FileNotFoundError(f"Taxon-screen input reads do not exist: {reads_path}")
    if not index_path.is_file():
        raise FileNotFoundError(f"Taxon-screen index does not exist: {index_path}")
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    output_path = outdir / "taxon_screened_reads.fastq.gz"
    summary_path = outdir / "taxon_screen_summary.json"
    executable_path = check_deacon(executable)
    command = build_deacon_filter_command(
        index_path,
        reads_path,
        output_path,
        summary_path,
        threads,
        absolute_threshold,
        relative_threshold,
        executable_path,
    )
    progress = progress or ProgressReporter(enabled=False)
    progress.step("Screening reads against the target-taxon Deacon index")
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        detail = result.stderr.strip() or result.stdout.strip() or "no diagnostic output"
        raise RuntimeError(
            f"Deacon target-taxon screen failed (exit {result.returncode}): {detail}"
        )
    if not output_path.is_file():
        raise RuntimeError("Deacon completed without writing the screened FASTQ")
    if not summary_path.is_file():
        raise RuntimeError("Deacon completed without writing its JSON summary")
    try:
        summary = json.loads(summary_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Could not read Deacon summary {summary_path}: {exc}") from exc
    progress.step(
        "Taxon screen retained "
        f"{int(summary.get('seqs_out', 0)):,}/{int(summary.get('seqs_in', 0)):,} reads"
    )
    return output_path, summary_path, summary
