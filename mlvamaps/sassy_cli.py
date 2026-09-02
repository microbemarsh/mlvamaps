from __future__ import annotations

import csv
import io
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Match:
    """Subset of a Sassy match consumed by mlvamaps."""

    text_start: int
    text_end: int
    cost: int
    strand: str
    cigar: str
    pattern_start: int
    pattern_end: int


def _resolve_executable() -> str:
    configured = os.environ.get("SASSY_BIN", "sassy")
    resolved = shutil.which(configured)
    if resolved is None:
        raise RuntimeError(
            "Sassy is required for approximate sequence matching. Install the "
            "Bioconda package with `conda install -c bioconda sassy`, or set "
            "SASSY_BIN to the executable path."
        )
    return resolved


class Searcher:
    """Compatibility adapter for Bioconda's ``sassy search`` executable.

    The CLI accepts FASTA/FASTQ rather than in-memory sequences. Each call uses
    private temporary FASTA files, making the adapter safe across processes and
    threads. Reverse-complement searching is disabled unless explicitly asked
    for because mlvamaps applies its own strand semantics.
    """

    def __init__(
        self,
        alphabet: str = "iupac",
        rc: bool = True,
        max_n_frac: float | None = None,
    ) -> None:
        if alphabet not in {"ascii", "dna", "iupac"}:
            raise ValueError(f"Unsupported Sassy alphabet: {alphabet!r}")
        # mlvamaps' ASCII searches contain concrete DNA primers, but assembly
        # targets may contain N. The CLI does not expose the library's ASCII
        # mode; IUPAC accepts those targets and max_n_frac=0 rejects matches
        # overlapping N, preserving the binding's effective PCR semantics.
        self.alphabet = "iupac" if alphabet == "ascii" else alphabet
        self.rc = rc
        self.max_n_frac = 1.0 if max_n_frac is None else max_n_frac

    def search(self, pattern: bytes, text: bytes, k: int) -> list[Match]:
        pattern_text = pattern.decode("ascii")
        sequence_text = text.decode("ascii")
        if not pattern_text or not sequence_text:
            return []

        with tempfile.TemporaryDirectory(prefix="mlvamaps-sassy-") as directory:
            directory_path = Path(directory)
            patterns_path = directory_path / "pattern.fasta"
            texts_path = directory_path / "text.fasta"
            patterns_path.write_text(f">pattern\n{pattern_text}\n")
            texts_path.write_text(f">text\n{sequence_text}\n")
            command = [
                _resolve_executable(),
                "search",
                "--pattern-fasta",
                str(patterns_path),
                "-k",
                str(k),
                "--alphabet",
                self.alphabet,
                "--max-n-frac",
                str(self.max_n_frac),
                "--threads",
                "1",
            ]
            if not self.rc:
                command.append("--no-rc")
            command.append(str(texts_path))
            completed = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False,
            )
        if completed.returncode:
            detail = completed.stderr.strip() or completed.stdout.strip() or "no diagnostic output"
            raise RuntimeError(f"Sassy search failed ({completed.returncode}): {detail}")

        matches: list[Match] = []
        for row in csv.DictReader(io.StringIO(completed.stdout), delimiter="\t"):
            if not row or not row.get("start"):
                continue
            matches.append(
                Match(
                    text_start=int(row["start"]),
                    text_end=int(row["end"]),
                    cost=int(row["cost"]),
                    strand=row["strand"],
                    cigar=row["cigar"],
                    pattern_start=0,
                    pattern_end=len(pattern_text),
                )
            )
        return matches