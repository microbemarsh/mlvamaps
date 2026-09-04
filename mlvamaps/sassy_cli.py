from __future__ import annotations

import csv
import io
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


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
        self._cache: dict[tuple[bytes, bytes], tuple[int, list[Match]]] = {}

    def prime(self, patterns: Iterable[bytes], text: bytes, k: int) -> None:
        """Search many primers against one sequence in a single Sassy process.

        In-silico PCR revisits every primer at mismatch thresholds ``0..k``.
        Caching the maximum-threshold result is safe because Sassy reports each
        match's edit cost; :meth:`search` filters that result for each round.
        """
        unique_patterns = tuple(dict.fromkeys(pattern for pattern in patterns if pattern))
        missing = [
            pattern
            for pattern in unique_patterns
            if (pattern, text) not in self._cache or self._cache[(pattern, text)][0] < k
        ]
        if not missing or not text:
            return

        pattern_texts = [pattern.decode("ascii") for pattern in missing]
        sequence_text = text.decode("ascii")
        with tempfile.TemporaryDirectory(prefix="mlvamaps-sassy-") as directory:
            directory_path = Path(directory)
            patterns_path = directory_path / "patterns.fasta"
            texts_path = directory_path / "text.fasta"
            patterns_path.write_text(
                "".join(f">p{index}\n{pattern}\n" for index, pattern in enumerate(pattern_texts))
            )
            texts_path.write_text(f">text\n{sequence_text}\n")
            completed = subprocess.run(
                self._command(patterns_path, texts_path, k),
                capture_output=True,
                text=True,
                check=False,
            )
        if completed.returncode:
            detail = completed.stderr.strip() or completed.stdout.strip() or "no diagnostic output"
            raise RuntimeError(f"Sassy search failed ({completed.returncode}): {detail}")

        matches_by_pattern: dict[str, list[Match]] = {
            f"p{index}": [] for index in range(len(missing))
        }
        for row in csv.DictReader(io.StringIO(completed.stdout), delimiter="\t"):
            if not row or not row.get("start"):
                continue
            pattern_id = row.get("pat_id", "")
            # ``pattern`` is retained for compatibility with the previous
            # one-pattern adapter and lightweight downstream test doubles.
            if pattern_id == "pattern" and len(missing) == 1:
                pattern_id = "p0"
            if pattern_id not in matches_by_pattern:
                continue
            pattern_index = int(pattern_id[1:])
            matches_by_pattern[pattern_id].append(
                self._match_from_row(row, len(pattern_texts[pattern_index]))
            )
        for index, pattern in enumerate(missing):
            self._cache[(pattern, text)] = (k, matches_by_pattern[f"p{index}"])

    def _command(self, patterns_path: Path, texts_path: Path, k: int) -> list[str]:
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
        return command

    @staticmethod
    def _match_from_row(row: dict[str, str], pattern_length: int) -> Match:
        return Match(
            text_start=int(row["start"]),
            text_end=int(row["end"]),
            cost=int(row["cost"]),
            strand=row["strand"],
            cigar=row["cigar"],
            pattern_start=0,
            pattern_end=pattern_length,
        )

    def search(self, pattern: bytes, text: bytes, k: int) -> list[Match]:
        if not pattern or not text:
            return []
        cached = self._cache.get((pattern, text))
        if cached is None or cached[0] < k:
            self.prime([pattern], text, k)
            cached = self._cache[(pattern, text)]
        return [match for match in cached[1] if match.cost <= k]