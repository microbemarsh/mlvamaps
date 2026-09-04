"""Shared competitive minimap2 mapping for all FASTQ technologies."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pysam

from .alignment_evidence import CandidateAlignment
from .candidate_contexts import CandidateContext
from .io import normalize_read_id


PRESETS: dict[str, tuple[str, ...]] = {
    # ``sr`` disables secondary output in its preset; explicitly turn it back on
    # below because adjacent allele contexts are expected competitors.
    "illumina": ("-x", "sr"),
    # MLVA contexts are much shorter than genomes. The trailing seed/chaining
    # overrides retain the technology scoring model while permitting short
    # primer-bounded molecules to align competitively.
    "ont": ("-x", "map-ont", "-k", "9", "-w", "5", "-m", "10", "-s", "10", "-n", "2"),
    "ont-hq": ("-x", "lr:hq", "-k", "11", "-w", "5", "-m", "10", "-s", "10", "-n", "2"),
    "hifi": ("-x", "map-hifi", "-k", "11", "-w", "5", "-m", "10", "-s", "10", "-n", "2"),
}

# Index parameters are deliberately technology-specific even though both
# indexes contain exactly the same biological candidate sequences.
INDEX_PARAMETERS: dict[str, tuple[int, int]] = {
    "short": (21, 11),
    "long": (11, 5),
}


def minimap2_version(executable: str) -> str:
    result = subprocess.run([executable, "--version"], capture_output=True, text=True, check=False)
    if result.returncode:
        raise RuntimeError(f"Could not run minimap2 at {executable}")
    return result.stdout.strip()


def build_minimap2_index(
    fasta: str | Path,
    index: str | Path,
    executable: str = "minimap2",
    *,
    kmer_size: int | None = None,
    window_size: int | None = None,
) -> Path:
    resolved = shutil.which(executable) or (executable if Path(executable).is_file() else None)
    if resolved is None:
        raise RuntimeError(f"minimap2 executable {executable!r} was not found")
    index = Path(index)
    index.parent.mkdir(parents=True, exist_ok=True)
    command = [resolved]
    if kmer_size is not None:
        command.extend(["-k", str(kmer_size)])
    if window_size is not None:
        command.extend(["-w", str(window_size)])
    command.extend(["-d", str(index), str(fasta)])
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode:
        raise RuntimeError(f"minimap2 indexing failed: {result.stderr.strip()}")
    return index


def build_competitive_indexes(
    fasta: str | Path,
    directory: str | Path,
    executable: str = "minimap2",
) -> dict[str, Path]:
    """Build reusable short- and long-read indexes over one candidate FASTA."""
    directory = Path(directory)
    paths: dict[str, Path] = {}
    for technology, (kmer_size, window_size) in INDEX_PARAMETERS.items():
        paths[technology] = build_minimap2_index(
            fasta,
            directory / f"{technology}.mmi",
            executable,
            kmer_size=kmer_size,
            window_size=window_size,
        )
    return paths


def minimap2_competitive_command(
    reference: str | Path,
    reads1: str | Path,
    reads2: str | Path | None,
    threads: int,
    technology: str,
    executable: str = "minimap2",
    max_secondary: int = 100,
) -> list[str]:
    if technology not in PRESETS:
        raise ValueError(f"Unknown minimap2 technology preset: {technology!r}")
    command = [executable, "-a"]
    command.extend(PRESETS[technology])
    command.extend(["--cs=long", "--secondary=yes", "-N", str(max_secondary)])
    command.extend(["-t", str(threads), str(reference), str(reads1)])
    if reads2 is not None:
        command.append(str(reads2))
    return command


def run_minimap2_competitive(command: list[str], sam_path: str | Path) -> None:
    with Path(sam_path).open("w") as output:
        result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE, text=True)
    if result.returncode:
        raise RuntimeError(
            f"minimap2 competitive candidate mapping failed ({result.returncode}): "
            f"{result.stderr.strip()}"
        )


def run_minimap2_competitive_bam(command: list[str], bam_path: str | Path) -> None:
    """Stream minimap2 SAM directly through htslib into compressed BAM."""
    bam_path = Path(bam_path)
    bam_path.parent.mkdir(parents=True, exist_ok=True)
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
        assert process.stdout is not None
        try:
            with pysam.AlignmentFile(process.stdout, "r") as source, pysam.AlignmentFile(
                str(bam_path), "wb", template=source
            ) as output:
                for alignment in source.fetch(until_eof=True):
                    output.write(alignment)
        except Exception:
            process.kill()
            process.wait()
            bam_path.unlink(missing_ok=True)
            raise
        stderr = (process.stderr.read() if process.stderr is not None else b"").decode(
            errors="replace"
        )
        returncode = process.wait()
    if returncode:
        bam_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"minimap2 competitive candidate mapping failed ({returncode}): {stderr.strip()}"
        )


def parse_candidate_alignments(
    sam_path: str | Path,
    contexts: list[CandidateContext],
) -> list[CandidateAlignment]:
    context_by_id = {context.candidate_id: context for context in contexts}
    rows: list[CandidateAlignment] = []
    with pysam.AlignmentFile(str(sam_path), "r") as sam:
        for alignment in sam.fetch(until_eof=True):
            if alignment.is_unmapped or alignment.is_supplementary:
                continue
            reference_name = sam.get_reference_name(alignment.reference_id)
            context = context_by_id.get(reference_name)
            if context is None or context.repeat_count is None:
                continue
            molecule_id, named_mate = normalize_read_id(alignment.query_name)
            mate = 1 if alignment.is_read1 else 2 if alignment.is_read2 else named_mate
            aligned = max(alignment.query_alignment_length, 0)
            nm = int(alignment.get_tag("NM")) if alignment.has_tag("NM") else 0
            denominator = max(aligned + sum(
                length for operation, length in (alignment.cigartuples or []) if operation == 2
            ), 1)
            identity = max(0.0, 1.0 - nm / denominator)
            rows.append(CandidateAlignment(
                molecule_id=molecule_id,
                read_id=alignment.query_name,
                mate=mate,
                locus_id=context.locus_id,
                candidate_id=context.candidate_id,
                repeat_count=context.repeat_count,
                reference_id=context.reference_accession or context.reference_id,
                alignment_score=float(alignment.get_tag("AS")) if alignment.has_tag("AS") else float(aligned),
                mapping_quality=int(alignment.mapping_quality),
                alignment_identity=identity,
                query_coverage=aligned / max(alignment.query_length or aligned, 1),
                reference_coverage=(alignment.reference_end - alignment.reference_start) / max(len(context.sequence), 1),
                query_start=alignment.query_alignment_start,
                query_end=alignment.query_alignment_end,
                reference_start=alignment.reference_start,
                reference_end=alignment.reference_end,
                cigar=alignment.cigarstring or "",
                cs=str(alignment.get_tag("cs")) if alignment.has_tag("cs") else "",
                primary=not alignment.is_secondary,
                secondary=alignment.is_secondary,
                reverse=alignment.is_reverse,
                query_sequence=alignment.get_forward_sequence() or "",
                query_quality=(
                    "".join(chr(value + 33) for value in alignment.get_forward_qualities())
                    if alignment.query_qualities is not None else None
                ),
                template_length=alignment.template_length,
                next_reference_start=alignment.next_reference_start,
            ))
    return rows


def map_reads_to_candidates(
    reference: str | Path,
    reads1: str | Path,
    reads2: str | Path | None,
    contexts: list[CandidateContext],
    sam_path: str | Path,
    threads: int,
    technology: str,
    executable: str = "minimap2",
    max_secondary: int = 100,
) -> list[CandidateAlignment]:
    resolved = shutil.which(executable) or (executable if Path(executable).is_file() else None)
    if resolved is None:
        raise RuntimeError(f"minimap2 executable {executable!r} was not found")
    command = minimap2_competitive_command(
        reference, reads1, reads2, threads, technology, resolved, max_secondary
    )
    run_minimap2_competitive(command, sam_path)
    return parse_candidate_alignments(sam_path, contexts)


def map_reads_to_candidates_bam(
    reference: str | Path,
    reads1: str | Path,
    reads2: str | Path | None,
    contexts: list[CandidateContext],
    bam_path: str | Path,
    threads: int,
    technology: str,
    executable: str = "minimap2",
    max_secondary: int = 100,
) -> list[CandidateAlignment]:
    """Map without materializing text SAM; pysam/htslib decodes CIGAR and tags."""
    resolved = shutil.which(executable) or (executable if Path(executable).is_file() else None)
    if resolved is None:
        raise RuntimeError(f"minimap2 executable {executable!r} was not found")
    command = minimap2_competitive_command(
        reference, reads1, reads2, threads, technology, resolved, max_secondary
    )
    run_minimap2_competitive_bam(command, bam_path)
    return parse_candidate_alignments(bam_path, contexts)