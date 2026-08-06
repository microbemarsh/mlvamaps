from __future__ import annotations

import math
import shutil
import statistics
import subprocess
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

from .calling import normalize_allele, repeat_unit_length
from .concurrency import DEFAULT_THREADS, resolve_threads
from .io import (
    open_text,
    gzip_output_file,
    normalize_read_id,
    read_fasta,
    read_fastq_pairs,
    read_profiles,
    write_fasta,
    write_fastq,
    write_tsv,
)
from .mapping import check_minimap2
from .locus_measurement import find_anchor, measure_locus_product
from .models import Locus, LocusMeasurement, ReadPair, ReadRecord
from .profile_matching import (
    PROFILE_MATCH_LOCUS_FIELDS,
    build_fingerprint,
    match_profiles,
    profile_match_locus_rows,
)
from .pipeline import (
    ALLELE_DISTRIBUTION_FIELDS,
    MATCH_FIELDS,
    REPEAT_COUNT_FIELDS,
    SIMPLE_CALL_FIELDS,
    allele_distribution_rows,
)
from .primers import read_loci_or_primers
from .recruitment import build_recruitment_references
from .report import write_report
from .sample_metadata import MYOGA_SAMPLE_FIELDS, myoga_sample_row, write_csv
from .sequence import mean_qscore, revcomp


EVIDENCE_CLASSES = {
    "COMPLETE_ASSEMBLED_PRODUCT",
    "BOUNDARY_SPANNING_READ_PAIR",
    "BOUNDARY_SPANNING_SINGLE_READ",
    "PARTIAL_REPEAT_EVIDENCE",
    "PRESENCE_ONLY",
    "AMBIGUOUS_ASSEMBLY",
    "MULTIPLE_ALLELES",
    "LOW_DEPTH",
    "NOT_FOUND",
}

SHORT_CALL_EXTRA_FIELDS = [
    "read_technology",
    "evidence_class",
    "informative_molecule_count",
    "boundary_1_support",
    "boundary_2_support",
    "both_boundary_support",
    "local_assembly_status",
    "repeat_count_min",
    "repeat_count_max",
    "repeat_count_interval_reason",
    "confidence_reason",
    "short_read_warning",
    "recruited_read_pairs",
    "assembled_contig_length",
    "assembly_depth",
    "failure_reason",
    "primary_allele_support",
    "secondary_allele_support",
    "informative_allele_reads",
    "uninformative_locus_reads",
    "estimated_primary_fraction",
    "estimated_secondary_fraction",
    "mixture_status",
    "evidence_sources",
]
SHORT_CALL_FIELDS = SIMPLE_CALL_FIELDS + SHORT_CALL_EXTRA_FIELDS

SHORT_QC_FIELDS = ["sample_id", "metric", "value"]
SHORT_RECRUITMENT_FIELDS = [
    "sample_id",
    "locus_id",
    "read_pairs_examined",
    "read_pairs_recruited",
    "uniquely_recruited_pairs",
    "ambiguous_pairs",
    "discordant_pairs",
    "orphan_reads",
    "mean_mapping_quality",
    "mean_alignment_identity",
]
SHORT_ASSEMBLY_FIELDS = [
    "sample_id",
    "locus_id",
    "recruited_reads",
    "recruited_pairs",
    "assembler",
    "assembly_status",
    "contig_count",
    "longest_contig_bp",
    "total_contig_bp",
    "estimated_depth",
    "failure_reason",
]
SAMPLE_SUMMARY_FIELDS = [
    "sample_id",
    "input_read_1",
    "input_read_2",
    "read_technology",
    "sample_mode",
    "total_reads",
    "total_read_pairs",
    "retained_reads",
    "retained_pairs",
    "callable_loci",
    "complete_loci",
    "partial_loci",
    "presence_only_loci",
    "mixed_loci",
    "missing_loci",
    "best_profile_id",
    "best_profile_distance",
    "profile_confidence",
    "run_status",
    "warnings",
]


@dataclass(frozen=True)
class RecruitedPair:
    pair: ReadPair
    outcome: str
    locus_id: str = ""
    score: int = 0
    identity: float = 0.0


@dataclass(frozen=True)
class LocusAssembly:
    locus_id: str
    contigs: tuple[str, ...]
    merged_reads: tuple[ReadRecord, ...]
    estimated_depth: float
    status: str
    failure_reason: str = ""
    assembler: str = "skesa"


def _trim_read(read: ReadRecord, trim_quality: int) -> ReadRecord:
    if trim_quality <= 0 or read.quality is None:
        return read
    end = len(read.sequence)
    while end and ord(read.quality[end - 1]) - 33 < trim_quality:
        end -= 1
    return ReadRecord(read.read_id, read.sequence[:end], read.quality[:end])


def _read_passes(read: ReadRecord, min_length: int, min_mean_quality: float) -> bool:
    return len(read.sequence) >= min_length and mean_qscore(read.quality) >= min_mean_quality


def qc_read_pairs(
    pairs: list[ReadPair],
    min_length: int,
    min_mean_quality: float,
    trim_quality: int,
    min_pair_retention: float,
) -> tuple[list[ReadPair], dict[str, int]]:
    if min_length < 1 or min_mean_quality < 0 or trim_quality < 0:
        raise ValueError("short-read QC thresholds must be non-negative")
    if not 0 <= min_pair_retention <= 1:
        raise ValueError("short-min-pair-retention must be between 0 and 1")
    retained: list[ReadPair] = []
    metrics = Counter()
    for pair in pairs:
        metrics["input_pairs"] += 1
        mates = [pair.read1] + ([pair.read2] if pair.read2 is not None else [])
        metrics["input_reads"] += len(mates)
        trimmed = [_trim_read(read, trim_quality) for read in mates]
        passing = [
            read
            for read in trimmed
            if _read_passes(read, min_length, min_mean_quality)
        ]
        if len(passing) / len(mates) < min_pair_retention or not passing:
            metrics["rejected_pairs"] += 1
            metrics["rejected_reads"] += len(mates)
            continue
        read1 = passing[0]
        read2 = passing[1] if len(passing) == 2 else None
        if len(mates) == 2 and len(passing) == 1:
            metrics["orphan_reads"] += 1
            if passing[0].read_id == pair.read2.read_id:  # type: ignore[union-attr]
                read1 = passing[0]
        retained.append(ReadPair(pair.molecule_id, read1, read2))
        metrics["retained_pairs"] += 1
        metrics["retained_reads"] += len(passing)
    return retained, dict(metrics)


def _canonical_kmer(kmer: str) -> str:
    reverse = revcomp(kmer)
    return min(kmer, reverse)


def _kmers(sequence: str, k: int) -> set[str]:
    if len(sequence) < k:
        return set()
    return {
        _canonical_kmer(sequence[index : index + k])
        for index in range(len(sequence) - k + 1)
        if set(sequence[index : index + k]) <= set("ACGT")
    }


def _short_read_reference_sequences(
    loci: list[Locus], database_path: str | Path | None
) -> dict[str, list[str]]:
    references = build_recruitment_references(loci, database_path)
    by_locus: dict[str, list[str]] = defaultdict(list)
    for reference in references:
        by_locus[str(reference["locus_id"])].append(str(reference["sequence"]))
    return dict(by_locus)


def build_short_read_references(
    loci: list[Locus], database_path: str | Path | None
) -> dict[str, str]:
    by_locus = _short_read_reference_sequences(loci, database_path)
    representative = {
        locus_id: min(sequences, key=lambda value: (abs(len(value) - statistics.median(map(len, sequences))), value))
        for locus_id, sequences in by_locus.items()
    }
    # Primer-only panels do not have a synthesizable complete product. Keep
    # their two primer targets separated by Ns so native recruitment retains
    # the former locus-specific primer-seed behavior without inventing a
    # primer-junction sequence.
    for locus in loci:
        if locus.locus_id in representative:
            continue
        primer_targets = [
            sequence
            for sequence in (locus.forward_primer, revcomp(locus.reverse_primer))
            if sequence and set(sequence) <= set("ACGT")
        ]
        if primer_targets:
            representative[locus.locus_id] = ("N" * 32).join(primer_targets)
    return representative


def build_short_read_reference_index(
    loci: list[Locus], database_path: str | Path | None, k: int = 15
) -> tuple[dict[str, set[str]], dict[str, str]]:
    by_locus = _short_read_reference_sequences(loci, database_path)
    representative = {
        locus_id: min(sequences, key=lambda value: (abs(len(value) - statistics.median(map(len, sequences))), value))
        for locus_id, sequences in by_locus.items()
    }
    owners: dict[str, set[str]] = defaultdict(set)
    for locus_id, sequences in by_locus.items():
        for sequence in sequences:
            for kmer in _kmers(sequence, k):
                owners[kmer].add(locus_id)
    unique = {
        locus_id: {kmer for kmer, loci_for_kmer in owners.items() if loci_for_kmer == {locus_id}}
        for locus_id in by_locus
    }
    # Primer-only panels can lack A/C/G/T reference products when their legacy
    # repeat motif is unknown (N). Preserve unique primer evidence in that case.
    for locus in loci:
        primer_sequence = locus.forward_primer + revcomp(locus.reverse_primer)
        unique.setdefault(locus.locus_id, set()).update(
            kmer
            for kmer in _kmers(primer_sequence, min(k, max(5, min(len(locus.forward_primer), len(locus.reverse_primer), k))))
            if sum(kmer in _kmers(other.forward_primer + revcomp(other.reverse_primer), len(kmer)) for other in loci) == 1
        )
    return unique, representative


def _read_locus_scores(read: ReadRecord, index: dict[str, set[str]], k: int) -> dict[str, int]:
    observed = _kmers(read.sequence, k)
    return {locus_id: len(observed & seeds) for locus_id, seeds in index.items()}


def _best_ungapped_identity(sequence: str, reference: str) -> float:
    if not sequence or not reference:
        return 0.0
    if len(sequence) > len(reference):
        sequence, reference = reference, sequence
    best = 0
    for query in (sequence, revcomp(sequence)):
        for start in range(len(reference) - len(query) + 1):
            matches = sum(
                observed == expected
                for observed, expected in zip(query, reference[start : start + len(query)])
            )
            best = max(best, matches)
    return best / len(sequence)


def recruit_read_pairs(
    pairs: list[ReadPair],
    index: dict[str, set[str]],
    k: int = 15,
    min_seeds: int = 2,
    min_margin: int = 1,
    references: dict[str, str] | None = None,
) -> list[RecruitedPair]:
    recruited: list[RecruitedPair] = []
    for pair in pairs:
        mate_scores = [_read_locus_scores(pair.read1, index, k)]
        if pair.read2 is not None:
            mate_scores.append(_read_locus_scores(pair.read2, index, k))
        confident_mates: list[str] = []
        for scores in mate_scores:
            ranked = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
            if not ranked or ranked[0][1] < min_seeds:
                continue
            second = ranked[1][1] if len(ranked) > 1 else 0
            if ranked[0][1] - second >= min_margin:
                confident_mates.append(ranked[0][0])
        if len(set(confident_mates)) > 1:
            recruited.append(RecruitedPair(pair, "discordant"))
            continue
        combined = {
            locus_id: sum(scores.get(locus_id, 0) for scores in mate_scores)
            for locus_id in index
        }
        ranked = sorted(combined.items(), key=lambda item: (-item[1], item[0]))
        if confident_mates:
            locus_id = confident_mates[0]
            score = combined[locus_id]
        elif not ranked or ranked[0][1] < min_seeds:
            recruited.append(RecruitedPair(pair, "unassigned"))
            continue
        else:
            second = ranked[1][1] if len(ranked) > 1 else 0
            if ranked[0][1] - second < min_margin:
                recruited.append(RecruitedPair(pair, "ambiguous"))
                continue
            locus_id, score = ranked[0]
        if references and references.get(locus_id):
            identity = statistics.mean(
                _best_ungapped_identity(read.sequence, references[locus_id])
                for read in (pair.read1, pair.read2)
                if read is not None
            )
        else:
            possible = sum(max(0, len(read.sequence) - k + 1) for read in (pair.read1, pair.read2) if read is not None)
            identity = min(1.0, score / max(possible, 1))
        recruited.append(RecruitedPair(pair, "unique", locus_id, score, identity))
    return recruited


def _short_read_minimap2_command(
    executable: str,
    reference_path: Path,
    reads1_path: Path,
    reads2_path: Path | None,
    threads: int,
) -> list[str]:
    """Build the native competitive-recruitment command.

    PAF is intentional: unlike SAM, it emits only mapped queries by default,
    which avoids serializing every unassigned WGS read back through Python.
    """
    command = [
        executable,
        "-x",
        "sr",
        "-t",
        str(resolve_threads(threads)),
        "-k",
        "15",
        "-w",
        "5",
        "--secondary=yes",
        str(reference_path),
        str(reads1_path),
    ]
    if reads2_path is not None:
        command.append(str(reads2_path))
    return command


def _run_short_read_minimap2(command: list[str], output_path: Path) -> None:
    with output_path.open("w") as output:
        completed = subprocess.run(
            command,
            stdout=output,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    if completed.returncode:
        detail = (completed.stderr or "").strip()
        raise RuntimeError(
            f"minimap2 short-read recruitment failed with exit code "
            f"{completed.returncode}" + (f": {detail}" if detail else "")
        )


def _parse_short_read_paf(
    path: Path,
    reference_loci: dict[str, str],
) -> dict[str, dict[int, dict[str, tuple[int, float, int]]]]:
    """Collect the best native alignment per molecule, mate, and locus."""
    hits: dict[str, dict[int, dict[str, tuple[int, float, int]]]] = defaultdict(
        lambda: defaultdict(dict)
    )
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise RuntimeError(
                    f"Malformed minimap2 PAF at {path}:{line_number}: expected at least 12 fields"
                )
            reference_name = fields[5]
            locus_id = reference_loci.get(reference_name)
            if locus_id is None:
                raise RuntimeError(
                    f"minimap2 returned unknown recruitment reference {reference_name!r}"
                )
            molecule_id, mate = normalize_read_id(fields[0])
            # In paired mode minimap2 appends /1 or /2. If the FASTQ name
            # already carried that suffix, PAF contains name/1/1; peel the
            # duplicated mate suffix while preserving the molecule ID.
            nested_molecule_id, nested_mate = normalize_read_id(molecule_id)
            if mate is not None and nested_mate == mate:
                molecule_id = nested_molecule_id
            mate_number = mate or 1
            matches = int(fields[9])
            aligned = max(1, int(fields[10]))
            mapq = int(fields[11])
            alignment_score = matches
            for tag in fields[12:]:
                if tag.startswith("AS:i:"):
                    alignment_score = int(tag[5:])
                    break
            candidate = (alignment_score, matches / aligned, mapq)
            current = hits[molecule_id][mate_number].get(locus_id)
            if current is None or candidate > current:
                hits[molecule_id][mate_number][locus_id] = candidate
    return hits


def _short_read_recruitment_decisions(
    hits: dict[str, dict[int, dict[str, tuple[int, float, int]]]],
) -> dict[str, tuple[str, str, int, float]]:
    """Resolve mate-level native mappings into conservative molecule calls."""
    decisions: dict[str, tuple[str, str, int, float]] = {}
    for molecule_id, mate_hits in hits.items():
        confident_mates: list[str] = []
        for locus_hits in mate_hits.values():
            ranked = sorted(
                locus_hits.items(), key=lambda item: (-item[1][0], item[0])
            )
            if ranked and (
                len(ranked) == 1 or ranked[0][1][0] > ranked[1][1][0]
            ):
                confident_mates.append(ranked[0][0])
        if len(set(confident_mates)) > 1:
            decisions[molecule_id] = ("discordant", "", 0, 0.0)
            continue

        combined_scores: Counter[str] = Counter()
        identities: dict[str, list[float]] = defaultdict(list)
        for locus_hits in mate_hits.values():
            for locus_id, (score, identity, _mapq) in locus_hits.items():
                combined_scores[locus_id] += score
                identities[locus_id].append(identity)
        ranked = sorted(combined_scores.items(), key=lambda item: (-item[1], item[0]))
        if not ranked or (len(ranked) > 1 and ranked[0][1] == ranked[1][1]):
            decisions[molecule_id] = ("ambiguous", "", 0, 0.0)
            continue
        locus_id, score = ranked[0]
        identity = statistics.mean(identities[locus_id])
        decisions[molecule_id] = ("unique", locus_id, score, identity)
    return decisions


def _collect_native_recruits(
    paired_reads: tuple[Path, Path] | None,
    single_reads: Path | None,
    decisions: dict[str, tuple[str, str, int, float]],
) -> tuple[list[RecruitedPair], Counter[str]]:
    recruited: list[RecruitedPair] = []
    outcomes: Counter[str] = Counter()

    def collect(pairs) -> None:
        for pair in pairs:
            outcome, locus_id, score, identity = decisions.get(
                pair.molecule_id, ("unassigned", "", 0, 0.0)
            )
            outcomes[outcome] += 1
            if outcome == "unique":
                recruited.append(
                    RecruitedPair(pair, outcome, locus_id, score, identity)
                )

    if paired_reads is not None:
        collect(read_fastq_pairs(*paired_reads))
    if single_reads is not None:
        collect(read_fastq_pairs(single_reads))
    return recruited, outcomes


def recruit_filtered_short_reads_native(
    references: dict[str, str],
    paired_reads: tuple[Path, Path] | None,
    single_reads: Path | None,
    workdir: Path,
    threads: int,
    minimap2_bin: str = "minimap2",
) -> tuple[list[RecruitedPair], Counter[str]]:
    """Recruit filtered Illumina reads with native, multithreaded minimap2."""
    workdir.mkdir(parents=True, exist_ok=True)
    reference_path = workdir / "references.fasta"
    reference_loci = {
        f"mlva_ref_{index + 1:04d}": locus_id
        for index, locus_id in enumerate(sorted(references))
        if references[locus_id]
    }
    write_fasta(
        (
            (reference_name, references[locus_id])
            for reference_name, locus_id in reference_loci.items()
        ),
        reference_path,
    )
    if not reference_loci:
        return _collect_native_recruits(paired_reads, single_reads, {})

    executable = check_minimap2(minimap2_bin)
    all_hits: dict[str, dict[int, dict[str, tuple[int, float, int]]]] = {}
    inputs = []
    if paired_reads is not None:
        inputs.append((paired_reads[0], paired_reads[1], workdir / "paired.paf"))
    if single_reads is not None:
        inputs.append((single_reads, None, workdir / "single.paf"))
    for reads1_path, reads2_path, paf_path in inputs:
        _run_short_read_minimap2(
            _short_read_minimap2_command(
                executable,
                reference_path,
                reads1_path,
                reads2_path,
                threads,
            ),
            paf_path,
        )
        parsed = _parse_short_read_paf(paf_path, reference_loci)
        for molecule_id, mate_hits in parsed.items():
            all_hits[molecule_id] = mate_hits
    decisions = _short_read_recruitment_decisions(all_hits)
    return _collect_native_recruits(
        paired_reads, single_reads, decisions
    )


def merge_read_pair(
    pair: ReadPair,
    min_overlap: int = 20,
    max_mismatch_fraction: float = 0.03,
) -> ReadRecord | None:
    if pair.read2 is None:
        return None
    left = pair.read1
    right_sequence = revcomp(pair.read2.sequence)
    right_quality = pair.read2.quality[::-1] if pair.read2.quality else "I" * len(right_sequence)
    left_quality = left.quality or "I" * len(left.sequence)
    best: tuple[int, int] | None = None
    maximum = min(len(left.sequence), len(right_sequence))
    for overlap in range(maximum, min_overlap - 1, -1):
        mismatches = sum(
            a != b for a, b in zip(left.sequence[-overlap:], right_sequence[:overlap])
        )
        if mismatches / overlap <= max_mismatch_fraction:
            best = overlap, mismatches
            break
    if best is None:
        return None
    overlap, _mismatches = best
    merged_sequence = list(left.sequence)
    merged_quality = list(left_quality)
    start = len(left.sequence) - overlap
    for offset in range(overlap):
        index = start + offset
        if merged_sequence[index] != right_sequence[offset] and right_quality[offset] > merged_quality[index]:
            merged_sequence[index] = right_sequence[offset]
        merged_quality[index] = max(merged_quality[index], right_quality[offset])
    merged_sequence.extend(right_sequence[overlap:])
    merged_quality.extend(right_quality[overlap:])
    return ReadRecord(
        f"{pair.molecule_id}/merged",
        "".join(merged_sequence),
        "".join(merged_quality),
    )


def _merge_overlap_is_repeat_only(
    pair: ReadPair, merged: ReadRecord, locus: Locus
) -> bool:
    if pair.read2 is None or not locus.repeat_motif:
        return False
    overlap = len(pair.read1.sequence) + len(pair.read2.sequence) - len(merged.sequence)
    if overlap <= 0:
        return False
    overlap_sequence = pair.read1.sequence[-overlap:].upper()
    motif = locus.repeat_motif.upper()
    if not motif or set(motif) - set("ACGT"):
        return False
    periodic = motif * math.ceil((overlap + len(motif) * 2) / len(motif))
    return overlap_sequence in periodic


@lru_cache(maxsize=None)
def check_skesa(executable: str = "skesa") -> str:
    """Resolve the required SKESA executable or fail before processing reads."""
    resolved = shutil.which(executable)
    if resolved is None:
        raise RuntimeError(
            f"SKESA executable {executable!r} was not found. Install the skesa "
            "Bioconda package or pass --skesa-bin with an executable path."
        )
    try:
        completed = subprocess.run(
            [resolved, "--version"],
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError as exc:
        raise RuntimeError(f"Unable to execute SKESA at {resolved!r}: {exc}") from exc
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        raise RuntimeError(
            f"SKESA at {resolved!r} failed its version check"
            + (f": {detail}" if detail else "")
        )
    return resolved


def _skesa_command(
    executable: str,
    paired_reads: tuple[Path, Path] | None,
    single_reads: Path | None,
    contigs_path: Path,
    cores: int,
) -> list[str]:
    command = [
        executable,
        "--cores",
        str(max(1, cores)),
        "--memory",
        "1",
        "--hash_count",
        "--skip_bloom_filter",
        "--estimated_kmers",
        "1",
        # These are already locus-enriched reads. SKESA's default vector
        # detector would mistake shared target 19-mers for adapter sequence
        # because they occur in far more than 5% of this small read set.
        "--vector_percent",
        "1",
        "--min_contig",
        "40",
        "--contigs_out",
        str(contigs_path),
    ]
    if paired_reads is not None:
        command.extend(["--reads", f"{paired_reads[0]},{paired_reads[1]}"])
    if single_reads is not None:
        command.extend(["--reads", str(single_reads)])
    return command


def _assemble_one_locus(
    locus: Locus,
    assignments: list[RecruitedPair],
    reference: str,
    executable: str,
    workdir: Path,
    cores: int,
) -> LocusAssembly:
    locus_id = locus.locus_id
    merged = [merged for item in assignments if (merged := merge_read_pair(item.pair)) is not None]
    sequences = [read.sequence for item in assignments for read in (item.pair.read1, item.pair.read2) if read is not None]
    if not sequences:
        return LocusAssembly(locus_id, (), tuple(merged), 0.0, "NO_READS", "no uniquely recruited reads")

    workdir.mkdir(parents=True, exist_ok=True)
    complete_pairs = [item.pair for item in assignments if item.pair.read2 is not None]
    single_reads = [item.pair.read1 for item in assignments if item.pair.read2 is None]
    paired_paths: tuple[Path, Path] | None = None
    if complete_pairs:
        reads1_path = workdir / "reads_1.fastq.gz"
        reads2_path = workdir / "reads_2.fastq.gz"
        write_fastq((pair.read1 for pair in complete_pairs), reads1_path)
        write_fastq(
            (pair.read2 for pair in complete_pairs if pair.read2 is not None),
            reads2_path,
        )
        paired_paths = reads1_path, reads2_path
    singles_path: Path | None = None
    if single_reads:
        singles_path = workdir / "single_reads.fastq.gz"
        write_fastq(single_reads, singles_path)

    contigs_path = workdir / "contigs.fasta"
    command = _skesa_command(
        executable,
        paired_paths,
        singles_path,
        contigs_path,
        cores,
    )
    try:
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError as exc:
        raise RuntimeError(f"SKESA failed to start for locus {locus_id!r}: {exc}") from exc
    (workdir / "skesa.log").write_text(
        "command: " + " ".join(command) + "\n\nstdout:\n" + completed.stdout
        + "\n\nstderr:\n" + completed.stderr
    )
    compressed_contigs_path = (
        gzip_output_file(contigs_path) if contigs_path.exists() else None
    )
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout).strip()
        raise RuntimeError(
            f"SKESA failed for locus {locus_id!r} with exit code "
            f"{completed.returncode}" + (f": {detail}" if detail else "")
        )
    contigs = (
        sorted(
            {
                sequence
                for _name, sequence in read_fasta(compressed_contigs_path)
            },
            key=lambda sequence: (-len(sequence), sequence),
        )
        if compressed_contigs_path is not None
        else []
    )
    total_bp = sum(map(len, sequences))
    depth = total_bp / max(len(reference), 1)
    if not contigs:
        status = "NO_CONTIGS"
        reason = "SKESA produced no contigs"
    elif len(contigs) > 1:
        status = "AMBIGUOUS_GRAPH"
        reason = "SKESA produced multiple contigs; assembly cannot create an exact call"
    else:
        status = "ASSEMBLED"
        reason = ""
    return LocusAssembly(locus_id, tuple(contigs), tuple(merged), round(depth, 3), status, reason)


def assemble_recruited_loci(
    recruited: list[RecruitedPair],
    loci: list[Locus],
    references: dict[str, str],
    threads: int,
    workdir: Path,
    skesa_bin: str = "skesa",
) -> dict[str, LocusAssembly]:
    by_locus: dict[str, list[RecruitedPair]] = defaultdict(list)
    for item in recruited:
        if item.outcome == "unique":
            by_locus[item.locus_id].append(item)
    executable = check_skesa(skesa_bin)
    items = [
        (
            locus,
            by_locus.get(locus.locus_id, []),
            references.get(locus.locus_id, ""),
            executable,
            workdir / f"locus_{index + 1:04d}",
        )
        for index, locus in enumerate(loci)
    ]
    total_threads = resolve_threads(threads)
    loci_with_reads = sum(bool(item[1]) for item in items)
    workers = min(4, total_threads, max(1, loci_with_reads))
    cores_per_process = max(1, total_threads // workers)
    if workers == 1:
        assembled = [_assemble_one_locus(*item, cores_per_process) for item in items]
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            assembled = list(
                executor.map(
                    lambda item: _assemble_one_locus(*item, cores_per_process),
                    items,
                )
            )
    return {item.locus_id: item for item in assembled}


def estimate_insert_size_distribution(
    recruited: list[RecruitedPair], references: dict[str, str]
) -> tuple[float | None, float | None, int]:
    """Estimate fragment spans from exact, oppositely oriented locus mappings."""
    spans: list[int] = []
    for item in recruited:
        if item.outcome != "unique" or item.pair.read2 is None:
            continue
        reference = references.get(item.locus_id, "")
        if not reference:
            continue
        positions: list[tuple[int, bool, int]] = []
        for read in (item.pair.read1, item.pair.read2):
            forward = reference.find(read.sequence)
            reverse_sequence = revcomp(read.sequence)
            reverse = reference.find(reverse_sequence)
            if forward >= 0:
                positions.append((forward, False, len(read.sequence)))
            elif reverse >= 0:
                positions.append((reverse, True, len(read.sequence)))
        if len(positions) != 2 or positions[0][1] == positions[1][1]:
            continue
        start = min(position for position, _reverse, _length in positions)
        end = max(position + length for position, _reverse, length in positions)
        if end > start:
            spans.append(end - start)
    if not spans:
        return None, None, 0
    median = float(statistics.median(spans))
    mad = float(statistics.median(abs(span - median) for span in spans))
    return median, mad, len(spans)


def _boundary_patterns(locus: Locus) -> tuple[str, str]:
    left = locus.left_flank_sequence or locus.forward_primer
    right = locus.right_flank_sequence or revcomp(locus.reverse_primer)
    return left, right


def _boundary_support(sequence: str, locus: Locus) -> tuple[bool, bool]:
    left, right = _boundary_patterns(locus)
    return bool(left and find_anchor(left, sequence)), bool(right and find_anchor(right, sequence))


def _measurement_allele(measurement: LocusMeasurement) -> int | float | str | None:
    if measurement.raw_repeat_count is None or measurement.called_allele is None:
        return None
    return normalize_allele(measurement.called_allele)


def _measure_best_orientation(
    read: ReadRecord, locus: Locus, source: str
) -> tuple[str, LocusMeasurement]:
    candidates = []
    for sequence, quality in (
        (read.sequence, read.quality),
        (revcomp(read.sequence), read.quality[::-1] if read.quality else None),
    ):
        measurement = measure_locus_product(
            sequence,
            locus,
            quality,
            source=source,
            sequence_id=read.read_id,
        )
        rank = {"FULL_PRODUCT": 4, "REPEAT_INFORMATIVE": 3, "PRESENCE_ONLY": 2, "ANCHOR_FAILURE": 1}.get(measurement.status, 0)
        candidates.append((rank, measurement.confidence or 0.0, sequence, measurement))
    _rank, _confidence, sequence, measurement = max(candidates, key=lambda item: (item[0], item[1]))
    return sequence, measurement


def _call_locus(
    sample_id: str,
    locus: Locus,
    recruited: list[RecruitedPair],
    assembly: LocusAssembly,
    min_depth: int,
    insert_size: tuple[float | None, float | None, int] = (None, None, 0),
) -> dict:
    locus_pairs = [item for item in recruited if item.outcome == "unique" and item.locus_id == locus.locus_id]
    measurements: list[tuple[str, str, LocusMeasurement]] = []
    boundary1 = boundary2 = both = 0
    paired_boundaries = 0
    for item in locus_pairs:
        pair_left = pair_right = False
        for read in (item.pair.read1, item.pair.read2):
            if read is None:
                continue
            oriented_sequence, measurement = _measure_best_orientation(
                read, locus, "illumina_single_read"
            )
            left, right = _boundary_support(oriented_sequence, locus)
            boundary1 += int(left)
            boundary2 += int(right)
            both += int(left and right)
            pair_left |= left
            pair_right |= right
            measurements.append((item.pair.molecule_id, "original_single_read", measurement))
        paired_boundaries += int(
            item.pair.read2 is not None and pair_left and pair_right
        )
    for merged in assembly.merged_reads:
        source_pair = next(
            (
                item.pair
                for item in locus_pairs
                if item.pair.molecule_id == merged.read_id.removesuffix("/merged")
            ),
            None,
        )
        if source_pair is not None and _merge_overlap_is_repeat_only(
            source_pair, merged, locus
        ):
            continue
        _oriented, measurement = _measure_best_orientation(
            merged, locus, "illumina_merged_pair"
        )
        measurements.append((merged.read_id.removesuffix("/merged"), "merged_pair", measurement))
    if assembly.status == "ASSEMBLED":
        for index, contig in enumerate(assembly.contigs):
            measurement = measure_locus_product(contig, locus, source="illumina_local_assembly", sequence_id=f"{locus.locus_id}_contig_{index + 1}")
            measurements.append((f"contig:{index + 1}", "local_assembly", measurement))

    exact: list[tuple[str, str, LocusMeasurement]] = []
    partial = []
    for molecule_id, source, measurement in measurements:
        if measurement.status in {"FULL_PRODUCT", "REPEAT_INFORMATIVE"} and _measurement_allele(measurement) is not None:
            exact.append((molecule_id, source, measurement))
        elif measurement.status == "PRESENCE_ONLY":
            partial.append((molecule_id, source, measurement))
    allele_molecules: dict[int | float | str, set[str]] = defaultdict(set)
    sources_by_allele: dict[int | float | str, set[str]] = defaultdict(set)
    measurement_by_allele: dict[int | float | str, LocusMeasurement] = {}
    for molecule_id, source, measurement in exact:
        allele = _measurement_allele(measurement)
        if allele is None:
            continue
        allele_molecules.setdefault(allele, set())
        if source != "local_assembly":
            allele_molecules[allele].add(molecule_id)
        sources_by_allele[allele].add(source)
        measurement_by_allele.setdefault(allele, measurement)
    ranked = sorted(allele_molecules, key=lambda allele: (-len(allele_molecules[allele]), float(allele)))
    primary = ranked[0] if ranked else None
    primary_support = len(allele_molecules.get(primary, set())) if primary is not None else 0
    secondary_support = sum(len(allele_molecules[allele]) for allele in ranked[1:])
    informative = len({molecule for molecules in allele_molecules.values() for molecule in molecules})
    uninformative = max(0, len(locus_pairs) - informative)
    mixed = len(ranked) > 1
    effective_informative = (
        max(informative, len(locus_pairs))
        if primary is not None and "local_assembly" in sources_by_allele[primary]
        else informative
    )
    if mixed:
        evidence_class = "MULTIPLE_ALLELES"
    elif primary is not None and "local_assembly" in sources_by_allele[primary]:
        evidence_class = "COMPLETE_ASSEMBLED_PRODUCT"
    elif primary is not None and "merged_pair" in sources_by_allele[primary]:
        evidence_class = "BOUNDARY_SPANNING_READ_PAIR"
    elif primary is not None:
        evidence_class = "BOUNDARY_SPANNING_SINGLE_READ"
    elif paired_boundaries:
        evidence_class = "BOUNDARY_SPANNING_READ_PAIR"
    elif len(assembly.contigs) > 1 and partial:
        evidence_class = "AMBIGUOUS_ASSEMBLY"
    elif partial:
        evidence_class = "PARTIAL_REPEAT_EVIDENCE"
    elif locus_pairs:
        evidence_class = "PRESENCE_ONLY"
    else:
        evidence_class = "NOT_FOUND"
    if primary is not None and effective_informative < min_depth and not mixed:
        displayed_evidence_class = "LOW_DEPTH"
    else:
        displayed_evidence_class = evidence_class

    repeat_min: int | float | str = ""
    repeat_max: int | float | str = ""
    interval_reason = ""
    if primary is not None:
        repeat_min = repeat_max = primary
    elif paired_boundaries:
        median_insert, mad_insert, insert_support = insert_size
        unit = repeat_unit_length(locus)
        reference_fixed_bp = (
            len(locus.forward_primer)
            + len(locus.left_flank_sequence)
            + len(locus.right_flank_sequence)
            + len(locus.reverse_primer)
        )
        if median_insert is not None and unit > 0 and insert_support >= 2:
            center = (median_insert - reference_fixed_bp) / unit
            uncertainty = max(1, math.ceil(3 * max(mad_insert or 0.0, 1.0) / unit))
            repeat_min = max(locus.expected_min_repeats, math.floor(center - uncertainty))
            repeat_max = min(locus.expected_max_repeats, math.ceil(center + uncertainty))
            if float(repeat_min) > float(repeat_max):
                repeat_min, repeat_max = locus.expected_min_repeats, locus.expected_max_repeats
            interval_reason = f"opposite boundaries plus empirical insert span (median {median_insert:.1f} bp, MAD {mad_insert or 0.0:.1f}, n={insert_support})"
        else:
            repeat_min = locus.expected_min_repeats
            repeat_max = locus.expected_max_repeats
            interval_reason = "opposite repeat boundaries observed in a pair; insert span does not uniquely determine an allele"
    confidence = 0.0
    confidence_reason = ""
    if primary is not None:
        confidence_support = (
            max(primary_support, len(locus_pairs))
            if "local_assembly" in sources_by_allele[primary]
            else primary_support
        )
        confidence = 0.99 if "local_assembly" in sources_by_allele[primary] and confidence_support >= min_depth else 0.9 if primary_support >= min_depth else 0.6
        confidence_reason = (
            "primer/flank-bounded local assembly supported by informative molecules"
            if "local_assembly" in sources_by_allele[primary]
            else "repeat boundaries observed on one or more sequenced molecules"
        )
        if mixed:
            confidence = min(confidence, 0.75)
            confidence_reason += "; conflicting defensible alleles retained"
    elif paired_boundaries:
        confidence_reason = "pair spans both boundaries but repeat length is not directly observed"
    elif locus_pairs:
        confidence_reason = "locus recruited without both repeat boundaries"
    else:
        confidence_reason = "no unique locus recruitment"
    warning = ""
    if primary is None and locus_pairs:
        warning = "short reads do not contain enough boundary-spanning sequence for an exact repeat count"
    total_support = primary_support + secondary_support
    primary_fraction = primary_support / total_support if total_support else ""
    secondary_fraction = secondary_support / total_support if total_support else ""
    product_measurement = measurement_by_allele.get(primary) if primary is not None else None
    present = bool(locus_pairs or primary is not None)
    status = (
        "MULTIPLE_ALLELES" if mixed else
        "LOW_DEPTH" if primary is not None and effective_informative < min_depth else
        "PASS" if primary is not None else
        "INTERVAL_ONLY" if paired_boundaries else
        "PRESENCE_ONLY" if present else
        "NOT_FOUND"
    )
    secondary = ";".join(f"{allele}:{len(allele_molecules[allele])}" for allele in ranked[1:])
    distribution = ";".join(
        f"{allele}:{len(allele_molecules[allele]) / max(total_support, 1):.6f}"
        for allele in ranked
    )
    return {
        "sample_id": sample_id,
        "locus_id": locus.locus_id,
        "present": "yes" if present else "no",
        "repeat_count": "" if primary is None else primary,
        "repeat_count_raw": "" if product_measurement is None else product_measurement.raw_repeat_count,
        "product_size_bp": "" if product_measurement is None or not product_measurement.product_sequence else len(product_measurement.product_sequence),
        "read_depth": len(locus_pairs),
        "primary_read_depth": primary_support,
        "mean_coverage": assembly.estimated_depth if locus_pairs else "",
        "allele_confidence": confidence,
        "second_best_repeat_count": ranked[1] if len(ranked) > 1 else "",
        "second_best_probability": secondary_fraction,
        "inference_method": "short_read_boundary_evidence",
        "dominant_variant_fraction": primary_fraction,
        "num_candidate_variants": len(ranked),
        "num_confirmed_secondary_variants": max(0, len(ranked) - 1),
        "secondary_alleles": secondary,
        "allele_distribution": distribution,
        "status": status,
        "evidence": evidence_class,
        "read_technology": "illumina",
        "evidence_class": displayed_evidence_class,
        "informative_molecule_count": informative,
        "boundary_1_support": boundary1,
        "boundary_2_support": boundary2,
        "both_boundary_support": both,
        "local_assembly_status": assembly.status,
        "repeat_count_min": repeat_min,
        "repeat_count_max": repeat_max,
        "repeat_count_interval_reason": interval_reason,
        "confidence_reason": confidence_reason,
        "short_read_warning": warning,
        "recruited_read_pairs": len(locus_pairs),
        "assembled_contig_length": max(map(len, assembly.contigs), default=0),
        "assembly_depth": assembly.estimated_depth,
        "failure_reason": assembly.failure_reason if not present else "",
        "primary_allele_support": primary_support,
        "secondary_allele_support": secondary_support,
        "informative_allele_reads": informative,
        "uninformative_locus_reads": uninformative,
        "estimated_primary_fraction": primary_fraction,
        "estimated_secondary_fraction": secondary_fraction,
        "mixture_status": "mixed" if mixed else "single" if primary is not None else "unavailable",
        "evidence_sources": ",".join(sorted(sources_by_allele.get(primary, set()))) if primary is not None else "paired_end_span" if paired_boundaries else "recruitment",
    }


def _allele_rows(call_rows: list[dict]) -> list[dict]:
    return [
        {
            "sample_id": row["sample_id"],
            "locus_id": row["locus_id"],
            "called_repeat_count": row["repeat_count"],
            "posterior_probability": row["allele_confidence"],
            "second_best_repeat_count": row["second_best_repeat_count"],
            "second_best_posterior": row["second_best_probability"],
            "read_depth": row["read_depth"],
            "primary_read_depth": row["primary_read_depth"],
            "num_vntr_asvs": row["num_candidate_variants"],
            "num_meaningful_variants": 1 + row["num_confirmed_secondary_variants"] if row["repeat_count"] != "" else 0,
            "num_candidate_variants": row["num_candidate_variants"],
            "num_confirmed_secondary_variants": row["num_confirmed_secondary_variants"],
            "dominant_vntr_asv": row["repeat_count"],
            "dominant_variant_fraction": row["estimated_primary_fraction"],
            "secondary_alleles": row["secondary_alleles"],
            "allele_distribution": row["allele_distribution"],
            "call_status": row["status"],
            "primary_product_size_bp": row["product_size_bp"],
            "primary_repeat_count_raw": row["repeat_count_raw"],
            "primary_measurement_source": row["evidence_sources"],
            "evidence_status": row["evidence_class"],
        }
        for row in call_rows
    ]


def recruitment_summary(
    sample_id: str,
    loci: list[Locus],
    recruited: list[RecruitedPair],
    read_pairs_examined: int | None = None,
    outcome_counts: Counter[str] | None = None,
) -> list[dict]:
    outcome_counts = outcome_counts or Counter(item.outcome for item in recruited)
    examined = len(recruited) if read_pairs_examined is None else read_pairs_examined
    rows = []
    for locus in loci:
        unique = [item for item in recruited if item.outcome == "unique" and item.locus_id == locus.locus_id]
        identities = [item.identity for item in unique]
        rows.append(
            {
                "sample_id": sample_id,
                "locus_id": locus.locus_id,
                "read_pairs_examined": examined,
                "read_pairs_recruited": len(unique),
                "uniquely_recruited_pairs": len(unique),
                "ambiguous_pairs": outcome_counts["ambiguous"],
                "discordant_pairs": outcome_counts["discordant"],
                "orphan_reads": sum(item.pair.read2 is None for item in unique),
                "mean_mapping_quality": round(60 * statistics.mean(identities), 3) if identities else "",
                "mean_alignment_identity": round(statistics.mean(identities), 6) if identities else "",
                "presence_status": "PRESENT_UNTYPED" if unique else "NO_EVIDENCE",
                "mapped_reads": sum(1 + int(item.pair.read2 is not None) for item in unique),
                "full_product_reads": 0,
                "genotype_informative_reads": 0,
                "candidate_alleles": "",
                "reference_source": "short_read_minimap2_sr",
            }
        )
    return rows


def run_short_read_call(
    reads1_path: str,
    reads2_path: str | None,
    loci_path: str | None,
    outdir: str,
    sample_id: str,
    primers_path: str | None = None,
    profiles_path: str | None = None,
    database_path: str | None = None,
    sample_metadata: dict[str, str] | None = None,
    short_min_read_length: int = 40,
    short_min_mean_quality: float = 15.0,
    short_trim_quality: int = 0,
    short_min_pair_retention: float = 0.5,
    min_depth: int = 3,
    threads: int = DEFAULT_THREADS,
    keep_intermediates: bool = False,
    sample_mode: str = "isolate",
    skesa_bin: str = "skesa",
    minimap2_bin: str = "minimap2",
    show_progress: bool = True,
) -> dict[str, Path]:
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    resolved_skesa = check_skesa(skesa_bin)
    resolved_minimap2 = check_minimap2(minimap2_bin)
    loci = read_loci_or_primers(loci_path, primers_path)
    profiles = read_profiles(profiles_path)
    references = build_short_read_references(loci, database_path)
    qc_counter: Counter[str] = Counter()
    filtered1_path = outdir_path / "filtered_reads_1.fastq.gz"
    filtered2_path = outdir_path / "filtered_reads_2.fastq.gz"
    filtered_orphans_path = outdir_path / "filtered_orphan_reads.fastq.gz"

    def write_record(handle, read: ReadRecord) -> None:
        quality = read.quality or "I" * len(read.sequence)
        handle.write(f"@{read.read_id}\n{read.sequence}\n+\n{quality}\n")

    def process_chunk(
        chunk: list[ReadPair],
        handle1,
        handle2,
        orphan_handle,
    ) -> tuple[int, int]:
        retained, metrics = qc_read_pairs(
            chunk,
            short_min_read_length,
            short_min_mean_quality,
            short_trim_quality,
            short_min_pair_retention,
        )
        qc_counter.update(metrics)
        for pair in retained:
            if reads2_path and pair.read2 is None:
                write_record(
                    orphan_handle,
                    ReadRecord(
                        pair.molecule_id,
                        pair.read1.sequence,
                        pair.read1.quality,
                    ),
                )
            else:
                write_record(handle1, pair.read1)
                if handle2 is not None and pair.read2 is not None:
                    write_record(handle2, pair.read2)
        return (
            sum(pair.read2 is not None for pair in retained),
            sum(pair.read2 is None for pair in retained),
        )

    def filter_and_recruit(recruitment_workdir: Path):
        filter_started = time.perf_counter()
        if show_progress:
            print(f"[{sample_id}] Filtering and validating Illumina read pairs")
        recruitment_workdir.mkdir(parents=True, exist_ok=True)
        paired_count = 0
        single_count = 0
        with ExitStack() as stack:
            handle1 = stack.enter_context(open_text(filtered1_path, "wt"))
            handle2 = (
                stack.enter_context(open_text(filtered2_path, "wt"))
                if reads2_path
                else None
            )
            orphan_handle = stack.enter_context(
                open_text(filtered_orphans_path, "wt")
            )
            chunk: list[ReadPair] = []
            for pair in read_fastq_pairs(reads1_path, reads2_path):
                chunk.append(pair)
                if len(chunk) >= 5000:
                    new_paired, new_single = process_chunk(
                        chunk,
                        handle1,
                        handle2,
                        orphan_handle,
                    )
                    paired_count += new_paired
                    single_count += new_single
                    chunk.clear()
            if chunk:
                new_paired, new_single = process_chunk(
                    chunk,
                    handle1,
                    handle2,
                    orphan_handle,
                )
                paired_count += new_paired
                single_count += new_single
        if show_progress:
            elapsed = time.perf_counter() - filter_started
            print(
                f"[{sample_id}] QC retained {paired_count:,} pairs and "
                f"{single_count:,} single/orphan reads in {elapsed:.1f}s"
            )
            print(
                f"[{sample_id}] Recruiting reads with minimap2 "
                f"({resolve_threads(threads)} threads)"
            )
        recruitment_started = time.perf_counter()
        result = recruit_filtered_short_reads_native(
            references,
            (filtered1_path, filtered2_path) if reads2_path and paired_count else None,
            (
                filtered_orphans_path
                if reads2_path and single_count
                else filtered1_path
                if not reads2_path and single_count
                else None
            ),
            recruitment_workdir,
            threads,
            resolved_minimap2,
        )
        if show_progress:
            elapsed = time.perf_counter() - recruitment_started
            print(
                f"[{sample_id}] minimap2 recruited {len(result[0]):,} "
                f"molecules in {elapsed:.1f}s"
            )
        return result

    if keep_intermediates:
        recruitment_workdir = outdir_path / "short_read_recruitment_intermediates"
        recruited, outcome_counts = filter_and_recruit(recruitment_workdir)
        reference_file = recruitment_workdir / "references.fasta"
        if reference_file.exists():
            gzip_output_file(reference_file)
    else:
        with tempfile.TemporaryDirectory(
            prefix=".mlvamaps-recruitment-", dir=outdir_path
        ) as temporary_recruitment:
            recruited, outcome_counts = filter_and_recruit(
                Path(temporary_recruitment)
            )
    qc = dict(qc_counter)
    qc_rows = [{"sample_id": sample_id, "metric": metric, "value": value} for metric, value in sorted(qc.items())]
    write_tsv(qc_rows, outdir_path / "short_read_qc_summary.tsv", SHORT_QC_FIELDS)
    write_tsv(qc_rows, outdir_path / "qc_summary.tsv", SHORT_QC_FIELDS)
    recruitment_rows = recruitment_summary(
        sample_id,
        loci,
        recruited,
        read_pairs_examined=qc.get("retained_pairs", 0),
        outcome_counts=outcome_counts,
    )
    write_tsv(recruitment_rows, outdir_path / "short_read_recruitment_summary.tsv", SHORT_RECRUITMENT_FIELDS)
    write_tsv(recruitment_rows, outdir_path / "locus_presence.tsv", SHORT_RECRUITMENT_FIELDS)
    if show_progress:
        loci_with_reads = len({item.locus_id for item in recruited})
        print(
            f"[{sample_id}] Assembling {loci_with_reads:,} recruited loci "
            f"with parallel SKESA"
        )
    assembly_started = time.perf_counter()
    if keep_intermediates:
        assembly_workdir = outdir_path / "short_read_assembly_intermediates"
        assembly_workdir.mkdir(parents=True, exist_ok=True)
        assemblies = assemble_recruited_loci(
            recruited,
            loci,
            references,
            threads,
            assembly_workdir,
            resolved_skesa,
        )
    else:
        with tempfile.TemporaryDirectory(
            prefix=".mlvamaps-skesa-", dir=outdir_path
        ) as temporary_workdir:
            assemblies = assemble_recruited_loci(
                recruited,
                loci,
                references,
                threads,
                Path(temporary_workdir),
                resolved_skesa,
            )
    if show_progress:
        print(
            f"[{sample_id}] SKESA assembly finished in "
            f"{time.perf_counter() - assembly_started:.1f}s"
        )
    insert_size = estimate_insert_size_distribution(recruited, references)
    if insert_size[0] is not None:
        qc_rows.extend(
            [
                {"sample_id": sample_id, "metric": "empirical_insert_size_median_bp", "value": insert_size[0]},
                {"sample_id": sample_id, "metric": "empirical_insert_size_mad_bp", "value": insert_size[1]},
                {"sample_id": sample_id, "metric": "empirical_insert_size_pairs", "value": insert_size[2]},
            ]
        )
        write_tsv(qc_rows, outdir_path / "short_read_qc_summary.tsv", SHORT_QC_FIELDS)
        write_tsv(qc_rows, outdir_path / "qc_summary.tsv", SHORT_QC_FIELDS)
    contig_records = [
        (f"{locus_id}_contig_{index + 1}", sequence)
        for locus_id, assembly in sorted(assemblies.items())
        for index, sequence in enumerate(assembly.contigs)
    ]
    write_fasta(contig_records, outdir_path / "local_locus_products.fasta.gz")
    assembly_rows = []
    for locus in loci:
        assembly = assemblies[locus.locus_id]
        locus_items = [item for item in recruited if item.outcome == "unique" and item.locus_id == locus.locus_id]
        assembly_rows.append(
            {
                "sample_id": sample_id,
                "locus_id": locus.locus_id,
                "recruited_reads": sum(1 + int(item.pair.read2 is not None) for item in locus_items),
                "recruited_pairs": len(locus_items),
                "assembler": assembly.assembler,
                "assembly_status": assembly.status,
                "contig_count": len(assembly.contigs),
                "longest_contig_bp": max(map(len, assembly.contigs), default=0),
                "total_contig_bp": sum(map(len, assembly.contigs)),
                "estimated_depth": assembly.estimated_depth,
                "failure_reason": assembly.failure_reason,
            }
        )
    write_tsv(assembly_rows, outdir_path / "short_read_assembly_summary.tsv", SHORT_ASSEMBLY_FIELDS)
    write_tsv(assembly_rows, outdir_path / "local_assembly_concordance.tsv", SHORT_ASSEMBLY_FIELDS)

    call_rows = [
        _call_locus(sample_id, locus, recruited, assemblies[locus.locus_id], min_depth, insert_size)
        for locus in loci
    ]
    calls_path = outdir_path / "calls.tsv"
    write_tsv(call_rows, calls_path, SHORT_CALL_FIELDS)
    write_tsv(call_rows, outdir_path / "locus_repeat_counts.tsv", REPEAT_COUNT_FIELDS + [field for field in SHORT_CALL_EXTRA_FIELDS if field not in REPEAT_COUNT_FIELDS])
    allele_rows = _allele_rows(call_rows)
    write_tsv(
        allele_distribution_rows(sample_id, call_rows),
        outdir_path / "allele_probability_distribution.tsv",
        ALLELE_DISTRIBUTION_FIELDS,
    )
    fingerprint_rows, probabilistic_rows = build_fingerprint(sample_id, allele_rows, loci)
    write_tsv(fingerprint_rows, outdir_path / "mlva_fingerprint.tsv", ["sample_id"] + [locus.locus_id for locus in loci])
    write_tsv(probabilistic_rows, outdir_path / "mlva_fingerprint_probabilistic.tsv", ["sample_id", "locus_id", "repeat_count", "posterior_probability"])
    matches = match_profiles(sample_id, fingerprint_rows[0], profiles, allele_rows=allele_rows)
    write_tsv(matches, outdir_path / "profile_matches.tsv", MATCH_FIELDS)
    write_tsv(
        profile_match_locus_rows(sample_id, fingerprint_rows[0], profiles, matches, allele_rows),
        outdir_path / "profile_match_loci.tsv",
        PROFILE_MATCH_LOCUS_FIELDS,
    )
    best = matches[0] if matches else {}
    complete_classes = {"COMPLETE_ASSEMBLED_PRODUCT", "BOUNDARY_SPANNING_READ_PAIR", "BOUNDARY_SPANNING_SINGLE_READ", "MULTIPLE_ALLELES", "LOW_DEPTH"}
    summary = {
        "sample_id": sample_id,
        "input_read_1": str(Path(reads1_path)),
        "input_read_2": "" if reads2_path is None else str(Path(reads2_path)),
        "read_technology": "illumina",
        "sample_mode": sample_mode,
        "total_reads": qc.get("input_reads", 0),
        "total_read_pairs": qc.get("input_pairs", 0),
        "retained_reads": qc.get("retained_reads", 0),
        "retained_pairs": qc.get("retained_pairs", 0),
        "callable_loci": sum(row["repeat_count"] != "" for row in call_rows),
        "complete_loci": sum(row["evidence_class"] in complete_classes and row["repeat_count"] != "" for row in call_rows),
        "partial_loci": sum(row["evidence_class"] in {"PARTIAL_REPEAT_EVIDENCE", "BOUNDARY_SPANNING_READ_PAIR"} and row["repeat_count"] == "" for row in call_rows),
        "presence_only_loci": sum(row["evidence_class"] == "PRESENCE_ONLY" for row in call_rows),
        "mixed_loci": sum(row["evidence_class"] == "MULTIPLE_ALLELES" for row in call_rows),
        "missing_loci": sum(row["evidence_class"] == "NOT_FOUND" for row in call_rows),
        "best_profile_id": best.get("best_profile_id", ""),
        "best_profile_distance": best.get("distance", ""),
        "profile_confidence": best.get("confidence", ""),
        "run_status": "success",
        "warnings": ";".join(sorted({row["short_read_warning"] for row in call_rows if row["short_read_warning"]})),
    }
    write_tsv([summary], outdir_path / "sample_summary.tsv", SAMPLE_SUMMARY_FIELDS)
    myoga_row = myoga_sample_row(sample_id, sample_metadata, summary, fingerprint_rows[0], len(loci))
    write_csv([myoga_row], outdir_path / "myoga_samples.csv", MYOGA_SAMPLE_FIELDS)
    write_csv(
        [
            {
                "genome_id": sample_id,
                "sample_id": sample_id,
                "locus_id": row["locus_id"],
                "repeat_count": row["repeat_count"],
                "repeat_count_min": row["repeat_count_min"],
                "repeat_count_max": row["repeat_count_max"],
                "evidence_class": row["evidence_class"],
                "confidence": row["allele_confidence"],
            }
            for row in call_rows
        ],
        outdir_path / "myoga_loci.csv",
        ["genome_id", "sample_id", "locus_id", "repeat_count", "repeat_count_min", "repeat_count_max", "evidence_class", "confidence"],
    )
    write_report(
        outdir_path,
        sample_id,
        allele_rows,
        loci,
        matches,
        profiles,
        presence_rows=recruitment_rows,
        local_assembly_rows=[],
        short_read_rows=call_rows,
    )
    return {
        "outdir": outdir_path,
        "calls": calls_path,
        "repeat_counts": outdir_path / "locus_repeat_counts.tsv",
        "allele_distribution": outdir_path / "allele_probability_distribution.tsv",
        "fingerprint": outdir_path / "mlva_fingerprint.tsv",
        "profile_matches": outdir_path / "profile_matches.tsv",
        "profile_match_loci": outdir_path / "profile_match_loci.tsv",
        "report": outdir_path / "report.html",
        "sample_summary": outdir_path / "sample_summary.tsv",
        "myoga_samples": outdir_path / "myoga_samples.csv",
        "myoga_loci": outdir_path / "myoga_loci.csv",
        "short_read_qc": outdir_path / "short_read_qc_summary.tsv",
        "short_read_recruitment": outdir_path / "short_read_recruitment_summary.tsv",
        "short_read_assembly": outdir_path / "short_read_assembly_summary.tsv",
    }
