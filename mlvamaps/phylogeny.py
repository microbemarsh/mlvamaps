from __future__ import annotations

import csv
import gzip
import json
import math
import re
import shutil
import statistics
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from .calling import estimate_repeat_count_from_product_length, normalize_allele, repeat_unit_length
from .concurrency import resolve_threads
from .models import Locus, RepeatFeature
from .progress import ProgressReporter


PLACEMENT_FIELDS = [
    "sample_id",
    "locus_id",
    "reference_id",
    "phylogenetic_distance",
    "likelihood_weighted_phylogenetic_distance",
    "placement_edge",
    "like_weight_ratio",
    "pendant_length",
    "distal_length",
    "placement_count",
    "placement_entropy",
    "reference_tree_scale",
]

SUMMARY_FIELDS = [
    "sample_id",
    "reference_id",
    "total_phylogenetic_distance",
    "total_likelihood_weighted_distance",
    "compared_loci",
    "mean_phylogenetic_distance",
    "mean_likelihood_weighted_distance",
    "distance_gap_to_next",
    "relative_distance_gap_to_next",
    "rank",
]

STATUS_FIELDS = ["locus_id", "reference_sequences", "query_sequence", "status"]

MARKER_COMPONENT_FIELDS = [
    "sample_id",
    "locus_id",
    "record_type",
    "reference_id",
    "repeat_count_raw",
    "repeat_count",
    "repeat_sequence",
    "repeat_haplotype",
    "repeat_region_start",
    "repeat_region_end",
    "snp_sequence_length",
    "masking_method",
]

LOCUS_MARKER_DISTANCE_FIELDS = [
    "sample_id",
    "locus_id",
    "reference_id",
    "query_repeat_count",
    "reference_repeat_count",
    "repeat_count_delta",
    "normalized_repeat_distance",
    "likelihood_weighted_snp_distance",
    "reference_tree_scale",
    "normalized_snp_distance",
    "placement_entropy",
]

COMBINED_MARKER_FIELDS = [
    "sample_id",
    "reference_id",
    "total_likelihood_weighted_snp_distance",
    "total_normalized_snp_distance",
    "total_repeat_count_distance",
    "total_normalized_repeat_distance",
    "snp_weight",
    "repeat_weight",
    "combined_marker_distance",
    "compared_loci",
    "repeat_compared_loci",
    "distance_gap_to_next",
    "relative_distance_gap_to_next",
    "collection_date",
    "latitude",
    "longitude",
    "location",
    "source",
    "rank",
]

_FASTA_SUFFIXES = {".fa", ".fas", ".fasta", ".fna", ".ffn"}
_SAFE_FILE = re.compile(r"[^A-Za-z0-9_.-]+")
_IUPAC = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"),
    "T": frozenset("T"), "R": frozenset("AG"), "Y": frozenset("CT"),
    "S": frozenset("CG"), "W": frozenset("AT"), "K": frozenset("GT"),
    "M": frozenset("AC"), "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"), "N": frozenset("ACGT"),
}
_RC = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


@dataclass(frozen=True)
class MarkerComponents:
    oriented_sequence: str
    snp_sequence: str
    repeat_sequence: str
    repeat_count_raw: float | None
    repeat_count: int | float | None
    repeat_region_start: int | None
    repeat_region_end: int | None
    masking_method: str


def _reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(_RC)[::-1]


def _iupac_find(pattern: str, sequence: str, start: int = 0) -> int:
    pattern = pattern.upper()
    sequence = sequence.upper()
    if not pattern:
        return -1
    for index in range(start, len(sequence) - len(pattern) + 1):
        if all(base in _IUPAC.get(code, frozenset(code)) for code, base in zip(pattern, sequence[index:])):
            return index
    return -1


def _longest_motif_run(sequence: str, motif: str) -> tuple[int, int] | None:
    motif = motif.upper()
    if not motif or set(motif) == {"N"}:
        return None
    motif_length = len(motif)
    best: tuple[int, int] | None = None
    for offset in range(len(sequence) - motif_length + 1):
        end = offset
        while end + motif_length <= len(sequence):
            chunk = sequence[end : end + motif_length]
            mismatches = sum(
                base not in _IUPAC.get(code, frozenset(code))
                for code, base in zip(motif, chunk)
            )
            if mismatches > motif_length // 8:
                break
            end += motif_length
        if end > offset and (best is None or end - offset > best[1] - best[0]):
            best = (offset, end)
    return best


def decompose_marker_sequence(locus: Locus, sequence: str) -> MarkerComponents:
    """Separate explicit VNTR characters from sequence used for SNP placement."""
    sequence = sequence.upper().replace("-", "")
    reverse_primer_site = _reverse_complement(locus.reverse_primer)
    forward = _iupac_find(locus.forward_primer, sequence)
    reverse = _iupac_find(reverse_primer_site, sequence, max(forward, 0))
    if forward < 0 or reverse < 0 or reverse <= forward:
        reverse_sequence = _reverse_complement(sequence)
        reverse_forward = _iupac_find(locus.forward_primer, reverse_sequence)
        reverse_reverse = _iupac_find(
            reverse_primer_site, reverse_sequence, max(reverse_forward, 0)
        )
        if reverse_forward >= 0 and reverse_reverse > reverse_forward:
            sequence = reverse_sequence
            forward, reverse = reverse_forward, reverse_reverse

    inner_start = forward + len(locus.forward_primer) if forward >= 0 else 0
    inner_end = reverse if reverse > inner_start else len(sequence)
    repeat_start: int | None = None
    repeat_end: int | None = None
    method = "unmasked"
    if locus.left_flank_sequence and locus.right_flank_sequence:
        left = _iupac_find(locus.left_flank_sequence, sequence, inner_start)
        right_start = left + len(locus.left_flank_sequence) if left >= 0 else inner_start
        right = _iupac_find(locus.right_flank_sequence, sequence, right_start)
        if left >= 0 and right >= right_start and right <= inner_end:
            repeat_start = right_start
            repeat_end = right
            method = "flank_bounded"
    if repeat_start is None:
        motif_run = _longest_motif_run(sequence[inner_start:inner_end], locus.repeat_motif)
        if motif_run is not None:
            repeat_start = inner_start + motif_run[0]
            repeat_end = inner_start + motif_run[1]
            method = "motif_run"

    repeat_sequence = (
        sequence[repeat_start:repeat_end]
        if repeat_start is not None and repeat_end is not None
        else ""
    )
    unit_length = repeat_unit_length(locus)
    if repeat_sequence and unit_length:
        raw_count = len(repeat_sequence) / unit_length
    else:
        raw_count = estimate_repeat_count_from_product_length(locus, len(sequence))
    repeat_count = normalize_allele(raw_count) if raw_count is not None else None
    snp_sequence = (
        sequence[:repeat_start] + sequence[repeat_end:]
        if repeat_start is not None and repeat_end is not None
        else sequence
    )
    if not snp_sequence:
        snp_sequence = sequence
        method += "_empty_mask_fallback"
    return MarkerComponents(
        sequence,
        snp_sequence,
        repeat_sequence,
        raw_count,
        repeat_count,
        repeat_start,
        repeat_end,
        method,
    )


def _repeat_haplotype(sequence: str, unit_length: int) -> str:
    if not sequence or unit_length <= 0:
        return ""
    return "|".join(
        sequence[index : index + unit_length]
        for index in range(0, len(sequence), unit_length)
    )


def _marker_component_row(
    sample_id: str,
    locus: Locus,
    record_type: str,
    reference_id: str,
    components: MarkerComponents,
) -> dict:
    return {
        "sample_id": sample_id,
        "locus_id": locus.locus_id,
        "record_type": record_type,
        "reference_id": reference_id,
        "repeat_count_raw": ""
        if components.repeat_count_raw is None
        else f"{components.repeat_count_raw:.6f}",
        "repeat_count": "" if components.repeat_count is None else components.repeat_count,
        "repeat_sequence": components.repeat_sequence,
        "repeat_haplotype": _repeat_haplotype(
            components.repeat_sequence, repeat_unit_length(locus)
        ),
        "repeat_region_start": ""
        if components.repeat_region_start is None
        else components.repeat_region_start,
        "repeat_region_end": ""
        if components.repeat_region_end is None
        else components.repeat_region_end,
        "snp_sequence_length": len(components.snp_sequence),
        "masking_method": components.masking_method,
    }


def _read_fasta(path: str | Path) -> list[tuple[str, str]]:
    path = Path(path)
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    records: list[tuple[str, str]] = []
    name: str | None = None
    sequence: list[str] = []
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(sequence).upper()))
                name = line[1:].split()[0]
                sequence = []
            elif name is None:
                raise ValueError(f"Sequence appeared before a FASTA header in {path}")
            else:
                sequence.append(line)
    if name is not None:
        records.append((name, "".join(sequence).upper()))
    return records


def _write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")


def _write_tsv(rows: list[dict], path: Path, fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


@dataclass
class _Node:
    name: str | None
    children: list[tuple["_Node", float]]
    edge_num: int | None = None


def check_mafft(executable: str) -> str:
    path = shutil.which(executable)
    if path is None:
        raise RuntimeError(
            f"MAFFT executable {executable!r} was not found. Install mafft from "
            "Bioconda or pass --mafft-bin."
        )
    result = subprocess.run(
        [path, "--version"], capture_output=True, text=True, check=False
    )
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(f"Could not run MAFFT at {path}: {detail}")
    return path


def build_mafft_reference_command(
    reference_fasta: str | Path, threads: int, executable: str = "mafft"
) -> list[str]:
    return [executable, "--auto", "--thread", str(threads), str(reference_fasta)]


def build_mafft_add_command(
    query_fasta: str | Path,
    reference_alignment: str | Path,
    threads: int,
    executable: str = "mafft",
) -> list[str]:
    return [
        executable,
        "--add",
        str(query_fasta),
        "--keeplength",
        "--thread",
        str(threads),
        str(reference_alignment),
    ]


def _run_mafft(command: list[str], output_path: Path, stage: str) -> None:
    with output_path.open("w") as output:
        result = subprocess.run(
            command,
            stdout=output,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
    if result.returncode:
        output_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"MAFFT {stage} failed (exit {result.returncode}): "
            f"{(result.stderr or '').strip()}"
        )


def check_raxml_ng(executable: str) -> str:
    path = shutil.which(executable)
    if path is None:
        raise RuntimeError(
            f"RAxML-NG executable {executable!r} was not found. Install raxml-ng "
            "from Bioconda or pass --raxml-ng-bin."
        )
    result = subprocess.run(
        [path, "--version"], capture_output=True, text=True, check=False
    )
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(f"Could not run RAxML-NG at {path}: {detail}")
    return path


def build_raxml_ng_command(
    alignment_path: str | Path,
    prefix: str | Path,
    threads: int,
    executable: str = "raxml-ng",
    model: str = "DNA",
) -> list[str]:
    return [
        executable,
        "--search",
        "--msa",
        str(alignment_path),
        "--model",
        model,
        "--prefix",
        str(prefix),
        "--seed",
        "12345",
        "--threads",
        str(threads),
        "--redo",
    ]


def _run_raxml_ng(
    command: list[str],
    prefix: Path,
    output_tree: Path,
    stage: str,
    progress: ProgressReporter | None = None,
) -> None:
    """Run RAxML-NG, reducing only its threads for low-pattern alignments."""
    attempted_command = list(command)
    thread_index = attempted_command.index("--threads") + 1
    raxml_threads = int(attempted_command[thread_index])
    log_path = Path(f"{prefix}.mlvamaps.raxml.log")
    attempts: list[str] = []

    while True:
        attempted_command[thread_index] = str(raxml_threads)
        result = subprocess.run(
            attempted_command, capture_output=True, text=True, check=False
        )
        detail = "\n".join(
            part.strip() for part in (result.stdout, result.stderr) if part.strip()
        )
        attempts.append(
            f"$ {' '.join(attempted_command)}\n{detail}".rstrip()
        )
        log_path.write_text("\n\n".join(attempts) + "\n")
        if not result.returncode:
            break
        if "Too few patterns per thread" not in detail or raxml_threads == 1:
            detail_tail = "\n".join(detail.splitlines()[-40:])
            raise RuntimeError(
                f"RAxML-NG {stage} failed (exit {result.returncode}). Full output: "
                f"{log_path}\n{detail_tail}"
            )
        next_threads = max(1, raxml_threads // 2)
        if progress is not None:
            progress.step(
                f"RAxML-NG found too few alignment patterns for {raxml_threads} "
                f"threads; retrying {stage} with {next_threads}"
            )
        raxml_threads = next_threads

    best_tree = Path(f"{prefix}.raxml.bestTree")
    if not best_tree.exists():
        raise RuntimeError(
            f"RAxML-NG {stage} completed without producing {best_tree}"
        )
    shutil.copyfile(best_tree, output_tree)


def check_epa_ng(executable: str) -> str:
    path = shutil.which(executable)
    if path is None:
        raise RuntimeError(
            f"EPA-ng executable {executable!r} was not found. Install epa-ng "
            "from Bioconda or pass --epa-ng-bin."
        )
    result = subprocess.run(
        [path, "--version"], capture_output=True, text=True, check=False
    )
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(f"Could not run EPA-ng at {path}: {detail}")
    return path


def build_epa_ng_command(
    reference_alignment: str | Path,
    reference_tree: str | Path,
    query_alignment: str | Path,
    model_path: str | Path,
    outdir: str | Path,
    threads: int,
    executable: str = "epa-ng",
) -> list[str]:
    return [
        executable,
        "--ref-msa",
        str(reference_alignment),
        "--tree",
        str(reference_tree),
        "--query",
        str(query_alignment),
        "--model",
        str(model_path),
        "--outdir",
        str(outdir),
        "--threads",
        str(threads),
    ]


def _run_epa_ng(command: list[str], outdir: Path, stage: str) -> Path:
    if outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True)
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(f"EPA-ng {stage} failed (exit {result.returncode}): {detail}")
    jplace = outdir / "epa_result.jplace"
    if not jplace.exists():
        raise RuntimeError(f"EPA-ng {stage} completed without producing {jplace}")
    return jplace


@dataclass(frozen=True)
class _PlacementJob:
    locus_id: str
    query_name: str
    query_fasta: Path
    reference_alignment: Path
    reference_tree: Path
    reference_model: Path
    placed_alignment: Path
    query_alignment: Path
    epa_outdir: Path


def _run_placement_job(
    job: _PlacementJob,
    mafft: str,
    epa_ng: str,
    native_threads: int = 1,
) -> Path:
    """Add and place one query; jobs are parallelized across independent loci."""
    _run_mafft(
        build_mafft_add_command(
            job.query_fasta, job.reference_alignment, native_threads, mafft
        ),
        job.placed_alignment,
        f"--add placement for {job.locus_id}",
    )
    aligned_query_records = [
        (name, sequence)
        for name, sequence in _read_fasta(job.placed_alignment)
        if name == job.query_name
    ]
    if len(aligned_query_records) != 1:
        raise RuntimeError(
            f"MAFFT placement alignment for {job.locus_id} did not contain exactly "
            f"one {job.query_name!r}"
        )
    _write_fasta(aligned_query_records, job.query_alignment)
    return _run_epa_ng(
        build_epa_ng_command(
            job.reference_alignment,
            job.reference_tree,
            job.query_alignment,
            job.reference_model,
            job.epa_outdir,
            native_threads,
            epa_ng,
        ),
        job.epa_outdir,
        f"placement for {job.locus_id}",
    )


def _read_database_fasta(path: Path, locus_ids: set[str]) -> dict[str, list[tuple[str, str]]]:
    records = _read_fasta(path)
    if path.parent != path and path.stem in locus_ids:
        return {path.stem: records}
    by_locus: dict[str, list[tuple[str, str]]] = {}
    for header, sequence in records:
        parts = header.split("|")
        matching = [part for part in parts if part in locus_ids]
        if len(matching) != 1:
            raise ValueError(
                f"Could not identify one panel locus in FASTA header {header!r}. "
                "Use LOCUS.fasta files or headers such as reference_id|LOCUS."
            )
        locus_id = matching[0]
        reference_id = next((part for part in parts if part != locus_id), "")
        if not reference_id:
            raise ValueError(f"No reference id was present in FASTA header {header!r}")
        by_locus.setdefault(locus_id, []).append((reference_id, sequence))
    return by_locus


def read_sequence_database(
    database_path: str | Path, locus_ids: set[str]
) -> dict[str, list[tuple[str, str]]]:
    """Read per-locus references from a directory, FASTA, or long-form TSV."""
    path = Path(database_path)
    if not path.exists():
        raise ValueError(f"Sequence database path does not exist: {path}")
    by_locus: dict[str, list[tuple[str, str]]] = {}
    if path.is_dir():
        fasta_paths = sorted(
            item for item in path.iterdir() if item.is_file() and item.suffix.lower() in _FASTA_SUFFIXES
        )
        if not fasta_paths:
            raise ValueError(f"Sequence database directory contains no FASTA files: {path}")
        for fasta_path in fasta_paths:
            if fasta_path.stem not in locus_ids:
                continue
            by_locus[fasta_path.stem] = _read_fasta(fasta_path)
    elif path.suffix.lower() in _FASTA_SUFFIXES:
        by_locus = _read_database_fasta(path, locus_ids)
    else:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"reference_id", "locus_id", "sequence"}
            if not reader.fieldnames or not required.issubset(reader.fieldnames):
                raise ValueError(
                    "Sequence database TSV requires reference_id, locus_id, and sequence columns"
                )
            for row in reader:
                locus_id = row["locus_id"]
                if locus_id in locus_ids:
                    by_locus.setdefault(locus_id, []).append(
                        (row["reference_id"], row["sequence"].upper())
                    )
    for locus_id, records in by_locus.items():
        names = [name for name, _sequence in records]
        if not records:
            raise ValueError(f"No reference sequences found for locus {locus_id!r}")
        if len(names) != len(set(names)):
            raise ValueError(f"Duplicate reference ids found for locus {locus_id!r}")
        if any(not name or not sequence for name, sequence in records):
            raise ValueError(f"Blank reference id or sequence found for locus {locus_id!r}")
    if not by_locus:
        raise ValueError("Sequence database contains no loci matching the supplied panel")
    return by_locus


def _resolve_reference_database_layout(
    database_path: str | Path,
) -> tuple[Path, Path | None]:
    """Resolve sequence data and optional reusable trees from a reference build."""
    path = Path(database_path)
    bundled_database = path / "database"
    bundled_phylogeny = path / "phylogeny"
    if path.is_dir() and bundled_database.is_dir():
        return bundled_database, bundled_phylogeny if bundled_phylogeny.is_dir() else None
    sibling_phylogeny = path.parent / "phylogeny"
    if path.is_dir() and path.name == "database" and sibling_phylogeny.is_dir():
        return path, sibling_phylogeny
    return path, None


def _copy_reference_artifact(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() != destination.resolve():
        shutil.copyfile(source, destination)


def read_reference_metadata(path: str | Path | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    metadata_path = Path(path)
    if not metadata_path.exists():
        raise ValueError(f"Reference metadata path does not exist: {metadata_path}")
    with metadata_path.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        first_line = sample.splitlines()[0] if sample.splitlines() else ""
        delimiter = "\t" if "\t" in first_line else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames or "reference_id" not in reader.fieldnames:
            raise ValueError("Reference metadata requires a reference_id column")
        rows = list(reader)
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        reference_id = str(row.get("reference_id", "")).strip()
        if not reference_id:
            continue
        if reference_id in result:
            raise ValueError(f"Duplicate reference metadata for {reference_id!r}")
        result[reference_id] = {str(key): str(value or "") for key, value in row.items()}
    return result


def _parse_newick(text: str) -> _Node:
    """Parse the branch lengths and leaf names needed from a RAxML Newick tree."""
    text = text.strip()
    position = 0

    def skip_space() -> None:
        nonlocal position
        while position < len(text) and text[position].isspace():
            position += 1

    def label() -> str:
        nonlocal position
        skip_space()
        if position < len(text) and text[position] in {"'", '"'}:
            quote = text[position]
            position += 1
            parts = []
            while position < len(text):
                if text[position] == quote:
                    if position + 1 < len(text) and text[position + 1] == quote:
                        parts.append(quote)
                        position += 2
                        continue
                    position += 1
                    return "".join(parts)
                parts.append(text[position])
                position += 1
            raise ValueError("Unterminated quoted name in RAxML-NG Newick tree")
        start = position
        while position < len(text) and text[position] not in ",():;":
            position += 1
        return text[start:position].strip()

    def branch_length() -> tuple[float, int | None]:
        nonlocal position
        skip_space()
        if position >= len(text) or text[position] != ":":
            return 0.0, None
        position += 1
        start = position
        while position < len(text) and text[position] not in ",);{":
            position += 1
        value = text[start:position].strip()
        try:
            length = float(value)
        except ValueError as exc:
            raise ValueError(f"Invalid branch length {value!r} in RAxML-NG tree") from exc
        edge_num = None
        if position < len(text) and text[position] == "{":
            position += 1
            edge_start = position
            while position < len(text) and text[position] != "}":
                position += 1
            if position >= len(text):
                raise ValueError("Unterminated edge number in jplace tree")
            try:
                edge_num = int(text[edge_start:position])
            except ValueError as exc:
                raise ValueError("Invalid edge number in jplace tree") from exc
            position += 1
        return length, edge_num

    def subtree() -> tuple[_Node, float]:
        nonlocal position
        skip_space()
        if position < len(text) and text[position] == "(":
            position += 1
            children: list[tuple[_Node, float]] = []
            while True:
                children.append(subtree())
                skip_space()
                if position < len(text) and text[position] == ",":
                    position += 1
                    continue
                if position >= len(text) or text[position] != ")":
                    raise ValueError("Malformed RAxML-NG Newick tree")
                position += 1
                break
            # RAxML may put a support value or internal-node label here. It is
            # irrelevant for patristic distance but must still be consumed.
            internal_name = label()
            node = _Node(internal_name or None, children)
            length, node.edge_num = branch_length()
            return node, length
        name = label()
        if not name:
            raise ValueError("Blank leaf name in RAxML-NG Newick tree")
        node = _Node(name, [])
        length, node.edge_num = branch_length()
        return node, length

    root, _root_length = subtree()
    skip_space()
    if position < len(text) and text[position] == ";":
        position += 1
    skip_space()
    if position != len(text):
        raise ValueError("Unexpected trailing content in RAxML-NG Newick tree")
    return root


def _placement_patristic_distances(
    root: _Node,
    placement_edge: int,
    distal_length: float,
    pendant_length: float,
) -> dict[str, float]:
    adjacency: dict[int, list[tuple[_Node, float]]] = {}
    target: tuple[_Node, _Node, float] | None = None

    def connect(node: _Node) -> None:
        nonlocal target
        adjacency.setdefault(id(node), [])
        for child, length in node.children:
            adjacency[id(node)].append((child, length))
            adjacency.setdefault(id(child), []).append((node, length))
            if child.edge_num == placement_edge:
                target = (node, child, length)
            connect(child)

    connect(root)
    if target is None:
        raise ValueError(f"Placement edge {placement_edge} is absent from the jplace tree")
    parent, child, edge_length = target
    distal = min(max(0.0, distal_length), edge_length)
    distances: dict[str, float] = {}

    def walk(node: _Node, blocked: _Node, total: float) -> None:
        if not node.children and node.name is not None:
            distances[node.name] = total
        for neighbor, length in adjacency[id(node)]:
            if neighbor is not blocked:
                walk(neighbor, node, total + length)

    walk(child, parent, pendant_length + distal)
    walk(parent, child, pendant_length + edge_length - distal)
    return distances


def _reference_tree_scale(root: _Node) -> float:
    """Median non-zero patristic distance among reference tips."""
    adjacency: dict[int, list[tuple[_Node, float]]] = {}
    leaves: list[_Node] = []

    def connect(node: _Node) -> None:
        adjacency.setdefault(id(node), [])
        if not node.children and node.name is not None:
            leaves.append(node)
        for child, length in node.children:
            adjacency[id(node)].append((child, length))
            adjacency.setdefault(id(child), []).append((node, length))
            connect(child)

    connect(root)

    def distance_between(source: _Node, target: _Node) -> float:
        stack = [(source, None, 0.0)]
        while stack:
            node, parent, distance = stack.pop()
            if node is target:
                return distance
            for neighbor, length in adjacency[id(node)]:
                if neighbor is not parent:
                    stack.append((neighbor, node, distance + length))
        raise ValueError("Reference tree is disconnected")

    distances = [
        distance_between(leaves[left], leaves[right])
        for left in range(len(leaves))
        for right in range(left + 1, len(leaves))
    ]
    positive = [distance for distance in distances if distance > 0]
    return statistics.median(positive) if positive else 1.0


def _tip_patristic_distances(root: _Node) -> dict[tuple[str, str], float]:
    adjacency: dict[int, list[tuple[_Node, float]]] = {}
    leaves: list[_Node] = []

    def connect(node: _Node) -> None:
        adjacency.setdefault(id(node), [])
        if not node.children and node.name is not None:
            leaves.append(node)
        for child, length in node.children:
            adjacency[id(node)].append((child, length))
            adjacency.setdefault(id(child), []).append((node, length))
            connect(child)

    connect(root)

    def walk(source: _Node, target: _Node) -> float:
        stack = [(source, None, 0.0)]
        while stack:
            node, parent, distance = stack.pop()
            if node is target:
                return distance
            for neighbor, length in adjacency[id(node)]:
                if neighbor is not parent:
                    stack.append((neighbor, node, distance + length))
        raise ValueError("Reference tree is disconnected")

    return {
        tuple(sorted((str(left.name), str(right.name)))): walk(left, right)
        for left_index, left in enumerate(leaves)
        for right in leaves[left_index + 1 :]
    }


def _newick_label(label: str) -> str:
    return "'" + str(label).replace("'", "_") + "'"


def neighbor_joining_tree(
    labels: list[str], distances: dict[tuple[str, str], float]
) -> str:
    """Build a deterministic Newick neighbor-joining tree from a distance matrix."""
    if not labels:
        return ";\n"
    if len(labels) == 1:
        return f"({_newick_label(labels[0])}:0.00000000);\n"
    active = list(labels)
    clusters = {label: _newick_label(label) for label in labels}
    matrix: dict[tuple[str, str], float] = {}

    def key(left: str, right: str) -> tuple[str, str]:
        return tuple(sorted((left, right)))

    for left_index, left in enumerate(labels):
        for right in labels[left_index + 1 :]:
            pair = key(left, right)
            if pair not in distances:
                raise ValueError(f"Missing distance between {left!r} and {right!r}")
            matrix[pair] = max(0.0, float(distances[pair]))

    node_number = 0
    while len(active) > 2:
        size = len(active)
        row_sums = {
            label: sum(
                matrix[key(label, other)] for other in active if other != label
            )
            for label in active
        }
        left, right = min(
            (
                (left, right)
                for left_index, left in enumerate(active)
                for right in active[left_index + 1 :]
            ),
            key=lambda pair: (
                (size - 2) * matrix[key(*pair)]
                - row_sums[pair[0]]
                - row_sums[pair[1]],
                pair,
            ),
        )
        pair_distance = matrix[key(left, right)]
        left_length = 0.5 * pair_distance + (
            row_sums[left] - row_sums[right]
        ) / (2 * (size - 2))
        right_length = pair_distance - left_length
        left_length = max(0.0, left_length)
        right_length = max(0.0, right_length)
        node_number += 1
        joined = f"__NJ_{node_number}"
        clusters[joined] = (
            f"({clusters[left]}:{left_length:.8f},"
            f"{clusters[right]}:{right_length:.8f})"
        )
        for other in active:
            if other in {left, right}:
                continue
            matrix[key(joined, other)] = max(
                0.0,
                (
                    matrix[key(left, other)]
                    + matrix[key(right, other)]
                    - pair_distance
                )
                / 2,
            )
        active = [label for label in active if label not in {left, right}]
        active.append(joined)

    left, right = active
    final_distance = matrix[key(left, right)] / 2
    return (
        f"({clusters[left]}:{max(0.0, final_distance):.8f},"
        f"{clusters[right]}:{max(0.0, final_distance):.8f});\n"
    )


def read_epa_ng_placement(jplace_path: str | Path) -> tuple[dict[str, float], dict[str, float | int]]:
    best_distances, _expected_distances, metadata = read_epa_ng_placement_statistics(
        jplace_path
    )
    return best_distances, metadata


def read_epa_ng_placement_statistics(
    jplace_path: str | Path,
) -> tuple[dict[str, float], dict[str, float], dict[str, float | int]]:
    """Return best-placement and LWR-weighted expected tip distances."""
    with Path(jplace_path).open() as handle:
        document = json.load(handle)
    fields = list(document.get("fields", []))
    placements = document.get("placements", [])
    if not placements:
        raise ValueError(f"EPA-ng produced no placements in {jplace_path}")
    required = {"edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length"}
    if not required.issubset(fields):
        raise ValueError(f"EPA-ng jplace is missing required fields: {sorted(required - set(fields))}")
    indexes = {field: fields.index(field) for field in required}
    candidates = placements[0].get("p", [])
    if not candidates:
        raise ValueError(f"EPA-ng produced an empty placement record in {jplace_path}")
    raw_weights = [max(0.0, float(row[indexes["like_weight_ratio"]])) for row in candidates]
    weight_total = sum(raw_weights)
    if weight_total <= 0:
        weights = [1.0 / len(candidates)] * len(candidates)
    else:
        weights = [weight / weight_total for weight in raw_weights]
    best_index = max(range(len(candidates)), key=lambda index: weights[index])
    best = candidates[best_index]
    edge_num = int(best[indexes["edge_num"]])
    distal_length = float(best[indexes["distal_length"]])
    pendant_length = float(best[indexes["pendant_length"]])
    root = _parse_newick(str(document["tree"]))
    best_distances = _placement_patristic_distances(
        root, edge_num, distal_length, pendant_length
    )
    expected_distances = {reference_id: 0.0 for reference_id in best_distances}
    for row, weight in zip(candidates, weights):
        candidate_distances = _placement_patristic_distances(
            root,
            int(row[indexes["edge_num"]]),
            float(row[indexes["distal_length"]]),
            float(row[indexes["pendant_length"]]),
        )
        if set(candidate_distances) != set(expected_distances):
            raise ValueError("EPA-ng placement candidates refer to inconsistent trees")
        for reference_id, distance in candidate_distances.items():
            expected_distances[reference_id] += weight * distance
    entropy = -sum(weight * math.log(weight) for weight in weights if weight > 0)
    metadata: dict[str, float | int] = {
        "placement_edge": edge_num,
        "like_weight_ratio": float(best[indexes["like_weight_ratio"]]),
        "pendant_length": pendant_length,
        "distal_length": distal_length,
        "placement_count": len(candidates),
        "placement_entropy": entropy,
        "reference_tree_scale": _reference_tree_scale(root),
    }
    return best_distances, expected_distances, metadata


def run_phylogenetic_placement(
    query_sequences: dict[str, str],
    database_path: str | Path,
    outdir: str | Path,
    sample_id: str,
    locus_ids: set[str] | list[Locus],
    threads: int,
    mafft_bin: str = "mafft",
    raxml_ng_bin: str = "raxml-ng",
    epa_ng_bin: str = "epa-ng",
    raxml_model: str = "DNA",
    snp_weight: float = 1.0,
    repeat_weight: float = 1.0,
    reference_metadata_path: str | Path | None = None,
    progress: ProgressReporter | None = None,
) -> dict[str, Path]:
    if snp_weight < 0 or repeat_weight < 0 or snp_weight + repeat_weight <= 0:
        raise ValueError("SNP and repeat weights must be non-negative with a positive total")
    locus_by_id = (
        {locus.locus_id: locus for locus in locus_ids}
        if not isinstance(locus_ids, set)
        else {}
    )
    requested_locus_ids = set(locus_by_id) if locus_by_id else set(locus_ids)
    sequence_database_path, reusable_phylogeny = _resolve_reference_database_layout(
        database_path
    )
    references = read_sequence_database(sequence_database_path, requested_locus_ids)
    if reference_metadata_path is None and sequence_database_path.is_dir():
        automatic_metadata = sequence_database_path / "reference_metadata.tsv"
        if automatic_metadata.exists():
            reference_metadata_path = automatic_metadata
    reference_metadata = read_reference_metadata(reference_metadata_path)
    mafft = check_mafft(mafft_bin)
    epa_ng = check_epa_ng(epa_ng_bin)
    raxml_ng: str | None = None
    output = Path(outdir) / "phylogeny"
    output.mkdir(parents=True, exist_ok=True)
    detail_rows: list[dict] = []
    status_rows: list[dict] = []
    safe_sample = _SAFE_FILE.sub("_", sample_id).strip("_") or "sample"
    query_name = f"QUERY__{safe_sample}"
    placed_loci: list[str] = []
    marker_rows: list[dict] = []
    query_components_by_locus: dict[str, MarkerComponents] = {}
    reference_components_by_locus: dict[str, dict[str, MarkerComponents]] = {}
    reference_pairwise_by_locus: dict[str, dict[tuple[str, str], float]] = {}
    placement_jobs: dict[str, _PlacementJob] = {}

    for locus_id in sorted(references):
        query_sequence = query_sequences.get(locus_id, "")
        locus = locus_by_id.get(locus_id)
        locus_references = references[locus_id]
        if locus is not None:
            reference_components = {
                reference_id: decompose_marker_sequence(locus, sequence)
                for reference_id, sequence in locus_references
            }
            reference_components_by_locus[locus_id] = reference_components
            for reference_id, components in reference_components.items():
                marker_rows.append(
                    _marker_component_row(
                        sample_id, locus, "reference", reference_id, components
                    )
                )
            if query_sequence:
                query_components = decompose_marker_sequence(locus, query_sequence)
                query_components_by_locus[locus_id] = query_components
                marker_rows.append(
                    _marker_component_row(
                        sample_id, locus, "query", query_name, query_components
                    )
                )
                query_sequence = query_components.snp_sequence
            locus_references = [
                (reference_id, reference_components[reference_id].snp_sequence)
                for reference_id, _sequence in locus_references
            ]
        if query_sequence and any(reference_id == query_name for reference_id, _sequence in references[locus_id]):
            raise ValueError(
                f"Reference id {query_name!r} is reserved for the query sequence"
            )
        safe_locus = _SAFE_FILE.sub("_", locus_id).strip("_") or "locus"
        locus_dir = output / safe_locus
        locus_dir.mkdir(parents=True, exist_ok=True)
        reference_fasta = locus_dir / "references.fasta"
        reference_alignment = locus_dir / "references.aligned.fasta"
        reference_tree = locus_dir / "reference_tree.nwk"
        reference_prefix = locus_dir / "reference"
        reference_model = Path(f"{reference_prefix}.raxml.bestModel")
        query_fasta = locus_dir / "query.fasta"
        placed_alignment = locus_dir / "query_placed.aligned.fasta"
        query_alignment = locus_dir / "query.aligned.fasta"
        epa_outdir = locus_dir / "epa-ng"
        _write_fasta(locus_references, reference_fasta)
        if reusable_phylogeny is not None:
            reusable_locus_dir = reusable_phylogeny / safe_locus
            reusable_alignment = reusable_locus_dir / "references.aligned.fasta"
            reusable_tree = reusable_locus_dir / "reference_tree.nwk"
            reusable_model = reusable_locus_dir / "reference.raxml.bestModel"
            missing_artifacts = [
                path
                for path in (reusable_alignment, reusable_tree, reusable_model)
                if not path.is_file()
            ]
            if missing_artifacts:
                status_rows.append(
                    {
                        "locus_id": locus_id,
                        "reference_sequences": len(references[locus_id]),
                        "query_sequence": "yes" if query_sequence else "no",
                        "status": "REFERENCE_TREE_UNAVAILABLE",
                    }
                )
                continue
            _copy_reference_artifact(reusable_alignment, reference_alignment)
            _copy_reference_artifact(reusable_tree, reference_tree)
            _copy_reference_artifact(reusable_model, reference_model)
        else:
            _run_mafft(
                build_mafft_reference_command(reference_fasta, threads, mafft),
                reference_alignment,
                f"reference alignment for {locus_id}",
            )
            if raxml_ng is None:
                raxml_ng = check_raxml_ng(raxml_ng_bin)
            _run_raxml_ng(
                build_raxml_ng_command(
                    reference_alignment,
                    reference_prefix,
                    threads,
                    raxml_ng,
                    raxml_model,
                ),
                reference_prefix,
                reference_tree,
                f"reference tree search for {locus_id}",
            )
        if not reference_model.exists():
            raise RuntimeError(
                f"RAxML-NG reference search did not produce model file {reference_model}"
            )
        if not query_sequence:
            status_rows.append(
                {
                    "locus_id": locus_id,
                    "reference_sequences": len(references[locus_id]),
                    "query_sequence": "no",
                    "status": "REFERENCE_TREE_ONLY",
                }
            )
            continue
        _write_fasta([(query_name, query_sequence)], query_fasta)
        placement_jobs[locus_id] = _PlacementJob(
            locus_id=locus_id,
            query_name=query_name,
            query_fasta=query_fasta,
            reference_alignment=reference_alignment,
            reference_tree=reference_tree,
            reference_model=reference_model,
            placed_alignment=placed_alignment,
            query_alignment=query_alignment,
            epa_outdir=epa_outdir,
        )

    placement_results: dict[str, Path] = {}
    if placement_jobs:
        cpu_budget = resolve_threads(threads)
        worker_count = min(cpu_budget, len(placement_jobs))
        native_threads = max(1, cpu_budget // worker_count)
        if progress is not None:
            progress.step(
                f"Running {len(placement_jobs):,} independent EPA-ng locus "
                f"placements with {worker_count} worker(s) and "
                f"{native_threads} native thread(s) per worker"
            )
        if worker_count == 1:
            for completed, (locus_id, job) in enumerate(
                placement_jobs.items(), start=1
            ):
                placement_results[locus_id] = _run_placement_job(
                    job, mafft, epa_ng, native_threads
                )
                if progress is not None:
                    progress.count(
                        "Completed EPA-ng loci",
                        completed,
                        len(placement_jobs),
                        force=True,
                    )
        else:
            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                futures = {
                    executor.submit(
                        _run_placement_job, job, mafft, epa_ng, native_threads
                    ): locus_id
                    for locus_id, job in placement_jobs.items()
                }
                for completed, future in enumerate(as_completed(futures), start=1):
                    locus_id = futures[future]
                    try:
                        placement_results[locus_id] = future.result()
                    except Exception as exc:
                        raise RuntimeError(
                            f"Phylogenetic placement failed for locus {locus_id!r}"
                        ) from exc
                    if progress is not None:
                        progress.count(
                            "Completed EPA-ng loci",
                            completed,
                            len(placement_jobs),
                            force=True,
                        )

    # Parse completed placements in locus order for deterministic tables.
    for locus_id in sorted(placement_results):
        jplace_path = placement_results[locus_id]
        (
            placement_distances,
            expected_placement_distances,
            placement_metadata,
        ) = read_epa_ng_placement_statistics(jplace_path)
        with Path(jplace_path).open() as jplace_handle:
            jplace_tree = _parse_newick(str(json.load(jplace_handle)["tree"]))
        reference_pairwise_by_locus[locus_id] = _tip_patristic_distances(
            jplace_tree
        )
        expected_references = {reference_id for reference_id, _sequence in references[locus_id]}
        if set(placement_distances) != expected_references:
            missing = sorted(expected_references - set(placement_distances))
            unexpected = sorted(set(placement_distances) - expected_references)
            raise RuntimeError(
                f"EPA-ng tree/reference mismatch for {locus_id}: missing={missing}, unexpected={unexpected}"
            )
        for reference_id, distance in placement_distances.items():
            detail_rows.append(
                {
                    "sample_id": sample_id,
                    "locus_id": locus_id,
                    "reference_id": reference_id,
                    "phylogenetic_distance": f"{distance:.8f}",
                    "likelihood_weighted_phylogenetic_distance": (
                        f"{expected_placement_distances[reference_id]:.8f}"
                    ),
                    "placement_edge": placement_metadata["placement_edge"],
                    "like_weight_ratio": f'{placement_metadata["like_weight_ratio"]:.8f}',
                    "pendant_length": f'{placement_metadata["pendant_length"]:.8f}',
                    "distal_length": f'{placement_metadata["distal_length"]:.8f}',
                    "placement_count": placement_metadata["placement_count"],
                    "placement_entropy": f'{placement_metadata["placement_entropy"]:.8f}',
                    "reference_tree_scale": f'{placement_metadata["reference_tree_scale"]:.8f}',
                }
            )
        placed_loci.append(locus_id)
        status_rows.append(
            {
                "locus_id": locus_id,
                "reference_sequences": len(references[locus_id]),
                "query_sequence": "yes",
                "status": "PLACED",
            }
        )

    detail_path = output / "locus_phylogenetic_distances.tsv"
    _write_tsv(detail_rows, detail_path, PLACEMENT_FIELDS)
    marker_components_path = output / "marker_components.tsv"
    _write_tsv(marker_rows, marker_components_path, MARKER_COMPONENT_FIELDS)

    repeat_scale_by_locus: dict[str, float] = {}
    for locus_id, components_by_reference in reference_components_by_locus.items():
        counts = [
            float(components.repeat_count_raw)
            for components in components_by_reference.values()
            if components.repeat_count_raw is not None
        ]
        repeat_scale_by_locus[locus_id] = max(
            0.5, statistics.pstdev(counts) if len(counts) > 1 else 0.0
        )

    locus_marker_rows: list[dict] = []
    for row in detail_rows:
        locus_id = str(row["locus_id"])
        reference_id = str(row["reference_id"])
        query_components = query_components_by_locus.get(locus_id)
        reference_components = reference_components_by_locus.get(locus_id, {}).get(
            reference_id
        )
        query_count = (
            query_components.repeat_count_raw if query_components is not None else None
        )
        reference_count = (
            reference_components.repeat_count_raw
            if reference_components is not None
            else None
        )
        if query_count is not None and reference_count is not None:
            repeat_delta: str | float = abs(query_count - reference_count)
            normalized_repeat: str | float = (
                repeat_delta / repeat_scale_by_locus[locus_id]
            )
        else:
            repeat_delta = ""
            normalized_repeat = ""
        weighted_snp = float(row["likelihood_weighted_phylogenetic_distance"])
        tree_scale = max(float(row["reference_tree_scale"]), 1e-12)
        locus_marker_rows.append(
            {
                "sample_id": sample_id,
                "locus_id": locus_id,
                "reference_id": reference_id,
                "query_repeat_count": "" if query_count is None else f"{query_count:.6f}",
                "reference_repeat_count": ""
                if reference_count is None
                else f"{reference_count:.6f}",
                "repeat_count_delta": ""
                if repeat_delta == ""
                else f"{float(repeat_delta):.6f}",
                "normalized_repeat_distance": ""
                if normalized_repeat == ""
                else f"{float(normalized_repeat):.8f}",
                "likelihood_weighted_snp_distance": f"{weighted_snp:.8f}",
                "reference_tree_scale": f"{tree_scale:.8f}",
                "normalized_snp_distance": f"{weighted_snp / tree_scale:.8f}",
                "placement_entropy": row["placement_entropy"],
            }
        )
    locus_marker_path = output / "locus_marker_distances.tsv"
    _write_tsv(
        locus_marker_rows, locus_marker_path, LOCUS_MARKER_DISTANCE_FIELDS
    )
    totals: dict[str, list[float]] = {}
    weighted_totals: dict[str, list[float]] = {}
    for row in detail_rows:
        reference_id = str(row["reference_id"])
        totals.setdefault(reference_id, []).append(float(row["phylogenetic_distance"]))
        weighted_totals.setdefault(reference_id, []).append(
            float(row["likelihood_weighted_phylogenetic_distance"])
        )
    # A partial reference must not win merely because absent loci contribute
    # nothing. Rank only references represented in every successfully placed
    # locus, so the total is comparable across candidates.
    complete_totals = {
        reference_id: values
        for reference_id, values in totals.items()
        if len(values) == len(placed_loci)
    }
    summary_rows = []
    ordered = sorted(
        complete_totals.items(),
        key=lambda item: (
            sum(weighted_totals[item[0]]),
            sum(item[1]),
            item[0],
        ),
    )
    for rank, (reference_id, values) in enumerate(ordered, start=1):
        weighted_values = weighted_totals[reference_id]
        weighted_total = sum(weighted_values)
        if rank < len(ordered):
            next_reference = ordered[rank][0]
            next_total = sum(weighted_totals[next_reference])
            distance_gap: str | float = f"{next_total - weighted_total:.8f}"
            relative_gap: str | float = f"{(next_total - weighted_total) / max(next_total, 1e-12):.8f}"
        else:
            distance_gap = ""
            relative_gap = ""
        summary_rows.append(
            {
                "sample_id": sample_id,
                "reference_id": reference_id,
                "total_phylogenetic_distance": f"{sum(values):.8f}",
                "total_likelihood_weighted_distance": f"{weighted_total:.8f}",
                "compared_loci": len(values),
                "mean_phylogenetic_distance": f"{sum(values) / len(values):.8f}",
                "mean_likelihood_weighted_distance": (
                    f"{weighted_total / len(weighted_values):.8f}"
                ),
                "distance_gap_to_next": distance_gap,
                "relative_distance_gap_to_next": relative_gap,
                "rank": rank,
            }
        )
    summary_path = output / "phylogenetic_matches.tsv"
    _write_tsv(summary_rows, summary_path, SUMMARY_FIELDS)

    complete_references = set(complete_totals)
    combined_totals: dict[str, dict[str, float | int]] = {}
    for row in locus_marker_rows:
        reference_id = str(row["reference_id"])
        if reference_id not in complete_references:
            continue
        values = combined_totals.setdefault(
            reference_id,
            {
                "snp": 0.0,
                "normalized_snp": 0.0,
                "repeat": 0.0,
                "normalized_repeat": 0.0,
                "loci": 0,
                "repeat_loci": 0,
            },
        )
        values["snp"] += float(row["likelihood_weighted_snp_distance"])
        values["normalized_snp"] += float(row["normalized_snp_distance"])
        values["loci"] += 1
        if row["repeat_count_delta"] != "":
            values["repeat"] += float(row["repeat_count_delta"])
            values["normalized_repeat"] += float(row["normalized_repeat_distance"])
            values["repeat_loci"] += 1

    required_repeat_loci = {
        locus_id
        for locus_id in placed_loci
        if query_components_by_locus.get(locus_id) is not None
        and query_components_by_locus[locus_id].repeat_count_raw is not None
    }
    eligible_combined_references = [
        reference_id
        for reference_id, values in combined_totals.items()
        if int(values["repeat_loci"]) == len(required_repeat_loci)
    ]
    combined_order = sorted(
        eligible_combined_references,
        key=lambda reference_id: (
            snp_weight * float(combined_totals[reference_id]["normalized_snp"])
            + repeat_weight
            * float(combined_totals[reference_id]["normalized_repeat"]),
            reference_id,
        ),
    )
    combined_rows: list[dict] = []
    combined_scores = {
        reference_id: (
            snp_weight * float(combined_totals[reference_id]["normalized_snp"])
            + repeat_weight
            * float(combined_totals[reference_id]["normalized_repeat"])
        )
        for reference_id in combined_order
    }
    for index, reference_id in enumerate(combined_order):
        values = combined_totals[reference_id]
        score = combined_scores[reference_id]
        if index + 1 < len(combined_order):
            next_score = combined_scores[combined_order[index + 1]]
            gap = next_score - score
            relative_gap: str | float = gap / max(next_score, 1e-12)
        else:
            gap = ""
            relative_gap = ""
        combined_rows.append(
            {
                "sample_id": sample_id,
                "reference_id": reference_id,
                "total_likelihood_weighted_snp_distance": f'{float(values["snp"]):.8f}',
                "total_normalized_snp_distance": f'{float(values["normalized_snp"]):.8f}',
                "total_repeat_count_distance": f'{float(values["repeat"]):.8f}',
                "total_normalized_repeat_distance": f'{float(values["normalized_repeat"]):.8f}',
                "snp_weight": f"{snp_weight:.6f}",
                "repeat_weight": f"{repeat_weight:.6f}",
                "combined_marker_distance": f"{score:.8f}",
                "compared_loci": values["loci"],
                "repeat_compared_loci": values["repeat_loci"],
                "distance_gap_to_next": "" if gap == "" else f"{float(gap):.8f}",
                "relative_distance_gap_to_next": ""
                if relative_gap == ""
                else f"{float(relative_gap):.8f}",
                "collection_date": reference_metadata.get(reference_id, {}).get(
                    "collection_date", ""
                ),
                "latitude": reference_metadata.get(reference_id, {}).get("latitude", ""),
                "longitude": reference_metadata.get(reference_id, {}).get("longitude", ""),
                "location": reference_metadata.get(reference_id, {}).get("location", ""),
                "source": reference_metadata.get(reference_id, {}).get("source", ""),
                "rank": index + 1,
            }
        )
    combined_path = output / "combined_marker_matches.tsv"
    _write_tsv(combined_rows, combined_path, COMBINED_MARKER_FIELDS)
    query_tree_label = sample_id if sample_id not in eligible_combined_references else query_name
    tree_labels = [*sorted(eligible_combined_references), query_tree_label]
    query_locus_rows = {
        (str(row["locus_id"]), str(row["reference_id"])): row
        for row in locus_marker_rows
    }
    tree_scale_by_locus = {
        str(row["locus_id"]): max(float(row["reference_tree_scale"]), 1e-12)
        for row in detail_rows
    }

    def components_for(label: str, locus_id: str) -> MarkerComponents | None:
        if label == query_tree_label:
            return query_components_by_locus.get(locus_id)
        return reference_components_by_locus.get(locus_id, {}).get(label)

    combined_pairwise: dict[tuple[str, str], float] = {}
    for left_index, left in enumerate(tree_labels):
        for right in tree_labels[left_index + 1 :]:
            snp_total = 0.0
            repeat_total = 0.0
            for locus_id in placed_loci:
                if query_tree_label in {left, right}:
                    reference_id = right if left == query_tree_label else left
                    marker_row = query_locus_rows[(locus_id, reference_id)]
                    snp_total += float(marker_row["normalized_snp_distance"])
                else:
                    pair = tuple(sorted((left, right)))
                    snp_total += (
                        reference_pairwise_by_locus[locus_id][pair]
                        / tree_scale_by_locus[locus_id]
                    )
                left_components = components_for(left, locus_id)
                right_components = components_for(right, locus_id)
                if (
                    left_components is not None
                    and right_components is not None
                    and left_components.repeat_count_raw is not None
                    and right_components.repeat_count_raw is not None
                ):
                    repeat_total += abs(
                        left_components.repeat_count_raw
                        - right_components.repeat_count_raw
                    ) / repeat_scale_by_locus[locus_id]
            combined_pairwise[tuple(sorted((left, right)))] = (
                snp_weight * snp_total + repeat_weight * repeat_total
            )
    combined_tree_path = output / "combined_markers.tree"
    combined_tree_path.write_text(
        neighbor_joining_tree(tree_labels, combined_pairwise)
    )
    status_path = output / "locus_status.tsv"
    _write_tsv(status_rows, status_path, STATUS_FIELDS)
    return {
        "phylogeny": output,
        "phylogenetic_distances": detail_path,
        "phylogenetic_matches": summary_path,
        "phylogenetic_status": status_path,
        "marker_components": marker_components_path,
        "locus_marker_distances": locus_marker_path,
        "combined_marker_matches": combined_path,
        "combined_marker_tree": combined_tree_path,
    }


def build_reference_phylogenies(
    database_path: str | Path,
    outdir: str | Path,
    loci: list[Locus],
    threads: int,
    *,
    min_references: int = 3,
    mafft_bin: str = "mafft",
    raxml_ng_bin: str = "raxml-ng",
    raxml_model: str = "DNA",
    progress: ProgressReporter | None = None,
) -> dict[str, Path]:
    """Build reusable, repeat-masked reference alignments and locus trees.

    The raw amplicons remain in the database.  Trees use the SNP component
    when a repeat can be bounded from panel flanks or motif; this prevents
    variable repeat length from being interpreted as many independent SNPs.
    """
    if min_references < 2:
        raise ValueError("min_references must be at least 2")
    locus_by_id = {locus.locus_id: locus for locus in loci}
    references = read_sequence_database(database_path, set(locus_by_id))
    mafft = check_mafft(mafft_bin)
    raxml_ng = check_raxml_ng(raxml_ng_bin)
    output = Path(outdir)
    output.mkdir(parents=True, exist_ok=True)
    status_rows: list[dict] = []
    component_rows: list[dict] = []

    sorted_locus_ids = sorted(locus_by_id)
    for locus_number, locus_id in enumerate(sorted_locus_ids, start=1):
        records = references.get(locus_id, [])
        safe_locus = _SAFE_FILE.sub("_", locus_id).strip("_") or "locus"
        locus_dir = output / safe_locus
        portable_tree = output / f"{safe_locus}.tree"
        portable_tree.unlink(missing_ok=True)
        if len(records) < min_references:
            status_rows.append(
                {
                    "locus_id": locus_id,
                    "reference_sequences": len(records),
                    "status": "INSUFFICIENT_REFERENCES",
                    "tree": "",
                }
            )
            if progress is not None:
                progress.count(
                    "Processed tree loci", locus_number, len(sorted_locus_ids), force=True
                )
            continue
        if progress is not None:
            progress.step(
                f"Building tree for {locus_id} ({locus_number}/{len(sorted_locus_ids)}; "
                f"{len(records):,} references)"
            )
        locus_dir.mkdir(parents=True, exist_ok=True)
        masked_records: list[tuple[str, str]] = []
        for reference_id, sequence in records:
            components = decompose_marker_sequence(locus_by_id[locus_id], sequence)
            masked_records.append((reference_id, components.snp_sequence))
            component_rows.append(
                _marker_component_row(
                    "REFERENCE_DATABASE",
                    locus_by_id[locus_id],
                    "reference",
                    reference_id,
                    components,
                )
            )
        reference_fasta = locus_dir / "references.fasta"
        alignment = locus_dir / "references.aligned.fasta"
        tree = locus_dir / "reference_tree.nwk"
        prefix = locus_dir / "reference"
        _write_fasta(masked_records, reference_fasta)
        _run_mafft(
            build_mafft_reference_command(reference_fasta, threads, mafft),
            alignment,
            f"reference alignment for {locus_id}",
        )
        _run_raxml_ng(
            build_raxml_ng_command(alignment, prefix, threads, raxml_ng, raxml_model),
            prefix,
            tree,
            f"reference tree search for {locus_id}",
            progress,
        )
        shutil.copyfile(tree, portable_tree)
        status_rows.append(
            {
                "locus_id": locus_id,
                "reference_sequences": len(records),
                "status": "BUILT",
                "tree": str(portable_tree),
            }
        )
        if progress is not None:
            progress.count(
                "Processed tree loci", locus_number, len(sorted_locus_ids), force=True
            )

    status_path = output / "reference_tree_status.tsv"
    components_path = output / "reference_marker_components.tsv"
    _write_tsv(
        status_rows,
        status_path,
        ["locus_id", "reference_sequences", "status", "tree"],
    )
    _write_tsv(component_rows, components_path, MARKER_COMPONENT_FIELDS)
    return {
        "phylogeny": output,
        "tree_status": status_path,
        "marker_components": components_path,
    }


def dominant_read_query_sequences(
    features: list[RepeatFeature], asv_rows: list[dict]
) -> dict[str, str]:
    """Return the dominant full amplicon observed at each FASTQ locus."""
    feature_by_key = {(feature.locus_id, feature.read_id): feature for feature in features}
    dominant: dict[str, dict] = {}
    for row in asv_rows:
        locus_id = str(row.get("locus_id", ""))
        current = dominant.get(locus_id)
        if current is None or int(row.get("support_reads") or 0) > int(
            current.get("support_reads") or 0
        ):
            dominant[locus_id] = row
    sequences = {}
    for locus_id, row in dominant.items():
        feature = feature_by_key.get((locus_id, str(row.get("representative_read_id") or "")))
        if feature is not None:
            sequence = (feature.amplicon_sequence or feature.repeat_sequence).upper()
            if sequence:
                sequences[locus_id] = sequence
    return sequences
