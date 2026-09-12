from __future__ import annotations

import gzip
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .calling import estimate_repeat_count_from_product_length, normalize_allele, repeat_unit_length
from .models import Locus

RAXML_NG_THREADS_PER_PROCESS = 1

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

def _aligned_snp_distance(query: str, reference: str) -> float:
    """Return an ambiguity-aware mismatch fraction from a shared alignment."""
    if len(query) != len(reference):
        raise ValueError("Aligned query and reference sequences have different lengths")
    mismatches = 0
    comparable = 0
    for query_base, reference_base in zip(query.upper(), reference.upper()):
        if query_base == "-" and reference_base == "-":
            continue
        comparable += 1
        if query_base == reference_base:
            continue
        if query_base == "-" or reference_base == "-":
            mismatches += 1
            continue
        query_states = _IUPAC.get(query_base)
        reference_states = _IUPAC.get(reference_base)
        if query_states is None or reference_states is None or query_states.isdisjoint(
            reference_states
        ):
            mismatches += 1
    return mismatches / comparable if comparable else 0.0

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
        str(RAXML_NG_THREADS_PER_PROCESS),
        "--redo",
    ]

def _run_raxml_ng(
    command: list[str],
    prefix: Path,
    output_tree: Path,
    stage: str,
    progress: ProgressReporter | None = None,
) -> None:
    """Run RAxML-NG with the invariant one-thread process configuration."""
    attempted_command = list(command)
    thread_index = attempted_command.index("--threads") + 1
    attempted_command[thread_index] = str(RAXML_NG_THREADS_PER_PROCESS)
    log_path = Path(f"{prefix}.mlvamaps.raxml.log")
    result = subprocess.run(
        attempted_command, capture_output=True, text=True, check=False
    )
    detail = "\n".join(
        part.strip() for part in (result.stdout, result.stderr) if part.strip()
    )
    rendered_command = " ".join(attempted_command)
    log_path.write_text(f"$ {rendered_command}\n{detail}".rstrip() + "\n")
    if result.returncode:
        detail_tail = "\n".join(detail.splitlines()[-40:])
        raise RuntimeError(
            f"RAxML-NG {stage} failed (exit {result.returncode}). "
            f"Command: {rendered_command}. Full output: {log_path}\n{detail_tail}"
        )

    best_tree = Path(f"{prefix}.raxml.bestTree")
    if not best_tree.exists():
        raise RuntimeError(
            f"RAxML-NG {stage} completed without producing {best_tree}"
        )
    shutil.copyfile(best_tree, output_tree)

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

def _tip_patristic_distance_matrix(root: _Node) -> tuple[list[str], np.ndarray]:
    """Build a dense tip-distance matrix using NumPy's compiled matrix kernels."""
    edge_lengths: list[float] = []
    names: list[str] = []
    tip_paths: list[list[int]] = []

    def collect(node: _Node, path: list[int]) -> None:
        if not node.children and node.name is not None:
            names.append(str(node.name))
            tip_paths.append(path)
        for child, length in node.children:
            edge_index = len(edge_lengths)
            edge_lengths.append(float(length))
            collect(child, [*path, edge_index])

    collect(root, [])
    if len(names) != len(set(names)):
        raise ValueError("Reference tree contains duplicate tip names")
    incidence = np.zeros((len(names), len(edge_lengths)), dtype=np.float64)
    for tip_index, path in enumerate(tip_paths):
        incidence[tip_index, path] = 1.0
    weighted_paths = incidence * np.asarray(edge_lengths, dtype=np.float64)
    root_distances = weighted_paths.sum(axis=1)
    shared_distances = weighted_paths @ incidence.T
    matrix = (
        root_distances[:, None]
        + root_distances[None, :]
        - 2.0 * shared_distances
    )
    np.maximum(matrix, 0.0, out=matrix)
    return names, matrix

def _newick_label(label: str) -> str:
    # Newick escapes a quote inside a quoted label by doubling it. Preserve
    # sample IDs exactly instead of silently rewriting meaningful characters.
    return "'" + str(label).replace("'", "''") + "'"

def neighbor_joining_tree_from_matrix(
    labels: list[str], distances: np.ndarray
) -> str:
    """Build deterministic Newick directly from a finite dense distance matrix."""
    return _neighbor_joining_tree_from_matrix(labels, distances)

def _neighbor_joining_tree_from_matrix(
    labels: list[str], distances: np.ndarray
) -> str:
    """Run neighbor joining on a dense matrix through compiled NumPy kernels."""
    if not labels:
        return ";\n"
    if len(labels) == 1:
        return f"({_newick_label(labels[0])}:0.00000000);\n"
    active = list(labels)
    clusters = {label: _newick_label(label) for label in labels}
    matrix = np.asarray(distances, dtype=np.float64).copy()
    if matrix.shape != (len(labels), len(labels)):
        raise ValueError("Neighbor-joining distance matrix shape does not match labels")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Neighbor-joining distances must be finite")
    np.maximum(matrix, 0.0, out=matrix)
    np.fill_diagonal(matrix, 0.0)

    node_number = 0
    while len(active) > 2:
        size = len(active)
        row_sums = matrix.sum(axis=1)
        q_matrix = (size - 2) * matrix - row_sums[:, None] - row_sums[None, :]
        q_matrix[np.tril_indices(size)] = np.inf
        left_index, right_index = np.unravel_index(
            int(np.argmin(q_matrix)), q_matrix.shape
        )
        left = active[left_index]
        right = active[right_index]
        pair_distance = float(matrix[left_index, right_index])
        left_length = 0.5 * pair_distance + (
            float(row_sums[left_index]) - float(row_sums[right_index])
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
        retained = [
            index for index in range(size) if index not in {left_index, right_index}
        ]
        joined_distances = np.maximum(
            0.0,
            (
                matrix[left_index, retained]
                + matrix[right_index, retained]
                - pair_distance
            )
            / 2,
        )
        next_matrix = np.zeros((size - 1, size - 1), dtype=float)
        next_matrix[:-1, :-1] = matrix[np.ix_(retained, retained)]
        next_matrix[-1, :-1] = joined_distances
        next_matrix[:-1, -1] = joined_distances
        matrix = next_matrix
        active = [active[index] for index in retained]
        active.append(joined)

    left, right = active
    final_distance = float(matrix[0, 1]) / 2
    return (
        f"({clusters[left]}:{max(0.0, final_distance):.8f},"
        f"{clusters[right]}:{max(0.0, final_distance):.8f});\n"
    )
