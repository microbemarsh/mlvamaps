from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .models import RepeatFeature


PLACEMENT_FIELDS = [
    "sample_id",
    "locus_id",
    "reference_id",
    "phylogenetic_distance",
    "placement_edge",
    "like_weight_ratio",
    "pendant_length",
    "distal_length",
]

SUMMARY_FIELDS = [
    "sample_id",
    "reference_id",
    "total_phylogenetic_distance",
    "compared_loci",
    "mean_phylogenetic_distance",
    "rank",
]

STATUS_FIELDS = ["locus_id", "reference_sequences", "query_sequence", "status"]

_FASTA_SUFFIXES = {".fa", ".fas", ".fasta", ".fna", ".ffn"}
_SAFE_FILE = re.compile(r"[^A-Za-z0-9_.-]+")


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
    model: str = "GTR+G",
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


def _run_raxml_ng(command: list[str], prefix: Path, output_tree: Path, stage: str) -> None:
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode:
        detail = (result.stderr or result.stdout or "").strip()
        raise RuntimeError(
            f"RAxML-NG {stage} failed (exit {result.returncode}): {detail}"
        )
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


def read_epa_ng_placement(jplace_path: str | Path) -> tuple[dict[str, float], dict[str, float | int]]:
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
    best = max(candidates, key=lambda row: float(row[indexes["like_weight_ratio"]]))
    edge_num = int(best[indexes["edge_num"]])
    distal_length = float(best[indexes["distal_length"]])
    pendant_length = float(best[indexes["pendant_length"]])
    root = _parse_newick(str(document["tree"]))
    distances = _placement_patristic_distances(root, edge_num, distal_length, pendant_length)
    metadata: dict[str, float | int] = {
        "placement_edge": edge_num,
        "like_weight_ratio": float(best[indexes["like_weight_ratio"]]),
        "pendant_length": pendant_length,
        "distal_length": distal_length,
    }
    return distances, metadata


def run_phylogenetic_placement(
    query_sequences: dict[str, str],
    database_path: str | Path,
    outdir: str | Path,
    sample_id: str,
    locus_ids: set[str],
    threads: int,
    mafft_bin: str = "mafft",
    raxml_ng_bin: str = "raxml-ng",
    epa_ng_bin: str = "epa-ng",
    raxml_model: str = "GTR+G",
) -> dict[str, Path]:
    references = read_sequence_database(database_path, locus_ids)
    mafft = check_mafft(mafft_bin)
    raxml_ng = check_raxml_ng(raxml_ng_bin)
    epa_ng = check_epa_ng(epa_ng_bin)
    output = Path(outdir) / "phylogeny"
    output.mkdir(parents=True, exist_ok=True)
    detail_rows: list[dict] = []
    status_rows: list[dict] = []
    safe_sample = _SAFE_FILE.sub("_", sample_id).strip("_") or "sample"
    query_name = f"QUERY__{safe_sample}"
    placed_loci: list[str] = []

    for locus_id in sorted(references):
        query_sequence = query_sequences.get(locus_id, "")
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
        _write_fasta(references[locus_id], reference_fasta)
        _run_mafft(
            build_mafft_reference_command(reference_fasta, threads, mafft),
            reference_alignment,
            f"reference alignment for {locus_id}",
        )
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
        _run_mafft(
            build_mafft_add_command(query_fasta, reference_alignment, threads, mafft),
            placed_alignment,
            f"--add placement for {locus_id}",
        )
        aligned_query_records = [
            (name, sequence)
            for name, sequence in _read_fasta(placed_alignment)
            if name == query_name
        ]
        if len(aligned_query_records) != 1:
            raise RuntimeError(
                f"MAFFT placement alignment for {locus_id} did not contain exactly one {query_name!r}"
            )
        _write_fasta(aligned_query_records, query_alignment)
        jplace_path = _run_epa_ng(
            build_epa_ng_command(
                reference_alignment,
                reference_tree,
                query_alignment,
                reference_model,
                epa_outdir,
                threads,
                epa_ng,
            ),
            epa_outdir,
            f"placement for {locus_id}",
        )
        placement_distances, placement_metadata = read_epa_ng_placement(jplace_path)
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
                    "placement_edge": placement_metadata["placement_edge"],
                    "like_weight_ratio": f'{placement_metadata["like_weight_ratio"]:.8f}',
                    "pendant_length": f'{placement_metadata["pendant_length"]:.8f}',
                    "distal_length": f'{placement_metadata["distal_length"]:.8f}',
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
    totals: dict[str, list[float]] = {}
    for row in detail_rows:
        totals.setdefault(str(row["reference_id"]), []).append(float(row["phylogenetic_distance"]))
    # A partial reference must not win merely because absent loci contribute
    # nothing. Rank only references represented in every successfully placed
    # locus, so the total is comparable across candidates.
    complete_totals = {
        reference_id: values
        for reference_id, values in totals.items()
        if len(values) == len(placed_loci)
    }
    summary_rows = []
    ordered = sorted(complete_totals.items(), key=lambda item: (sum(item[1]), item[0]))
    for rank, (reference_id, values) in enumerate(ordered, start=1):
        summary_rows.append(
            {
                "sample_id": sample_id,
                "reference_id": reference_id,
                "total_phylogenetic_distance": f"{sum(values):.8f}",
                "compared_loci": len(values),
                "mean_phylogenetic_distance": f"{sum(values) / len(values):.8f}",
                "rank": rank,
            }
        )
    summary_path = output / "phylogenetic_matches.tsv"
    _write_tsv(summary_rows, summary_path, SUMMARY_FIELDS)
    status_path = output / "locus_status.tsv"
    _write_tsv(status_rows, status_path, STATUS_FIELDS)
    return {
        "phylogeny": output,
        "phylogenetic_distances": detail_path,
        "phylogenetic_matches": summary_path,
        "phylogenetic_status": status_path,
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
