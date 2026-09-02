from __future__ import annotations

import csv
import gzip
import itertools
import re
import warnings
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
from pathlib import Path

import regex

from .concurrency import DEFAULT_THREADS, resolve_threads
from .io import open_text
from .models import Locus
from .primers import read_loci_or_primers
from .sassy_cli import Searcher


PCR_TSV_FIELDS = [
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

_IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}
_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_POSITION_RE = re.compile(r"\bpos=(\d+)-(\d+)\b")


@dataclass(frozen=True)
class _PrimerMatch:
    start: int
    end: int
    cost: int
    cigar: str

    @property
    def length(self) -> int:
        return self.end - self.start


def _reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(_COMPLEMENT)[::-1]


def _read_fasta(path: str | Path):
    """Stream FASTA without routing sequence data back through Python objects twice."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    name: str | None = None
    parts: list[str] = []
    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].split()[0]
                if not name:
                    raise ValueError(f"Empty FASTA identifier at {path}:{line_number}")
                parts = []
            elif name is None:
                raise ValueError(f"Expected a FASTA header at {path}:{line_number}")
            else:
                parts.append(line)
    if name is not None:
        yield name, "".join(parts).upper()


def _expand_degenerate(primer: str) -> tuple[str, ...]:
    """Expand an IUPAC primer exactly as MLVA_finder does, deterministically."""
    try:
        choices = [_IUPAC[base] for base in primer.upper()]
    except KeyError as exc:
        raise ValueError(f"Unsupported primer base {exc.args[0]!r} in {primer!r}") from exc
    return tuple("".join(parts) for parts in itertools.product(*choices))


def _new_searcher():
    # ASCII is intentional. MLVA_finder expands ambiguity in the primer, while
    # an N in an assembly is an error rather than an IUPAC wildcard.
    try:
        return Searcher("ascii", rc=False)
    except TypeError:  # pragma: no cover - compatibility with early bindings
        return Searcher("ascii")


def _match_key(match: _PrimerMatch) -> tuple[int, int, int, str]:
    return (match.start, match.end, match.cost, match.cigar)


def _sassy_matches(searcher, pattern: str, sequence: str, max_errors: int) -> list[_PrimerMatch]:
    """Discover with Sassy, then resolve legacy ``regex`` traceback ties."""
    pattern_bytes = pattern.encode("ascii")
    sequence_bytes = sequence.encode("ascii")
    if hasattr(searcher, "search_all_alignments"):
        raw = itertools.chain.from_iterable(
            searcher.search_all_alignments(pattern_bytes, sequence_bytes, k=max_errors)
        )
    elif hasattr(searcher, "search_all"):  # pragma: no cover - binding compatibility
        raw = searcher.search_all(pattern_bytes, sequence_bytes, k=max_errors)
    else:  # pragma: no cover - retained for a useful dependency error
        raw = searcher.search(pattern_bytes, sequence_bytes, k=max_errors)

    sassy_matches: dict[tuple[int, int], _PrimerMatch] = {}
    candidate_starts: set[int] = set()
    for item in raw:
        start = int(item.text_start)
        end = int(item.text_end)
        pattern_start = int(getattr(item, "pattern_start", 0))
        pattern_end = int(getattr(item, "pattern_end", len(pattern)))
        if start < 0 or end <= start or pattern_start != 0 or pattern_end != len(pattern):
            continue
        candidate = _PrimerMatch(start, end, int(item.cost), str(item.cigar))
        candidate_starts.add(start)
        key = (start, end)
        previous = sassy_matches.get(key)
        if previous is None or (candidate.cost, candidate.cigar) < (previous.cost, previous.cigar):
            sassy_matches[key] = candidate

    # MLVA_finder used the third-party ``regex`` module. Its choice between
    # equal-cost insertion/deletion tracebacks differs from Sassy's. Running
    # anchored matches only at SIMD-discovered starts preserves that observable
    # behavior without making regex scan the genome.
    expression = regex.compile(f"({pattern}){{e<={max_errors}}}")
    resolved: list[_PrimerMatch] = []
    legacy_starts = sorted({
        nearby
        for start in candidate_starts
        for nearby in range(max(0, start - max_errors), min(len(sequence), start + max_errors + 1))
    })
    intervals: list[tuple[int, int]] = []
    for start in legacy_starts:
        if intervals and start <= intervals[-1][1] + 1:
            intervals[-1] = (intervals[-1][0], start)
        else:
            intervals.append((start, start))
    seen_legacy_spans: set[tuple[int, int]] = set()
    for interval_start, interval_end in intervals:
        window_end = min(len(sequence), interval_end + len(pattern) + max_errors)
        window = sequence[interval_start:window_end]
        for legacy in expression.finditer(window, overlapped=True):
            actual_start = interval_start + legacy.start()
            actual_end = interval_start + legacy.end()
            if actual_start > interval_end or (actual_start, actual_end) in seen_legacy_spans:
                continue
            seen_legacy_spans.add((actual_start, actual_end))
            observed = legacy.group(1)
            # Preserve MLVA_finder.positionsOfMatches: repeated identical fuzzy
            # strings are assigned the position of their first occurrence.
            start = actual_start if max_errors == 0 else sequence.find(observed)
            end = start + len(observed)
            substitutions, insertions, deletions = legacy.fuzzy_counts
            cost = substitutions + insertions + deletions
            traced = sassy_matches.get((actual_start, actual_end))
            cigar = traced.cigar if traced is not None else f"{len(observed)}M"
            resolved.append(_PrimerMatch(start, end, cost, cigar))
    return sorted(resolved, key=_match_key)


def _legacy_primer_matches(
    searcher, primer: str, sequence: str, max_errors: int
) -> list[_PrimerMatch]:
    """Match an IUPAC primer with MLVA_finder's equal-length preference."""
    matches: dict[tuple[int, int], _PrimerMatch] = {}
    for concrete in _expand_degenerate(primer):
        concrete_matches = _sassy_matches(
            searcher, concrete, sequence, max_errors
        )
        equal_length = [
            match for match in concrete_matches if match.length == len(concrete)
        ]
        for match in equal_length or concrete_matches:
            key = (match.start, match.end)
            previous = matches.get(key)
            if previous is None or (match.cost, match.cigar) < (previous.cost, previous.cigar):
                matches[key] = match
    return sorted(matches.values(), key=_match_key)


def _legacy_first_primer_orientations(
    searcher,
    primer: str,
    original: str,
    reverse: str,
    max_errors: int,
    search_rc: bool,
) -> list[tuple[str, str, list[_PrimerMatch]]]:
    """Preserve MLVA_finder's reverse fallback for each IUPAC expansion."""
    forward: dict[tuple[int, int], _PrimerMatch] = {}
    reverse_matches: dict[tuple[int, int], _PrimerMatch] = {}
    for concrete in _expand_degenerate(primer):
        concrete_forward = _sassy_matches(
            searcher, concrete, original, max_errors
        )
        equal_forward = [
            match
            for match in concrete_forward
            if match.length == len(concrete)
        ]
        concrete_forward = equal_forward or concrete_forward
        target = forward
        concrete_matches = concrete_forward
        if not concrete_matches and search_rc:
            concrete_matches = _sassy_matches(
                searcher, concrete, reverse, max_errors
            )
            equal_reverse = [
                match
                for match in concrete_matches
                if match.length == len(concrete)
            ]
            concrete_matches = equal_reverse or concrete_matches
            target = reverse_matches
        for match in concrete_matches:
            key = (match.start, match.end)
            previous = target.get(key)
            if previous is None or (match.cost, match.cigar) < (
                previous.cost,
                previous.cigar,
            ):
                target[key] = match
    orientations = []
    if forward:
        orientations.append(
            ("+", original, sorted(forward.values(), key=_match_key))
        )
    if reverse_matches:
        orientations.append(
            ("-", reverse, sorted(reverse_matches.values(), key=_match_key))
        )
    return orientations


def write_primer_pairs(loci: list[Locus], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "forward", "reverse"])
        writer.writeheader()
        for locus in loci:
            writer.writerow(
                {"name": locus.locus_id, "forward": locus.forward_primer, "reverse": locus.reverse_primer}
            )
    return path


def expected_amplicon_bounds(loci: list[Locus]) -> tuple[int, int]:
    mins = [locus.expected_amplicon_min_bp for locus in loci if locus.expected_amplicon_min_bp > 0]
    maxes = [locus.expected_amplicon_max_bp for locus in loci if locus.expected_amplicon_max_bp > 0]
    return (min(mins) if mins else 50, max(maxes) if maxes else 5000)


def _locus_product_rows(
    source: str,
    reference_id: str,
    original: str,
    locus: Locus,
    max_errors: int,
    min_len: int,
    max_len: int,
    search_rc: bool,
    trim_primers: bool,
    max_n_fraction: float,
    searcher=None,
) -> list[dict[str, str | int]]:
    searcher = searcher or _new_searcher()
    rows: list[dict[str, str | int]] = []
    seen_products: set[tuple[str, int, int, int, int]] = set()
    original = original.upper()
    reverse = _reverse_complement(original)
    # The legacy program reruns unresolved loci at thresholds 0..k. Retain the
    # union and let the assembly caller choose the first successful round.
    for error_round in range(max_errors + 1):
        orientations = _legacy_first_primer_orientations(
            searcher,
            locus.forward_primer,
            original,
            reverse,
            error_round,
            search_rc,
        )
        for strand, sequence, first_matches in orientations:
            second_matches = _legacy_primer_matches(
                searcher,
                _reverse_complement(locus.reverse_primer),
                sequence,
                error_round,
            )
            for first in first_matches:
                for second in second_matches:
                    product_key = (
                        strand,
                        first.start,
                        first.end,
                        second.start,
                        second.end,
                    )
                    if product_key in seen_products:
                        continue
                    seen_products.add(product_key)
                    row = _paired_product_row(
                        source,
                        reference_id,
                        original,
                        locus,
                        strand,
                        sequence,
                        first,
                        second,
                        min_len,
                        max_len,
                        trim_primers,
                        max_n_fraction,
                        len(rows) + 1,
                    )
                    if row is not None:
                        rows.append(row)
    return rows


def _record_product_rows(
    task: tuple[
        str,
        str,
        str,
        list[Locus],
        int,
        int,
        int,
        bool,
        bool,
        float,
    ],
) -> list[dict[str, str | int]]:
    (
        source,
        reference_id,
        original,
        loci,
        max_errors,
        min_len,
        max_len,
        search_rc,
        trim_primers,
        max_n_fraction,
    ) = task
    searcher = _new_searcher()
    rows = []
    for locus in loci:
        rows.extend(
            _locus_product_rows(
                source,
                reference_id,
                original,
                locus,
                max_errors,
                min_len,
                max_len,
                search_rc,
                trim_primers,
                max_n_fraction,
                searcher,
            )
        )
    return rows


def _bounded_ordered_map(executor, function, iterable, max_pending: int):
    """Map with bounded input consumption while preserving FASTA order."""
    iterator = iter(enumerate(iterable))
    pending = {}
    buffered = {}
    next_output = 0

    def submit_next() -> bool:
        try:
            index, item = next(iterator)
        except StopIteration:
            return False
        pending[executor.submit(function, item)] = index
        return True

    for _ in range(max_pending):
        if not submit_next():
            break
    while pending:
        completed, _remaining = wait(pending, return_when=FIRST_COMPLETED)
        for future in completed:
            buffered[pending.pop(future)] = future.result()
            submit_next()
        while next_output in buffered:
            yield buffered.pop(next_output)
            next_output += 1


def _product_rows(
    input_path: str | Path,
    loci: list[Locus],
    max_errors: int,
    min_len: int,
    max_len: int,
    search_rc: bool,
    trim_primers: bool,
    max_n_fraction: float,
    threads: int,
) -> list[dict[str, str | int]]:
    source = str(input_path)
    thread_count = resolve_threads(threads)
    record_iterator = iter(_read_fasta(input_path))
    try:
        first_record = next(record_iterator)
    except StopIteration:
        return []
    sentinel = object()
    second_record = next(record_iterator, sentinel)

    def record_tasks():
        initial_records = (
            (first_record,)
            if second_record is sentinel
            else (first_record, second_record)
        )
        for reference_id, original in itertools.chain(
            initial_records, record_iterator
        ):
            yield (
                source,
                reference_id,
                original,
                loci,
                max_errors,
                min_len,
                max_len,
                search_rc,
                trim_primers,
                max_n_fraction,
            )

    rows = []
    if (
        thread_count > 1
        and second_record is sentinel
        and len(loci) > 1
    ):
        reference_id, original = first_record

        def locus_tasks():
            for locus in loci:
                yield (
                    source,
                    reference_id,
                    original,
                    [locus],
                    max_errors,
                    min_len,
                    max_len,
                    search_rc,
                    trim_primers,
                    max_n_fraction,
                )

        with ProcessPoolExecutor(
            max_workers=min(thread_count, len(loci))
        ) as executor:
            for locus_rows in _bounded_ordered_map(
                executor,
                _record_product_rows,
                locus_tasks(),
                max_pending=max(2, thread_count * 2),
            ):
                rows.extend(locus_rows)
    elif thread_count == 1:
        record_results = map(_record_product_rows, record_tasks())
        for record_rows in record_results:
            rows.extend(record_rows)
    else:
        with ProcessPoolExecutor(max_workers=thread_count) as executor:
            for record_rows in _bounded_ordered_map(
                executor,
                _record_product_rows,
                record_tasks(),
                max_pending=max(2, thread_count * 2),
            ):
                rows.extend(record_rows)
    for ordinal, row in enumerate(rows, start=1):
        row["amplicon_id"] = (
            f"{row['reference_id']}:{row['primer_name']}:{ordinal}"
            + ("_rc" if row["strand"] == "-" else "")
        )
    return rows


def _paired_product_row(
    source: str,
    reference_id: str,
    original: str,
    locus: Locus,
    strand: str,
    sequence: str,
    first: _PrimerMatch,
    second: _PrimerMatch,
    min_len: int,
    max_len: int,
    trim_primers: bool,
    max_n_fraction: float,
    ordinal: int,
) -> dict[str, str | int] | None:
    """Build one compatibility-shaped row for existing downstream table readers."""
    legacy_len = (
        second.start
        + len(locus.reverse_primer)
        - first.start
        - first.length
        + len(locus.forward_primer)
    )
    if legacy_len <= 0:
        return None
    observed_start = first.start
    observed_end = second.end
    if observed_end <= observed_start:
        return None
    full_product = sequence[observed_start:observed_end]
    if not min_len <= legacy_len <= max_len:
        return None
    if full_product and full_product.count("N") / len(full_product) > max_n_fraction:
        return None
    if strand == "+":
        original_start, original_end = observed_start, observed_end
    else:
        original_start = len(original) - observed_end
        original_end = len(original) - observed_start
    product = sequence[first.end:second.start] if trim_primers else full_product
    amplicon_id = (
        f"{reference_id}:{locus.locus_id}:{ordinal}" + ("_rc" if strand == "-" else "")
    )
    return {
        "amplicon_id": amplicon_id,
        "reference_id": reference_id,
        "source_file": source,
        "primer_name": locus.locus_id,
        "product_len": len(product),
        "full_len": observed_end - observed_start,
        "fwd_start": first.start,
        "fwd_end": first.end,
        "fwd_mismatches": first.cost,
        "fwd_identity": f"{1 - first.cost / max(len(locus.forward_primer), 1):.6f}",
        "fwd_cigar": first.cigar,
        "rev_start": second.start,
        "rev_end": second.end,
        "rev_mismatches": second.cost,
        "rev_identity": f"{1 - second.cost / max(len(locus.reverse_primer), 1):.6f}",
        "rev_cigar": second.cigar,
        "strand": strand,
        "is_circular_wrap": "false",
        "product_seq": product,
        "original_start": original_start,
        "original_end": original_end,
    }


def run_in_silico_pcr(
    input_path: str | Path,
    loci_path: str | Path | None,
    outdir: str | Path,
    primers_path: str | Path | None = None,
    **kwargs,
) -> dict[str, Path]:
    loci = read_loci_or_primers(loci_path, primers_path)
    return run_in_silico_pcr_loci(input_path, loci, outdir, **kwargs)


def run_in_silico_pcr_loci(
    input_path: str | Path,
    loci: list[Locus],
    outdir: str | Path,
    max_errors: int = 2,
    threads: int = DEFAULT_THREADS,
    circular: bool = False,
    search_rc: bool = True,
    trim_primers: bool = False,
    max_n_fraction: float = 0.0,
    amplicon_bounds: tuple[int, int] | None = None,
    **_compatibility_options,
) -> dict[str, Path]:
    """Run mlvamaps' Sassy-backed, MLVA_finder-compatible in silico PCR."""
    if not 0 <= max_n_fraction <= 1:
        raise ValueError("max_n_fraction must be between 0 and 1")
    if max_errors < 0:
        raise ValueError("max_errors must be non-negative")
    resolve_threads(threads)
    if circular:
        warnings.warn(
            "Circular-wrap products are not yet emitted by the MLVA_finder compatibility engine",
            RuntimeWarning,
            stacklevel=2,
        )

    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    primers_output = write_primer_pairs(loci, outdir_path / "primers.csv")
    products_path = outdir_path / "products.fasta.gz"
    stats_path = outdir_path / "matches.tsv"
    min_len, max_len = amplicon_bounds or expected_amplicon_bounds(loci)
    rows = _product_rows(
        input_path,
        loci,
        max_errors,
        min_len,
        max_len,
        search_rc,
        trim_primers,
        max_n_fraction,
        threads,
    )
    with open_text(products_path, "wt") as handle:
        for row in rows:
            handle.write(
                f">{row['amplicon_id']}\tpos={row['original_start']}-{row['original_end']}"
                f"\tstrand={row['strand']}\tlen={row['full_len']}\n{row['product_seq']}\n"
            )
    with stats_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PCR_TSV_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return {"primers": primers_output, "products": products_path, "stats": stats_path}


def read_pcr_results(stats_path: str | Path, products_path: str | Path) -> list[dict[str, str | int]]:
    coordinates: dict[str, tuple[int, int]] = {}
    with open_text(products_path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                match = _POSITION_RE.search(line)
                if match:
                    coordinates[line[1:].split("\t", 1)[0]] = (int(match.group(1)), int(match.group(2)))
    rows: list[dict[str, str | int]] = []
    with Path(stats_path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            amplicon_id = row.get("amplicon_id", "")
            if amplicon_id not in coordinates:
                raise RuntimeError(f"PCR match {amplicon_id!r} has no matching FASTA record")
            start, end = coordinates[amplicon_id]
            rows.append({**row, "original_start": start, "original_end": end})
    return rows


# Temporary source-compatible aliases for callers built against mlvamaps 0.1.
# They no longer invoke or require the Amplirust executable.
AMPLIRUST_TSV_FIELDS = PCR_TSV_FIELDS
write_amplirust_primers = write_primer_pairs
run_amplirust = run_in_silico_pcr
run_amplirust_loci = run_in_silico_pcr_loci
read_amplirust_results = read_pcr_results


def build_amplirust_command(*_args, **_kwargs):
    raise RuntimeError("Amplirust was replaced by mlvamaps' built-in Sassy-backed PCR engine")
