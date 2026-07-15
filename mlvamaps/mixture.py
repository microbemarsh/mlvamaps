from __future__ import annotations

from collections import defaultdict

import numpy as np

from .clustering import _alignment_metrics


MIXTURE_FIELDS = [
    "sample_id",
    "locus_id",
    "variant_id",
    "repeat_count",
    "observed_reads",
    "observed_fraction",
    "estimated_reads",
    "estimated_fraction",
    "abundance_class",
    "meaningful",
    "meaningful_threshold",
    "adaptive_em_floor",
    "estimated_error_rate",
    "em_iterations",
    "em_converged",
]


def _estimated_error_rate(rows: list[dict]) -> float:
    edits = sum(
        int(row.get("total_insertions") or 0)
        + int(row.get("total_deletions") or 0)
        + int(row.get("total_substitutions") or 0)
        for row in rows
    )
    aligned_bases = sum(
        int(row.get("support_reads") or 0)
        * max(int(row.get("representative_length_bp") or 0), 1)
        for row in rows
    )
    # Jeffreys-style smoothing prevents a zero error probability while allowing
    # clean, deeply sequenced loci to retain a suitably sharp likelihood model.
    rate = (edits + 0.5) / (aligned_bases + 1.0)
    return min(0.15, max(1e-5, rate))


def _log_likelihood_matrix(rows: list[dict], error_rate: float) -> np.ndarray:
    sequences = [str(row.get("representative_sequence") or "") for row in rows]
    size = len(sequences)
    matrix = np.zeros((size, size), dtype=np.float64)
    edit_log_odds = np.log(error_rate / (3.0 * (1.0 - error_rate)))
    for observed_index, observed in enumerate(sequences):
        for source_index, source in enumerate(sequences):
            if observed_index == source_index:
                continue
            metrics = _alignment_metrics(observed, source)
            edit_distance = int(metrics["edit_distance_to_representative"])
            matrix[observed_index, source_index] = edit_distance * edit_log_odds
    return matrix


def _run_em(
    counts: np.ndarray,
    log_likelihoods: np.ndarray,
    initial: np.ndarray,
    active: np.ndarray,
    tolerance: float,
    max_iterations: int,
) -> tuple[np.ndarray, int, bool]:
    frequencies = np.where(active, initial, 0.0).astype(np.float64)
    frequencies /= frequencies.sum(dtype=np.float64)
    previous_likelihood = -np.inf

    for iteration in range(1, max_iterations + 1):
        log_frequencies = np.full(frequencies.shape, -np.inf, dtype=np.float64)
        positive = frequencies > 0
        log_frequencies[positive] = np.log(frequencies[positive])
        log_weights = log_likelihoods + log_frequencies[np.newaxis, :]
        row_maxima = np.max(log_weights, axis=1)
        scaled = np.exp(log_weights - row_maxima[:, np.newaxis])
        denominators = scaled.sum(axis=1)
        responsibilities = scaled / denominators[:, np.newaxis]
        updated = (counts[:, np.newaxis] * responsibilities).sum(axis=0)
        updated = np.where(active, updated, 0.0)
        updated /= updated.sum(dtype=np.float64)
        log_likelihood = float(
            np.sum(counts * (row_maxima + np.log(denominators)), dtype=np.float64)
        )

        frequency_change = float(np.max(np.abs(updated - frequencies)))
        likelihood_change = log_likelihood - previous_likelihood
        frequencies = updated
        if (
            frequency_change <= tolerance
            and np.isfinite(previous_likelihood)
            and abs(likelihood_change) <= tolerance
        ):
            return frequencies, iteration, True
        previous_likelihood = log_likelihood

    return frequencies, max_iterations, False


def _adaptive_em_floor(read_count: int) -> float:
    if read_count <= 0:
        return 0.0
    if read_count > 1000:
        return min(1.0, 10.0 / read_count)
    return 1.0 / (read_count + 1.0)


def estimate_variant_mixtures(
    asv_rows: list[dict],
    min_fraction: float = 0.01,
    tolerance: float = 1e-7,
    max_iterations: int = 200,
) -> list[dict]:
    """Estimate per-locus variant abundance with an Emu-style count EM model.

    VSEARCH support counts are the observations. Pairwise distances between
    observed representatives define read-assignment likelihoods, and EM
    alternates between abundance-dependent assignments and abundance updates.
    """
    if not 0.0 <= min_fraction <= 1.0:
        raise ValueError("min_fraction must be between 0 and 1")
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    if max_iterations < 1:
        raise ValueError("max_iterations must be at least 1")

    by_locus: dict[str, list[dict]] = defaultdict(list)
    for row in asv_rows:
        by_locus[str(row["locus_id"])].append(row)

    output: list[dict] = []
    for locus_id, unsorted_rows in sorted(by_locus.items()):
        rows = sorted(
            unsorted_rows,
            key=lambda row: (-int(row.get("support_reads") or 0), str(row["variant_id"])),
        )
        counts = np.asarray(
            [int(row.get("support_reads") or 0) for row in rows], dtype=np.float64
        )
        total_reads = int(counts.sum(dtype=np.float64))
        if total_reads <= 0:
            continue

        observed_fractions = counts / total_reads
        error_rate = _estimated_error_rate(rows)
        adaptive_floor = _adaptive_em_floor(total_reads)
        if len(rows) == 1:
            frequencies = np.asarray([1.0], dtype=np.float64)
            iterations = 0
            converged = True
        else:
            log_likelihoods = _log_likelihood_matrix(rows, error_rate)
            initial = (counts + 0.5) / (total_reads + (0.5 * len(rows)))
            active = np.ones(len(rows), dtype=bool)
            frequencies, iterations, converged = _run_em(
                counts,
                log_likelihoods,
                initial,
                active,
                tolerance,
                max_iterations,
            )
            retained = frequencies >= adaptive_floor
            if not np.any(retained):
                retained[int(np.argmax(frequencies))] = True
            if not np.all(retained):
                frequencies, extra_iterations, reconverged = _run_em(
                    counts,
                    log_likelihoods,
                    frequencies,
                    retained,
                    tolerance,
                    max_iterations,
                )
                iterations += extra_iterations
                converged = converged and reconverged

        meaningful = frequencies >= max(min_fraction, adaptive_floor)
        if not np.any(meaningful):
            meaningful[int(np.argmax(frequencies))] = True
        ranking = np.argsort(-frequencies, kind="stable")
        rank_by_index = {int(index): rank for rank, index in enumerate(ranking, start=1)}

        for index in ranking:
            index = int(index)
            fraction = float(frequencies[index])
            is_meaningful = bool(meaningful[index])
            if rank_by_index[index] == 1:
                abundance_class = "DOMINANT"
            elif is_meaningful:
                abundance_class = "SECONDARY"
            else:
                abundance_class = "TRACE"
            row = rows[index]
            output.append(
                {
                    "sample_id": row.get("sample_id", ""),
                    "locus_id": locus_id,
                    "variant_id": row["variant_id"],
                    "repeat_count": row.get("repeat_count", ""),
                    "observed_reads": int(counts[index]),
                    "observed_fraction": round(float(observed_fractions[index]), 8),
                    "estimated_reads": round(fraction * total_reads, 3),
                    "estimated_fraction": round(fraction, 8),
                    "abundance_class": abundance_class,
                    "meaningful": "yes" if is_meaningful else "no",
                    "meaningful_threshold": round(
                        max(min_fraction, adaptive_floor), 8
                    ),
                    "adaptive_em_floor": round(adaptive_floor, 8),
                    "estimated_error_rate": round(error_rate, 8),
                    "em_iterations": iterations,
                    "em_converged": "yes" if converged else "no",
                }
            )
    return output
