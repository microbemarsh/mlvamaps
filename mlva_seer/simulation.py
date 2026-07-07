from __future__ import annotations

import random
from pathlib import Path

from .io import read_loci, read_profiles, write_fastq, write_tsv
from .models import Locus, ReadRecord
from .sequence import revcomp


def _profile_counts(loci: list[Locus], profiles: list[dict], profile_id: str | None) -> dict[str, int]:
    if profiles:
        selected = profiles[0]
        if profile_id:
            selected = next((row for row in profiles if row.get("profile_id") == profile_id), profiles[0])
        counts = {}
        for locus in loci:
            value = selected.get(locus.locus_id, "")
            counts[locus.locus_id] = int(float(value)) if value not in ("", None) else _midpoint(locus)
        return counts
    return {locus.locus_id: _midpoint(locus) for locus in loci}


def _midpoint(locus: Locus) -> int:
    return int(round((locus.expected_min_repeats + locus.expected_max_repeats) / 2))


def _amplicon(locus: Locus, repeat_count: int) -> str:
    return (
        locus.forward_primer
        + locus.left_flank_sequence
        + (locus.repeat_motif * repeat_count)
        + locus.right_flank_sequence
        + revcomp(locus.reverse_primer)
    )


def _mutate(sequence: str, quality: str, error_rate: float, rng: random.Random) -> tuple[str, str]:
    bases = "ACGT"
    out_seq = []
    out_qual = []
    for base, qchar in zip(sequence, quality):
        roll = rng.random()
        if roll < error_rate / 5:
            continue
        if roll < error_rate / 5 + error_rate / 5:
            inserted = rng.choice(bases.replace(base, "") or bases)
            out_seq.append(inserted)
            out_qual.append(qchar)
        if roll < error_rate:
            out_seq.append(rng.choice(bases.replace(base, "") or bases))
        else:
            out_seq.append(base)
        out_qual.append(qchar)
    return "".join(out_seq), "".join(out_qual)


def simulate_reads(
    loci_path: str,
    outdir: str,
    sample_id: str,
    profiles_path: str | None = None,
    profile_id: str | None = None,
    depth: int = 200,
    error_rate: float = 0.03,
    seed: int = 13,
) -> dict[str, Path]:
    rng = random.Random(seed)
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    loci = read_loci(loci_path)
    profiles = read_profiles(profiles_path)
    counts = _profile_counts(loci, profiles, profile_id)
    reads = []
    truth = {"sample_id": sample_id}
    for locus in loci:
        count = counts[locus.locus_id]
        truth[locus.locus_id] = count
        template = _amplicon(locus, count)
        for idx in range(depth):
            quality = "I" * len(template)
            seq, qual = _mutate(template, quality, error_rate, rng)
            if rng.random() < 0.5:
                seq = revcomp(seq)
                qual = qual[::-1]
            reads.append(ReadRecord(f"{sample_id}_{locus.locus_id}_{idx}", seq, qual))
    rng.shuffle(reads)
    reads_path = outdir_path / f"{sample_id}.fastq.gz"
    write_fastq(reads, reads_path)
    truth_path = outdir_path / "truth_profile.tsv"
    write_tsv([truth], truth_path, ["sample_id"] + [locus.locus_id for locus in loci])
    return {"reads": reads_path, "truth": truth_path}
