"""Deterministic, organism-independent synthetic FASTQ benchmark inputs."""
import json
import random


def prepare_challenge(directory, layout, copies=1):
    from mlvamaps.io import write_fastq, write_tsv
    from mlvamaps.models import Locus, ReadRecord
    from mlvamaps.sequence import revcomp
    from mlvamaps.candidate_contexts import generate_candidate_contexts, write_candidate_contexts
    from mlvamaps.minimap_mapping import build_competitive_indexes

    directory.mkdir(parents=True)
    database = directory / "database"
    database.mkdir()
    rng = random.Random(731)
    # name, motif, copies, coverage, per-base error, interruption, partial bases
    cases = [
        ("canonical", "ATGC", 5, 12, 0, False, 0),
        ("expansion", "ATGC", 40, 12, 0, False, 0),
        ("long_array", "ATGC", 400, 12, 0, False, 0),
        ("contraction", "ATGC", 1, 12, 0, False, 0),
        ("interrupted", "ATGC", 8, 12, 0, True, 0),
        ("partial_unit", "ATGC", 5, 12, 0, False, 1),
        ("homopolymer", "A", 16, 12, .02, False, 0),
        ("long_motif", "AGCTACGGA", 8, 12, .01, False, 0),
        ("noisy", "ATGC", 16, 20, .04, False, 0),
        ("low_coverage", "ATGC", 6, 1, 0, False, 0),
        ("zero_coverage", "ATGC", 6, 0, 0, False, 0),
        ("absent", "ATGC", 0, 0, 0, False, 0),
        ("mixture", "ATGC", 5, 20, 0, False, 0),
        ("ambiguous_a", "ATGC", 5, 10, 0, False, 0),
        ("ambiguous_b", "ATGC", 5, 10, 0, False, 0),
    ]
    loci, reads1, reads2, truth, read_truth = [], [], [], {}, {}
    previous_anchors = None

    def errors(sequence, rate):
        observed = []
        for base in sequence:
            if rng.random() >= rate:
                observed.append(base)
            else:
                operation = rng.choice(("substitute", "insert", "delete"))
                if operation == "substitute":
                    observed.append(rng.choice([value for value in "ACGT" if value != base]))
                elif operation == "insert":
                    observed.extend((base, rng.choice("ACGT")))
        return "".join(observed)

    for name, motif, count, coverage, rate, interrupted, partial in cases:
        anchors = ["".join(rng.choices("ACGT", k=size)) for size in (20, 50, 50, 20)]
        if name == "ambiguous_b":
            anchors = previous_anchors
        if name == "ambiguous_a":
            previous_anchors = anchors
        forward, left, right, reverse = anchors
        locus = Locus(name, forward_primer=forward, reverse_primer=reverse,
            left_flank_sequence=left, right_flank_sequence=right, repeat_motif=motif,
            repeat_unit_length_bp=len(motif), expected_min_repeats=0, expected_max_repeats=12,
            expected_product_size_bp=140 + 4 * len(motif), nominal_repeat_units=4)
        loci.append(locus)
        reference = forward + left + motif * 4 + right + revcomp(reverse)
        alternative = reference[:25] + ("A" if reference[25] != "A" else "C") + reference[26:]
        (database / f"{name}.fasta").write_text(f">R1\n{reference}\n>R1_copy\n{reference}\n>related_R2\n{alternative}\n")
        repeat = motif * count + motif[:partial]
        if interrupted:
            repeat = repeat[:9] + ("A" if repeat[9] != "A" else "C") + repeat[10:]
        product = forward + left + repeat + right + revcomp(reverse)
        other = forward + left + motif * 9 + right + revcomp(reverse)
        allele = count + .5 if partial else count
        truth[name] = dict(repeat_count=None if name == "absent" else allele,
            raw_repeat_count=count + partial / len(motif), coverage=coverage * copies,
            motif_length=len(motif), repeat_bases=len(repeat), error_rate=rate,
            products=[product, other] if name == "mixture" else [product],
            alleles=[5, 9] if name == "mixture" else [allele],
            expected_uncertainty=name.startswith("ambiguous") or name in {"zero_coverage", "absent", "mixture"})
        for index in range(coverage * copies):
            observed = other if name == "mixture" and index % 5 >= 3 else product
            identifier = f"{name}__{index}"
            read_truth[identifier] = name
            if layout == "long":
                first = errors(observed, rate)
                if index % 2:
                    first = revcomp(first)
                reads1.append(ReadRecord(identifier, first, ("5" if rate else "I") * len(first)))
            else:
                contig = "".join(rng.choices("ACGT", k=80)) + observed + "".join(rng.choices("ACGT", k=80))
                # Tile both boundaries even when the array exceeds a fragment.
                start = 40 + round(index % max(coverage, 1) * (len(contig) - 260) / max(1, coverage - 1))
                fragment = contig[start:start + 180]
                first, second = errors(fragment[:100], rate), errors(revcomp(fragment[-100:]), rate)
                reads1.append(ReadRecord(identifier + "/1", first, "I" * len(first)))
                reads2.append(ReadRecord(identifier + "/2", second, "I" * len(second)))
    panel = directory / "panel.tsv"
    write_tsv([vars(locus) for locus in loci], panel, list(vars(loci[0])))
    assets = write_candidate_contexts(generate_candidate_contexts(loci, database), database / "competitive_mapping")
    build_competitive_indexes(assets["fasta"], database / "competitive_mapping")
    (database / "manifest.json").write_text("{}")
    first_path, second_path = directory / "r1.fastq", directory / "r2.fastq"
    write_fastq(reads1, first_path)
    write_fastq(reads2, second_path)
    (directory / "truth.json").write_text(json.dumps(truth, indent=2) + "\n")
    (directory / "read_truth.json").write_text(json.dumps(read_truth, indent=2) + "\n")
    argv = ["call", "-p", str(panel), "--database", str(database), "-o", str(directory / "out"),
            "--sample-id", "sample", "--sample-mode", "isolate", "-t", "2", "--quiet"]
    if layout == "long":
        argv.extend(["-i", str(first_path), "--read-technology", "ont"])
    else:
        argv.extend(["-i", "sr", "--fq1", str(first_path), "--fq2", str(second_path)])
    return argv, [str(first_path)] + ([str(second_path)] if layout != "long" else []), truth
