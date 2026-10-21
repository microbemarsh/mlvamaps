"""Recruit mapping-supported molecules and retain observed query products."""

from collections import defaultdict
from dataclasses import replace
import re

from .models import ReadPair, ReadRecord
from .short_read_fragments import RecruitedFragment

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def aligned_flanks(row, context):
    """Count aligned query bases outside the repeat; deletions are not anchors."""
    position = row.reference_start
    left = right = 0
    for size, operation in _CIGAR.findall(row.cigar):
        size = int(size)
        end = position + size
        if operation in "M=X":
            left += max(0, min(end, context.repeat_start) - position)
            right += max(0, end - max(position, context.repeat_end))
        if operation in "MDN=X":
            position = end
    return left, right


def recruitment_identity(row, context):
    """Do not penalize locus recruitment for an observed repeat-length change."""
    position, deletions, repeat_gaps = row.reference_start, 0, 0
    for size, operation in _CIGAR.findall(row.cigar):
        size = int(size)
        if operation == "I" and context.repeat_start <= position <= context.repeat_end:
            repeat_gaps += size
        if operation == "D":
            deletions += size
            repeat_gaps += max(0, min(position + size, context.repeat_end) - max(position, context.repeat_start))
        if operation in "MDN=X":
            position += size
    denominator = row.query_end - row.query_start + deletions
    return min(1.0, row.alignment_identity * denominator / max(1, denominator - repeat_gaps))


def recruit_mapping_fragments(alignments, contexts, fragments=None, unmapped_mates=None,
                              *, minimum_identity=.9, minimum_margin=10):
    """Resolve locus identity across references, retaining each molecule once.

    Reference alleles within one locus remain competitors. Repeat-only matches
    cannot recruit, and conflicting anchored mates cannot create a chimeric pool.
    """
    pools = {locus: dict(pool) for locus, pool in (fragments or {}).items()}
    context_by_id = {context.candidate_id: context for context in contexts}
    grouped = defaultdict(list)
    for row in alignments:
        grouped[row.molecule_id].append(row)
    seed_owners = {name: locus for locus, pool in pools.items() for name in pool}
    for molecule, rows in grouped.items():
        seed_locus = seed_owners.get(molecule)
        if seed_locus and all(row.locus_id == seed_locus for row in rows):
            # Independent seed recruitment already retained the complete pair;
            # same-locus allele competitors cannot change its locus assignment.
            continue
        scores, anchors, reads = defaultdict(dict), defaultdict(dict), {}
        for row in rows:
            mate = row.mate or 1
            if row.query_sequence and (mate not in reads or row.primary):
                reads[mate] = ReadRecord(row.read_id, row.query_sequence, row.query_quality)
            context = context_by_id[row.candidate_id]
            if row.alignment_identity < minimum_identity and recruitment_identity(row, context) < minimum_identity:
                continue
            scores[row.locus_id][mate] = max(scores[row.locus_id].get(mate, -float("inf")), row.alignment_score)
            left, right = aligned_flanks(row, context)
            if max(left, right) >= 21:
                old = anchors[row.locus_id].get(mate)
                if old is None or row.alignment_score > old[0]:
                    anchors[row.locus_id][mate] = (row.alignment_score, left >= 21, right >= 21)
        if not anchors:
            continue  # Keep independent seed/PCR recovery when mapping is uninformative.
        ranked = sorted(scores, key=lambda locus: (-sum(scores[locus].values()), locus))
        winner = ranked[0]
        ambiguous = (winner not in anchors or (len(ranked) > 1 and
                     sum(scores[winner].values()) - sum(scores[ranked[1]].values()) < minimum_margin))
        for mate in (1, 2):
            anchored = {locus: values[mate][0] for locus, values in anchors.items() if mate in values}
            if anchored and winner not in {locus for locus, score in anchored.items() if score == max(anchored.values())}:
                ambiguous = True
        if seed_locus and (ambiguous or seed_locus != winner):
            pools[seed_locus].pop(molecule, None)
            continue
        if ambiguous:
            continue
        if seed_locus:
            continue  # Already have both original mates, including unmapped sequence.
        for mate in (1, 2):
            recovered = (unmapped_mates or {}).get((molecule, mate))
            if mate not in reads and recovered is not None:
                reads[mate] = recovered
        if not reads:
            continue
        ordered = [reads[mate] for mate in sorted(reads)]
        support = anchors[winner].values()
        pools.setdefault(winner, {})[molecule] = RecruitedFragment(
            ReadPair(molecule, ordered[0], ordered[1] if len(ordered) > 1 else None),
            any(item[1] for item in support), any(item[2] for item in support),
            len(reads) > len(anchors[winner]), "mapping",
        )
    return pools


def call_query_products(queries, sample_id, round_tolerance, max_primer_mismatches=2):
    """Batch sample contigs through the assembly PCR selection and calibration."""
    from .assembly_call import legacy_amplicon_bounds, legacy_assembly_call_rows, pcr_rows_to_products
    from .in_silico_pcr import _record_batch_product_rows

    if not queries:
        return {}
    loci = [locus for locus, sequence in queries]
    tasks = [
        ("mapped_query", locus.locus_id, sequence, [locus], max_primer_mismatches,
         *legacy_amplicon_bounds([locus]), True, False, 1.0)
        for locus, sequence in queries
    ]
    # Bound the primer/contig cross product for large panels.
    rows = [row for start in range(0, len(tasks), 32)
            for row in _record_batch_product_rows(tasks[start:start + 32])]
    products = pcr_rows_to_products(rows, loci, sample_id,
                                   enforce_locus_bounds=False, measure_products=False)
    calls = legacy_assembly_call_rows(loci, products, sample_id, round_tolerance=round_tolerance)
    sequences = {product["product_id"]: product["sequence"] for product in products}
    return {call["locus_id"]: (call, sequences.get(call["evidence"], "")) for call in calls}


def resolve_query_products(calls, evidence, loci, assembled, sample_id, thresholds, round_tolerance,
                           max_primer_mismatches=2):
    """Select supported sample products; recalibrate their original observations.

    No reference bases are used. Partial-read products come from the covered
    graph; complete products require agreement among independent molecules.
    """
    by_locus = defaultdict(lambda: defaultdict(set))
    spanning = defaultdict(dict)
    for row in evidence:
        sequence = row.metadata.get("product_sequence")
        if sequence and row.direct_product_measurement:
            by_locus[row.locus_id][row.molecule_id].add(sequence)
            spanning[row.locus_id].setdefault(row.molecule_id, (sequence, "", row.measured_repeat_count))
        if row.metadata.get("locus_sequence"):
            spanning[row.locus_id][row.molecule_id] = (
                row.metadata["locus_sequence"], row.metadata.get("locus_quality") or "",
                row.measured_repeat_count,
            )
    call_by_locus = {row["locus"]: row for row in calls}
    assembled = dict(assembled)
    products, audit, replacements, pending = {}, [], {}, {}
    for locus in loci:
        call = call_by_locus[locus.locus_id]
        molecules = by_locus[locus.locus_id]
        groups = defaultdict(set)
        for molecule, sequences in molecules.items():
            if len(sequences) == 1:
                groups[next(iter(sequences))].add(molecule)
        sequence, support, method = "", set(), ""
        if groups:
            best = sorted(groups, key=lambda item: (-len(groups[item]), item))[0]
            if (len(groups[best]) >= thresholds.minimum_molecules
                    and len(groups[best]) / len(molecules) >= thresholds.minimum_probability
                    and sum(len(names) == len(groups[best]) for names in groups.values()) == 1):
                sequence, support, method = best, groups[best], "complete_molecules"
        observed = spanning[locus.locus_id]
        if not sequence and len(observed) >= thresholds.minimum_molecules and call["status"] != "mixed":
            # Cluster measured alleles before correction. Unmeasurable reads
            # may join only a narrow length mode; never pool distant products.
            clusters = defaultdict(set)
            for name, (segment, quality, allele) in observed.items():
                clusters[allele if allele is not None else ("length", len(segment))].add(name)
            ranked = sorted(clusters.values(), key=lambda names: (-len(names), sorted(names)))
            support = ranked[0]
            if (len(support) >= thresholds.minimum_molecules
                    and len(support) / len(observed) >= thresholds.minimum_probability
                    and (len(ranked) == 1 or len(support) > len(ranked[1]))):
                from .local_assembly import sequence_consensus

                # ponytail: bound POA to 100 independent reads. All original
                # molecules remain in inference; deeper correction needs benchmarking.
                selected = sorted(support, key=lambda name: (
                    -sum(map(ord, observed[name][1])) / max(1, len(observed[name][1])), name))[:100]
                sequence = sequence_consensus(observed[name][0] for name in selected)
                method = "mapped_flank_consensus"
            else:
                support = set()
        if not sequence and locus.locus_id in assembled:
            sequence, method = assembled[locus.locus_id], "mapped_local_assembly"
        status = "NO_SUPPORTED_PRODUCT"
        if call["status"] == "mixed":
            sequence, status = "", "MIXED_LOCUS"
        if sequence:
            pending[locus.locus_id] = (locus, sequence, support)
        audit.append(dict(sample_id=sample_id, locus_id=locus.locus_id, status=status,
                          method=method, complete_molecules=len(molecules),
                          supporting_complete_molecules=len(support & set(molecules)), repeat_count="",
                          spanning_molecules=len(observed),
                          reconstruction_molecules=len(support),
                          reconstruction_fraction=len(support) / max(1, len(observed) if method == "mapped_flank_consensus" else len(molecules)),
                          raw_repeat_count="", sequence_length=0))
    pcr_results = call_query_products([(locus, sequence) for locus, sequence, support in pending.values()],
                                      sample_id, round_tolerance, max_primer_mismatches)
    for row in audit:
        locus_id = row["locus_id"]
        if locus_id not in pcr_results:
            continue
        pcr, product = pcr_results[locus_id]
        row["status"] = "PCR_UNRESOLVED"
        if pcr["status"] == "PASS" and product and set(product) <= set("ACGT"):
            allele = pcr["repeat_count"]
            products[locus_id] = (product, allele)
            row.update(status="RECONSTRUCTED", repeat_count=allele,
                       raw_repeat_count=pcr["repeat_count_raw"], sequence_length=len(product))
            # Exact observed products can recalibrate their own measurements.
            # A derived consensus is checked separately and creates no votes.
            for molecule in pending[locus_id][2] if row["method"] == "complete_molecules" else ():
                replacements[(locus_id, molecule)] = (allele, float(pcr["repeat_count_raw"]))
    changed = False
    updated = []
    for row in evidence:
        replacement = replacements.get((row.locus_id, row.molecule_id))
        if replacement and row.measured_repeat_count is not None:
            allele, raw = replacement
            if row.measured_repeat_count != allele or row.raw_repeat_count != raw:
                row = replace(row, measured_repeat_count=allele, raw_repeat_count=raw)
                changed = True
        updated.append(row)
    return products, audit, updated, changed


def apply_reconstructed_calls(calls, audit, thresholds):
    """Use supported PCR products without turning consensus into extra reads."""
    reconstructed = {row["locus_id"]: row for row in audit if row["status"] == "RECONSTRUCTED"}
    for call in calls:
        row = reconstructed.get(call["locus"])
        if not row or call["status"] == "mixed":
            continue
        if row["method"] != "mapped_flank_consensus":
            continue
        allele = row["repeat_count"]
        if call["status"] in {"called", "low_coverage"} and call["repeat_count"] != allele:
            call.update(status="ambiguous", repeat_count="", confidence="unresolved",
                        consensus_status="disagrees", estimation_method="",
                        estimation_reason="Supported PCR product conflicts with molecule measurements.")
            continue
        if row["reconstruction_molecules"] < thresholds.minimum_molecules:
            continue
        fraction = row["reconstruction_fraction"]
        call.update(status="called", repeat_count=allele, repeat_count_raw=row["raw_repeat_count"],
                    best_candidate_repeat=allele, confidence="moderate", consensus_status="agrees",
                    inference_method="mapped_consensus_assembly_pcr", estimation_method="", estimation_reason="",
                    probability_kind="reconstruction_support_fraction",
                    molecule_support=max(call["molecule_support"], row["reconstruction_molecules"]),
                    dominant_repeat=allele, dominant_fraction=fraction,
                    candidate_distribution=f"{allele}:{fraction:.8f}",
                    best_probability=fraction, second_best_probability=0.0, margin=0.0)
