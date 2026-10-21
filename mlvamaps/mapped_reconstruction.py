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
            if row.alignment_identity < minimum_identity:
                continue
            scores[row.locus_id][mate] = max(scores[row.locus_id].get(mate, -float("inf")), row.alignment_score)
            left, right = aligned_flanks(row, context_by_id[row.candidate_id])
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


def call_query_products(queries, sample_id, round_tolerance):
    """Batch sample contigs through the assembly PCR selection and calibration."""
    from .assembly_call import legacy_amplicon_bounds, legacy_assembly_call_rows, pcr_rows_to_products
    from .in_silico_pcr import _record_batch_product_rows

    if not queries:
        return {}
    loci = [locus for locus, sequence in queries]
    tasks = [
        ("mapped_query", locus.locus_id, sequence, [locus], 2,
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


def resolve_query_products(calls, evidence, loci, assembled, sample_id, thresholds, round_tolerance):
    """Select supported sample products; recalibrate their original observations.

    No reference bases are used. Partial-read products come from the covered
    graph; complete products require agreement among independent molecules.
    """
    by_locus = defaultdict(lambda: defaultdict(set))
    for row in evidence:
        sequence = row.metadata.get("product_sequence")
        if sequence and row.direct_product_measurement:
            by_locus[row.locus_id][row.molecule_id].add(sequence)
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
        # ponytail: this fast path retains an observed dominant haplotype. Noisy
        # complete products without one need the existing long-read POA path;
        # averaging unaligned sequences here could invent a chimeric product.
        if groups:
            best = sorted(groups, key=lambda item: (-len(groups[item]), item))[0]
            if (len(groups[best]) >= thresholds.minimum_molecules
                    and len(groups[best]) / len(molecules) >= thresholds.minimum_probability
                    and sum(len(names) == len(groups[best]) for names in groups.values()) == 1):
                sequence, support, method = best, groups[best], "complete_molecules"
        if not sequence and locus.locus_id in assembled:
            sequence, method = assembled[locus.locus_id], "mapped_local_assembly"
        status = "NO_SUPPORTED_PRODUCT"
        if call["status"] == "mixed":
            sequence, status = "", "MIXED_LOCUS"
        if sequence:
            pending[locus.locus_id] = (locus, sequence, support)
        audit.append(dict(sample_id=sample_id, locus_id=locus.locus_id, status=status,
                          method=method, complete_molecules=len(molecules),
                          supporting_complete_molecules=len(support), repeat_count="",
                          sequence_length=0))
    pcr_results = call_query_products([(locus, sequence) for locus, sequence, support in pending.values()],
                                      sample_id, round_tolerance)
    for row in audit:
        locus_id = row["locus_id"]
        if locus_id not in pcr_results:
            continue
        pcr, product = pcr_results[locus_id]
        row["status"] = "PCR_UNRESOLVED"
        if pcr["status"] == "PASS" and product and set(product) <= set("ACGT"):
            allele = pcr["repeat_count"]
            products[locus_id] = (product, allele)
            row.update(status="RECONSTRUCTED", repeat_count=allele, sequence_length=len(product))
            for molecule in pending[locus_id][2]:
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
