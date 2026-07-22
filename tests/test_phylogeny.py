from __future__ import annotations

import csv
import io
import json
from pathlib import Path
from threading import Barrier

import pytest

from mlvamaps import phylogeny as phylogeny_module
from mlvamaps.phylogeny import (
    _parse_newick,
    _placement_patristic_distances,
    _run_raxml_ng,
    _tip_patristic_distances,
    build_mafft_add_command,
    build_mafft_reference_command,
    build_epa_ng_command,
    build_raxml_ng_command,
    read_sequence_database,
    read_epa_ng_placement,
    read_epa_ng_placement_statistics,
    neighbor_joining_tree,
    run_phylogenetic_placement,
    decompose_marker_sequence,
)
from mlvamaps.models import Locus
from mlvamaps.progress import ProgressReporter


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _fake_mafft(tmp_path: Path) -> Path:
    executable = tmp_path / "mafft"
    executable.write_text(
        """#!/usr/bin/env python3
import pathlib, sys
if '--version' in sys.argv:
    print('v7.0', file=sys.stderr)
    raise SystemExit(0)
args = sys.argv[1:]
if '--add' in args:
    query = pathlib.Path(args[args.index('--add') + 1]).read_text()
    reference = pathlib.Path(args[-1]).read_text()
    print(reference.rstrip())
    print(query.rstrip())
else:
    print(pathlib.Path(args[-1]).read_text().rstrip())
"""
    )
    executable.chmod(0o755)
    return executable


def _fake_raxml_ng(tmp_path: Path) -> Path:
    executable = tmp_path / "raxml-ng"
    executable.write_text(
        """#!/usr/bin/env python3
import pathlib, sys
if '--version' in sys.argv:
    print('RAxML-NG v1.2.2')
    raise SystemExit(0)
args = sys.argv[1:]
def value(flag): return args[args.index(flag) + 1]
records = []
name = None
parts = []
for line in pathlib.Path(value('--msa')).read_text().splitlines():
    if line.startswith('>'):
        if name is not None: records.append((name, ''.join(parts)))
        name = line[1:].split()[0]; parts = []
    else: parts.append(line.strip())
if name is not None: records.append((name, ''.join(parts)))
query = next((sequence for name, sequence in records if name.startswith('QUERY__')), None)
branches = []
for name, sequence in records:
    if query is None or name.startswith('QUERY__'):
        length = 0.0
    else:
        length = sum(a != b for a, b in zip(sequence, query)) / max(len(query), 1)
    branches.append(f'{name}:{length:.8f}')
pathlib.Path(value('--prefix') + '.raxml.bestTree').write_text('(' + ','.join(branches) + ');\\n')
pathlib.Path(value('--prefix') + '.raxml.bestModel').write_text('GTR{1/1/1/1/1/1}+FU{0.25/0.25/0.25/0.25}+G4{1.0}\\n')
"""
    )
    executable.chmod(0o755)
    return executable


def _fake_epa_ng(tmp_path: Path) -> Path:
    executable = tmp_path / "epa-ng"
    executable.write_text(
        """#!/usr/bin/env python3
import json, pathlib, sys
if '--version' in sys.argv:
    print('EPA-ng v0.3.8')
    raise SystemExit(0)
args = sys.argv[1:]
def value(flag): return args[args.index(flag) + 1]
tree = '(R1:0.25{0},R2:1.0{1}){2};'
document = {
    'tree': tree,
    'placements': [{'p': [[0, -1.0, 1.0, 0.0, 0.0]], 'n': ['QUERY__sample']}],
    'fields': ['edge_num', 'likelihood', 'like_weight_ratio', 'distal_length', 'pendant_length'],
    'version': 3,
}
pathlib.Path(value('--outdir'), 'epa_result.jplace').write_text(json.dumps(document))
"""
    )
    executable.chmod(0o755)
    return executable


def _fake_skani(tmp_path: Path) -> Path:
    executable = tmp_path / "skani"
    executable.write_text(
        """#!/usr/bin/env python3
import pathlib, sys
if '--version' in sys.argv:
    print('skani 0.3.1')
    raise SystemExit(0)
args = sys.argv[1:]
if args[0] == 'sketch':
    output = pathlib.Path(args[args.index('-o') + 1])
    output.mkdir(parents=True)
    references = args[1:args.index('-o')]
    (output / 'fake_references.tsv').write_text('\\n'.join(references) + '\\n')
elif args[0] == 'search':
    database = pathlib.Path(args[args.index('-d') + 1])
    output = pathlib.Path(args[args.index('-o') + 1])
    query = args[args.index('-d') + 2]
    rows = ['Ref_file\\tQuery_file\\tANI\\tAlign_fraction_ref\\tAlign_fraction_query\\tRef_name\\tQuery_name']
    for reference in (database / 'fake_references.tsv').read_text().splitlines():
        ani = '100.00' if pathlib.Path(reference).stem == 'R3' else '99.90'
        rows.append(f'{reference}\\t{query}\\t{ani}\\t99.00\\t98.00\\tref\\tquery')
    output.write_text('\\n'.join(rows) + '\\n')
"""
    )
    executable.chmod(0o755)
    return executable


def _fake_epa_ng_prefers_r2(tmp_path: Path) -> Path:
    executable = tmp_path / "epa-ng-prefers-r2"
    executable.write_text(
        """#!/usr/bin/env python3
import json, pathlib, sys
if '--version' in sys.argv:
    print('EPA-ng v0.3.8')
    raise SystemExit(0)
args = sys.argv[1:]
def value(flag): return args[args.index(flag) + 1]
document = {
    'tree': '(R1:0.25{0},R2:1.0{1}){2};',
    'placements': [{'p': [[1, -1.0, 1.0, 0.0, 0.0]], 'n': ['QUERY__sample']}],
    'fields': ['edge_num', 'likelihood', 'like_weight_ratio', 'distal_length', 'pendant_length'],
    'version': 3,
}
pathlib.Path(value('--outdir'), 'epa_result.jplace').write_text(json.dumps(document))
"""
    )
    executable.chmod(0o755)
    return executable


def _fake_epa_ng_uses_supplied_tree(tmp_path: Path) -> Path:
    executable = tmp_path / "epa-ng-supplied-tree"
    executable.write_text(
        """#!/usr/bin/env python3
import json, pathlib, re, sys
if '--version' in sys.argv:
    print('EPA-ng v0.3.8')
    raise SystemExit(0)
args = sys.argv[1:]
def value(flag): return args[args.index(flag) + 1]
edge = -1
def annotate(match):
    global edge
    edge += 1
    return match.group(0) + '{' + str(edge) + '}'
tree = re.sub(r':[0-9.eE+-]+', annotate, pathlib.Path(value('--tree')).read_text().strip())
document = {
    'tree': tree,
    'placements': [{'p': [[0, -1.0, 1.0, 0.0, 0.0]], 'n': ['QUERY__sample']}],
    'fields': ['edge_num', 'likelihood', 'like_weight_ratio', 'distal_length', 'pendant_length'],
    'version': 3,
}
pathlib.Path(value('--outdir'), 'epa_result.jplace').write_text(json.dumps(document))
"""
    )
    executable.chmod(0o755)
    return executable


def test_sequence_database_directory_uses_locus_filenames(tmp_path):
    database = tmp_path / "database"
    database.mkdir()
    (database / "L1.fasta").write_text(">R1\nAAAA\n>R2\nAAAT\n")
    (database / "unrelated.fasta").write_text(">R1\nCCCC\n")
    assert read_sequence_database(database, {"L1"}) == {
        "L1": [("R1", "AAAA"), ("R2", "AAAT")]
    }


def test_mafft_commands_keep_reference_coordinates():
    assert build_mafft_reference_command("refs.fa", 4, "mafft-x") == [
        "mafft-x", "--auto", "--thread", "4", "refs.fa"
    ]
    assert build_mafft_add_command("query.fa", "refs.aln.fa", 4, "mafft-x") == [
        "mafft-x", "--add", "query.fa", "--keeplength", "--thread", "4", "refs.aln.fa"
    ]
    assert build_raxml_ng_command("placed.fa", "run", 4, "raxml-ng-x") == [
        "raxml-ng-x", "--search", "--msa", "placed.fa", "--model", "DNA",
        "--prefix", "run", "--seed", "12345", "--threads", "4", "--redo",
    ]
    assert build_epa_ng_command("refs.fa", "tree.nwk", "query.fa", "model", "epa", 4, "epa-x") == [
        "epa-x", "--ref-msa", "refs.fa", "--tree", "tree.nwk", "--query", "query.fa",
        "--model", "model", "--outdir", "epa", "--threads", "4",
    ]


def test_raxml_retries_low_pattern_alignment_with_fewer_threads(tmp_path):
    executable = tmp_path / "raxml-ng-retry"
    executable.write_text(
        """#!/usr/bin/env python3
import pathlib, sys
args = sys.argv[1:]
def value(flag): return args[args.index(flag) + 1]
prefix = value('--prefix')
threads = int(value('--threads'))
with open(prefix + '.attempts', 'a') as handle: handle.write(f'{threads}\\n')
if threads > 2:
    print('ERROR: Too few patterns per thread!', file=sys.stderr)
    raise SystemExit(1)
pathlib.Path(prefix + '.raxml.bestTree').write_text('(R1:0,R2:0);\\n')
pathlib.Path(prefix + '.raxml.bestModel').write_text('GTR+G\\n')
"""
    )
    executable.chmod(0o755)
    prefix = tmp_path / "reference"
    tree = tmp_path / "reference_tree.nwk"
    stream = io.StringIO()
    progress = ProgressReporter(stream=stream)

    _run_raxml_ng(
        build_raxml_ng_command("alignment.fa", prefix, 8, str(executable)),
        prefix,
        tree,
        "test tree",
        progress,
    )

    assert Path(f"{prefix}.attempts").read_text().splitlines() == ["8", "4", "2"]
    assert tree.read_text() == "(R1:0,R2:0);\n"
    assert "retrying test tree with 4" in stream.getvalue()
    assert "Too few patterns per thread" in Path(
        f"{prefix}.mlvamaps.raxml.log"
    ).read_text()


def test_phylogenetic_placement_ranks_all_locus_distance_sum(tmp_path):
    database = tmp_path / "database"
    database.mkdir()
    (database / "L1.fasta").write_text(">R1\nAAAA\n>R2\nTTTT\n")
    (database / "L2.fasta").write_text(">R1\nCCCC\n>R2\nGGGG\n")
    result = run_phylogenetic_placement(
        {"L1": "AAAA", "L2": "CCCC"},
        database,
        tmp_path / "out",
        "sample",
        {"L1", "L2"},
        2,
        str(_fake_mafft(tmp_path)),
        str(_fake_raxml_ng(tmp_path)),
        str(_fake_epa_ng(tmp_path)),
        exact_match_fast_path=False,
    )
    matches = _read_tsv(result["phylogenetic_matches"])
    assert matches[0]["reference_id"] == "R1"
    assert matches[0]["compared_loci"] == "2"
    assert matches[0]["total_likelihood_weighted_distance"] == matches[0]["total_phylogenetic_distance"]
    assert float(matches[0]["distance_gap_to_next"]) > 0
    assert (result["phylogeny"] / "L1" / "reference_tree.nwk").exists()
    assert (result["phylogeny"] / "L1" / "epa-ng" / "epa_result.jplace").exists()


@pytest.mark.parametrize("use_build_root", [True, False])
def test_phylogenetic_placement_reuses_reference_build_trees(
    tmp_path, use_build_root
):
    reference_build = tmp_path / "reference_build"
    database = reference_build / "database"
    reference_locus = reference_build / "phylogeny" / "L1"
    database.mkdir(parents=True)
    reference_locus.mkdir(parents=True)
    references = ">R1\nAAAA\n>R2\nTTTT\n"
    (database / "L1.fasta").write_text(references)
    (reference_locus / "references.aligned.fasta").write_text(references)
    saved_tree = "(R1:0.25,R2:1.0);\n"
    (reference_locus / "reference_tree.nwk").write_text(saved_tree)
    (reference_locus / "reference.raxml.bestModel").write_text("GTR+G\n")

    result = run_phylogenetic_placement(
        {"L1": "AAAA"},
        reference_build if use_build_root else database,
        tmp_path / ("root-out" if use_build_root else "database-out"),
        "sample",
        {"L1"},
        2,
        str(_fake_mafft(tmp_path)),
        str(tmp_path / "raxml-must-not-run"),
        str(_fake_epa_ng(tmp_path)),
    )

    locus_output = result["phylogeny"] / "L1"
    assert (locus_output / "reference_tree.nwk").read_text() == saved_tree
    assert not (locus_output / "reference.raxml.bestTree").exists()
    assert _read_tsv(result["phylogenetic_status"])[0]["status"] == "PLACED"


def test_exact_reference_fast_path_skips_external_placement_tools(
    tmp_path, monkeypatch
):
    locus = Locus(locus_id="L1")
    database = tmp_path / "database"
    database.mkdir()
    (database / "L1.fasta").write_text(
        ">R1\nAAAA\n>R2\nTTTT\n>R3\nAAAA\n"
    )
    index_rows = phylogeny_module.reference_sequence_index_rows(
        phylogeny_module.read_sequence_database(database, {"L1"}),
        [locus],
        database,
    )
    phylogeny_module._write_tsv(
        index_rows,
        database / "reference_sequence_index.tsv",
        phylogeny_module.REFERENCE_SEQUENCE_INDEX_FIELDS,
    )
    monkeypatch.setattr(
        phylogeny_module,
        "read_sequence_database",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("indexed exact lookup must not load the FASTA database")
        ),
    )

    result = run_phylogenetic_placement(
        {"L1": "AAAA"},
        database,
        tmp_path / "out",
        "sample",
        [locus],
        1,
        mafft_bin=str(tmp_path / "mafft-must-not-run"),
        raxml_ng_bin=str(tmp_path / "raxml-must-not-run"),
        epa_ng_bin=str(tmp_path / "epa-must-not-run"),
    )

    combined = _read_tsv(result["combined_marker_matches"])
    assert [row["reference_id"] for row in combined] == ["R1", "R3"]
    assert {row["rank"] for row in combined} == {"1"}
    assert {row["match_status"] for row in combined} == {"EXACT_AMPLICON_MATCH"}
    assert {row["combined_marker_distance"] for row in combined} == {"0.00000000"}
    assert {row["tie_break_status"] for row in combined} == {"NOT_APPLICABLE"}
    assert _read_tsv(result["phylogenetic_status"])[0]["status"] == "EXACT_AMPLICON_MATCH"
    assert not (result["phylogeny"] / "L1" / "epa-ng").exists()


def test_exact_reference_label_requires_every_panel_locus():
    loci = [Locus(locus_id="L1"), Locus(locus_id="L2")]
    references = {
        "L1": [("R1", "AAAA")],
        "L2": [("R1", "CCCC")],
    }
    index_rows = phylogeny_module.reference_sequence_index_rows(references, loci)

    match_type, reference_ids, locus_ids, _components = (
        phylogeny_module._exact_reference_group(
            {"L1": "AAAA"},
            {locus.locus_id: locus for locus in loci},
            index_rows,
        )
    )
    assert match_type == ""
    assert reference_ids == []
    assert locus_ids == []


def test_exact_amplicon_ties_use_skani_whole_genome_ani(tmp_path):
    locus = Locus(locus_id="L1")
    build = tmp_path / "reference"
    database = build / "database"
    database.mkdir(parents=True)
    (database / "L1.fasta").write_text(">R1\nAAAA\n>R3\nAAAA\n")
    references = phylogeny_module.read_sequence_database(database, {"L1"})
    index_rows = phylogeny_module.reference_sequence_index_rows(
        references, [locus], database
    )
    phylogeny_module._write_tsv(
        index_rows,
        database / "reference_sequence_index.tsv",
        phylogeny_module.REFERENCE_SEQUENCE_INDEX_FIELDS,
    )
    r1 = tmp_path / "R1.fa"
    r3 = tmp_path / "R3.fa"
    query = tmp_path / "query.fa"
    r1.write_text(">R1\nAAAAT\n")
    r3.write_text(">R3\nAAAAA\n")
    query.write_text(">query\nAAAAA\n")
    phylogeny_module._write_tsv(
        [
            {"reference_id": "R1", "assembly_file": str(r1.resolve())},
            {"reference_id": "R3", "assembly_file": str(r3.resolve())},
        ],
        database / "reference_assemblies.tsv",
        phylogeny_module.REFERENCE_ASSEMBLY_FIELDS,
    )
    fake_skani = _fake_skani(tmp_path)
    phylogeny_module.build_skani_reference_database(
        [("R1", r1), ("R3", r3)], build / "skani", 1, str(fake_skani)
    )

    result = run_phylogenetic_placement(
        {"L1": "AAAA"},
        database,
        tmp_path / "out",
        "sample",
        [locus],
        1,
        query_assembly_path=query,
        skani_bin=str(fake_skani),
    )

    combined = _read_tsv(result["combined_marker_matches"])
    assert [row["reference_id"] for row in combined] == ["R3", "R1"]
    assert [row["rank"] for row in combined] == ["1", "2"]
    assert [row["whole_genome_ani"] for row in combined] == [
        "100.00000000",
        "99.90000000",
    ]
    assert {row["combined_marker_distance"] for row in combined} == {"0.00000000"}
    assert {row["tie_break_status"] for row in combined} == {"APPLIED"}
    assert result["skani_tie_break"].exists()


def test_reference_phylogeny_build_persists_sequence_identity_index(tmp_path):
    reference_build = tmp_path / "reference_build"
    database = reference_build / "database"
    database.mkdir(parents=True)
    (database / "L1.fasta").write_text(
        ">R1\nAAAA\n>R2\nAAAA\n>R3\nAAAT\n"
    )

    result = phylogeny_module.build_reference_phylogenies(
        database,
        reference_build / "phylogeny",
        [Locus(locus_id="L1")],
        1,
        min_references=2,
        mafft_bin=str(_fake_mafft(tmp_path)),
        raxml_ng_bin=str(_fake_raxml_ng(tmp_path)),
    )

    assert result["sequence_index"] == database / "reference_sequence_index.tsv"
    index_rows = _read_tsv(result["sequence_index"])
    assert [row["reference_id"] for row in index_rows] == ["R1", "R2", "R3"]
    assert all(len(row["amplicon_sha256"]) == 64 for row in index_rows)
    assert index_rows[0]["amplicon_sha256"] == index_rows[1]["amplicon_sha256"]
    assert index_rows[0]["amplicon_sha256"] != index_rows[2]["amplicon_sha256"]
    assert [
        name
        for name, _sequence in phylogeny_module._read_fasta(
            result["phylogeny"] / "L1" / "references.fasta"
        )
    ] == ["R1", "R3"]
    haplotypes = _read_tsv(result["haplotype_groups"])
    assert [(row["reference_id"], row["haplotype_id"]) for row in haplotypes] == [
        ("R1", "R1"),
        ("R2", "R1"),
        ("R3", "R3"),
    ]

    placement = run_phylogenetic_placement(
        {"L1": "AATT"},
        reference_build,
        tmp_path / "out",
        "sample",
        [Locus(locus_id="L1")],
        1,
        mafft_bin=str(_fake_mafft(tmp_path)),
        raxml_ng_bin=str(tmp_path / "raxml-must-not-run"),
        epa_ng_bin=str(_fake_epa_ng_uses_supplied_tree(tmp_path)),
    )
    distances = {
        row["reference_id"]: row
        for row in _read_tsv(placement["phylogenetic_distances"])
    }
    assert set(distances) == {"R1", "R2", "R3"}
    assert distances["R1"]["phylogenetic_distance"] == distances["R2"][
        "phylogenetic_distance"
    ]


def test_phylogenetic_placement_parallelizes_across_loci(tmp_path, monkeypatch):
    database = tmp_path / "database"
    database.mkdir()
    (database / "L1.fasta").write_text(">R1\nAAAA\n>R2\nTTTT\n")
    (database / "L2.fasta").write_text(">R1\nCCCC\n>R2\nGGGG\n")
    barrier = Barrier(2)
    observed: list[tuple[str, int]] = []
    original = phylogeny_module._run_placement_job

    def synchronized_job(job, mafft, epa_ng, native_threads=1):
        observed.append((job.locus_id, native_threads))
        barrier.wait(timeout=2)
        return original(job, mafft, epa_ng, native_threads)

    monkeypatch.setattr(phylogeny_module, "_run_placement_job", synchronized_job)
    stream = io.StringIO()
    result = run_phylogenetic_placement(
        {"L1": "AAAA", "L2": "CCCC"},
        database,
        tmp_path / "out",
        "sample",
        {"L1", "L2"},
        2,
        str(_fake_mafft(tmp_path)),
        str(_fake_raxml_ng(tmp_path)),
        str(_fake_epa_ng(tmp_path)),
        progress=ProgressReporter(stream=stream),
    )

    assert sorted(observed) == [("L1", 1), ("L2", 1)]
    assert _read_tsv(result["phylogenetic_status"])[0]["status"] == "PLACED"
    assert "with 2 worker(s)" in stream.getvalue()
    assert "Completed EPA-ng loci: 2/2 (100.0%)" in stream.getvalue()
    assert "Computing reference-tip distances" in stream.getvalue()
    assert "Building combined neighbor-joining tree" in stream.getvalue()
    assert "Finished phylogenetic placement summaries" in stream.getvalue()


def test_missing_query_still_builds_reference_tree(tmp_path):
    database = tmp_path / "database"
    database.mkdir()
    (database / "L1.fasta").write_text(">R1\nAAAA\n>R2\nAAAT\n")
    result = run_phylogenetic_placement(
        {}, database, tmp_path / "out", "sample", {"L1"}, 1,
        str(_fake_mafft(tmp_path)), str(_fake_raxml_ng(tmp_path)),
        str(_fake_epa_ng(tmp_path)),
    )
    assert (result["phylogeny"] / "L1" / "reference_tree.nwk").exists()
    assert _read_tsv(result["phylogenetic_status"])[0]["status"] == "REFERENCE_TREE_ONLY"


def test_jplace_edge_lengths_drive_patristic_distances():
    root = _parse_newick("((R1:0.1{0},R2:0.2{1}):0.3{2},'R 3':0.4{3}){4};\n")
    distances = _placement_patristic_distances(
        root, placement_edge=2, distal_length=0.1, pendant_length=0.05
    )
    assert distances["R1"] == pytest.approx(0.25)
    assert distances["R2"] == pytest.approx(0.35)
    assert distances["R 3"] == pytest.approx(0.65)


def test_all_tip_distances_use_tree_branch_lengths():
    root = _parse_newick("((A:1,B:2):3,C:4);")
    assert _tip_patristic_distances(root) == {
        ("A", "B"): pytest.approx(3.0),
        ("A", "C"): pytest.approx(8.0),
        ("B", "C"): pytest.approx(9.0),
    }


def test_vectorized_neighbor_joining_is_deterministic():
    labels = ["A", "B", "C", "D"]
    distances = {
        ("A", "B"): 5.0,
        ("A", "C"): 9.0,
        ("A", "D"): 9.0,
        ("B", "C"): 10.0,
        ("B", "D"): 10.0,
        ("C", "D"): 8.0,
    }
    first = neighbor_joining_tree(labels, distances)
    assert neighbor_joining_tree(labels, distances) == first
    assert set(_tip_patristic_distances(_parse_newick(first))) == set(distances)


def test_epa_parser_uses_highest_likelihood_weight_placement(tmp_path):
    jplace = tmp_path / "result.jplace"
    jplace.write_text(
        json.dumps(
            {
                "tree": "(R1:0.25{0},R2:1.0{1}){2};",
                "fields": [
                    "edge_num", "likelihood", "like_weight_ratio",
                    "distal_length", "pendant_length",
                ],
                "placements": [
                    {"p": [[1, -2.0, 0.1, 0.0, 0.1], [0, -1.0, 0.9, 0.0, 0.05]]}
                ],
                "version": 3,
            }
        )
    )
    distances, metadata = read_epa_ng_placement(jplace)
    assert metadata["placement_edge"] == 0
    assert metadata["like_weight_ratio"] == pytest.approx(0.9)
    assert distances["R1"] == pytest.approx(0.05)
    best, expected, statistics = read_epa_ng_placement_statistics(jplace)
    assert best == distances
    assert expected["R1"] == pytest.approx(0.18)
    assert expected["R2"] == pytest.approx(1.18)
    assert statistics["placement_count"] == 2
    assert statistics["placement_entropy"] > 0


def test_repeat_masking_and_combined_marker_ranking(tmp_path):
    locus = Locus(
        locus_id="L1",
        forward_primer="ACG",
        reverse_primer="TTA",
        left_flank_sequence="TT",
        right_flank_sequence="CC",
        repeat_motif="GA",
        repeat_unit_length_bp=2,
        expected_min_repeats=1,
        expected_max_repeats=6,
    )
    query = "ACGTTGAGACCTAA"
    components = decompose_marker_sequence(locus, query)
    assert components.repeat_sequence == "GAGA"
    assert components.repeat_count == 2
    assert components.snp_sequence == "ACGTTCCTAA"
    assert components.masking_method == "flank_bounded"

    database = tmp_path / "database"
    database.mkdir()
    (database / "L1.fasta").write_text(
        ">R1\nACGTTGAGACCTAA\n>R2\nACGTTGAGAGACCTAA\n"
    )
    (database / "reference_metadata.tsv").write_text(
        "reference_id\tcollection_date\tlatitude\tlongitude\tlocation\tsource\n"
        "R1\t2024-01-02\t40.0\t-75.0\tSite A\tenvironment\n"
        "R2\t2023-05-01\t41.0\t-74.0\tSite B\tclinical\n"
    )
    result = run_phylogenetic_placement(
        {"L1": query},
        database,
        tmp_path / "out",
        "sample",
        [locus],
        1,
        str(_fake_mafft(tmp_path)),
        str(_fake_raxml_ng(tmp_path)),
        str(_fake_epa_ng(tmp_path)),
        exact_match_fast_path=False,
    )
    marker_rows = _read_tsv(result["marker_components"])
    query_row = next(row for row in marker_rows if row["record_type"] == "query")
    assert query_row["repeat_haplotype"] == "GA|GA"
    assert query_row["snp_sequence_length"] == "10"
    locus_rows = {
        row["reference_id"]: row
        for row in _read_tsv(result["locus_marker_distances"])
    }
    assert locus_rows["R1"]["repeat_count_delta"] == "0.000000"
    assert locus_rows["R2"]["repeat_count_delta"] == "1.000000"
    combined = _read_tsv(result["combined_marker_matches"])
    assert combined[0]["reference_id"] == "R1"
    assert combined[0]["repeat_compared_loci"] == "1"
    assert combined[0]["collection_date"] == "2024-01-02"
    assert combined[0]["location"] == "Site A"
    closest_bands = _read_tsv(result["closest_reference_bands"])
    assert closest_bands == [
        {
            "reference_id": "R1",
            "locus_id": "L1",
            "product_size_bp": str(len(query)),
            "repeat_count": "2",
        }
    ]
    combined_tree = result["combined_marker_tree"]
    assert combined_tree.name == "combined_markers.tree"
    root = _parse_newick(combined_tree.read_text())

    def tips(node):
        if not node.children:
            return {node.name}
        return set().union(*(tips(child) for child, _length in node.children))

    assert tips(root) == {"R1", "R2", "sample"}


def test_exact_marker_match_overrides_conflicting_epa_ranking(tmp_path):
    locus = Locus(
        locus_id="L1",
        forward_primer="ACG",
        reverse_primer="TTA",
        left_flank_sequence="TT",
        right_flank_sequence="CC",
        repeat_motif="GA",
        repeat_unit_length_bp=2,
        expected_min_repeats=1,
        expected_max_repeats=6,
    )
    query = "ACGTTGAGACCTAA"
    database = tmp_path / "database"
    database.mkdir()
    (database / "L1.fasta").write_text(
        f">R1\n{query}\n>R2\nCCGTTGAGACCTAA\n"
    )

    result = run_phylogenetic_placement(
        {"L1": query},
        database,
        tmp_path / "out",
        "sample",
        [locus],
        1,
        str(_fake_mafft(tmp_path)),
        str(_fake_raxml_ng(tmp_path)),
        str(_fake_epa_ng_prefers_r2(tmp_path)),
        exact_match_fast_path=False,
    )

    locus_rows = {
        row["reference_id"]: row
        for row in _read_tsv(result["locus_marker_distances"])
    }
    assert locus_rows["R1"]["exact_snp_match"] == "yes"
    assert locus_rows["R1"]["direct_snp_distance"] == "0.00000000"
    assert float(locus_rows["R1"]["placement_normalized_snp_distance"]) > 0
    assert locus_rows["R1"]["normalized_snp_distance"] == "0.00000000"

    combined = _read_tsv(result["combined_marker_matches"])
    assert combined[0]["reference_id"] == "R1"
    assert combined[0]["rank"] == "1"
    assert combined[0]["combined_marker_distance"] == "0.00000000"
    assert combined[0]["match_status"] == "EXACT_MARKER_MATCH"
    assert combined[0]["ranking_warning"] == "EXACT_MATCH_OVERRIDES_PLACEMENT"
    assert combined[1]["reference_id"] == "R2"
