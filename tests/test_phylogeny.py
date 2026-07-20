from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from mlvamaps.phylogeny import (
    _parse_newick,
    _placement_patristic_distances,
    build_mafft_add_command,
    build_mafft_reference_command,
    build_epa_ng_command,
    build_raxml_ng_command,
    read_sequence_database,
    read_epa_ng_placement,
    run_phylogenetic_placement,
)


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
        "raxml-ng-x", "--search", "--msa", "placed.fa", "--model", "GTR+G",
        "--prefix", "run", "--seed", "12345", "--threads", "4", "--redo",
    ]
    assert build_epa_ng_command("refs.fa", "tree.nwk", "query.fa", "model", "epa", 4, "epa-x") == [
        "epa-x", "--ref-msa", "refs.fa", "--tree", "tree.nwk", "--query", "query.fa",
        "--model", "model", "--outdir", "epa", "--threads", "4",
    ]


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
    )
    matches = _read_tsv(result["phylogenetic_matches"])
    assert matches[0]["reference_id"] == "R1"
    assert matches[0]["compared_loci"] == "2"
    assert (result["phylogeny"] / "L1" / "reference_tree.nwk").exists()
    assert (result["phylogeny"] / "L1" / "epa-ng" / "epa_result.jplace").exists()


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
