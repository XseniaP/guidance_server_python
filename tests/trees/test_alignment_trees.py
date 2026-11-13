import io

from Bio import Phylo

from guidance3.trees import alignment_trees as at


def _tree(nwk):
    return Phylo.read(io.StringIO(nwk), "newick")


LEAVES = ["A", "B", "C", "D", "E", "F"]


def test_splits_ignore_rooting():
    t1 = _tree("((A,B),(C,D),(E,F));")
    t2 = _tree("(A,(B,((C,D),(E,F))));")  # same unrooted topology, different root
    assert at.tree_splits(t1, LEAVES) == at.tree_splits(t2, LEAVES)
    assert at.tree_splits(t1, LEAVES) == {frozenset("CD"), frozenset("EF"), frozenset("CDEF")}


def test_extended_majority_adds_compatible_minority_split():
    trees = [_tree("((A,B),(C,D),(E,F));")] * 3 + [_tree("((A,C),(B,D),(E,F));")] * 2
    keep, counts = at.extended_majority_splits([at.tree_splits(t, LEAVES) for t in trees])
    assert counts[frozenset("EF")] == 5
    assert frozenset("CD") in keep and frozenset("EF") in keep
    assert frozenset("BD") not in keep  # incompatible with the majority split CD


def test_splits_to_newick_roundtrip():
    splits = {frozenset("CD"), frozenset("EF"), frozenset("CDEF")}
    t = _tree(at.splits_to_newick(splits, LEAVES))
    assert at.tree_splits(t, LEAVES) == splits


def test_splits_to_newick_leaves_polytomy_for_missing_splits():
    t = _tree(at.splits_to_newick({frozenset("EF")}, LEAVES))
    assert at.tree_splits(t, LEAVES) == {frozenset("EF")}
    assert sorted(x.name for x in t.get_terminals()) == LEAVES


def _counts(trees):
    _keep, counts = at.extended_majority_splits([at.tree_splits(t, LEAVES) for t in trees])
    return counts


def test_user_tree_alignment_support(tmp_path):
    trees = [_tree("((A,B),(C,D),(E,F));")] * 3 + [_tree("((A,C),(B,D),(E,F));")]
    user = tmp_path / "user.nwk"
    user.write_text("((A,B),(C,D),(E,F));\n")
    codes = {x: x for x in LEAVES}
    prefix = str(tmp_path / "out")
    assert at.evaluate_user_tree(str(user), _counts(trees), 4, codes, LEAVES, prefix, str(tmp_path / "log"))
    lines = [l for l in open(prefix + "_UserTree_alignment_support.tsv").read().splitlines() if not l.startswith("#")]
    header, rows = lines[0].split("\t"), [dict(zip(lines[0].split("\t"), l.split("\t"))) for l in lines[1:]]
    assert header[:5] == ["branch", "alignment_support", "abayes", "sh_alrt", "combined"]
    support = {r["taxa_smaller_side"]: int(r["alignment_support"]) for r in rows}
    assert support == {"A,B": 75, "C,D": 75, "E,F": 100}
    assert all(r["abayes"] == "NA" for r in rows)  # no IQ-TREE run in this test


def test_user_tree_with_wrong_names_is_reported(tmp_path):
    trees = [_tree("((A,B),(C,D),(E,F));")]
    user = tmp_path / "user.nwk"
    user.write_text("((A,B),(C,D),(E,X));\n")
    prefix = str(tmp_path / "out")
    ok = at.evaluate_user_tree(str(user), _counts(trees), 1, {x: x for x in LEAVES}, LEAVES, prefix,
                               str(tmp_path / "log"))
    assert not ok
    msg = open(prefix + "_UserTree_ERROR.txt").read()
    assert "X" in msg and "F" in msg
