import pytest

from guidance3.trees import cli

TREES = "((A,B),(C,D),(E,F));\n" * 3 + "((A,C),(B,D),(E,F));\n"


def _alt_trees(tmp_path):
    p = tmp_path / "alt.nwk"
    p.write_text(TREES)
    return str(p)


def test_user_tree_from_alt_trees_without_base_msa(tmp_path):
    user = tmp_path / "user.nwk"
    user.write_text("((A,B),(C,D),(E,F));\n")
    out = tmp_path / "out"
    rc = cli.main(["--alt-trees", _alt_trees(tmp_path), "--no-consensus", "--user-tree", str(user),
                   "--out-dir", str(out)])
    assert rc == 0
    lines = [l for l in (out / "MSA_UserTree_alignment_support.tsv").read_text().splitlines()
             if not l.startswith("#")]
    rows = [dict(zip(lines[0].split("\t"), l.split("\t"))) for l in lines[1:]]
    assert {r["taxa_smaller_side"]: int(r["alignment_support"]) for r in rows} == {"A,B": 75, "C,D": 75, "E,F": 100}
    assert all(r["abayes"] == "NA" for r in rows)
    assert not any(p.name.startswith("guidance3_tree_") for p in out.iterdir())  # temp dir removed


@pytest.mark.parametrize("extra, message", [
    (["--no-consensus"], "needs --user-tree"),
    ([], "needs --base-msa"),
    (["--n-trees", "2", "--base-msa", "x.fasta"], "at least 3"),
])
def test_argument_errors(tmp_path, capsys, extra, message):
    with pytest.raises(SystemExit) as e:
        cli.main(["--alt-trees", _alt_trees(tmp_path), "--out-dir", str(tmp_path / "out"), *extra])
    assert message in str(e.value)


def test_alt_trees_with_different_leaves(tmp_path):
    p = tmp_path / "alt.nwk"
    p.write_text("((A,B),(C,D),(E,F));\n" * 2 + "((A,B),(C,D),(E,X));\n")
    user = tmp_path / "user.nwk"
    user.write_text("((A,B),(C,D),(E,F));\n")
    with pytest.raises(SystemExit) as e:
        cli.main(["--alt-trees", str(p), "--no-consensus", "--user-tree", str(user), "--out-dir", str(tmp_path / "o")])
    assert "tree 3" in str(e.value)
