"""Phylogenetic inference from the GUIDANCE3 alternative MSAs.

One ML tree is inferred from each of k randomly chosen alternative MSAs (IQ-TREE -fast). From these trees:

* the GUIDANCE3 tree: the extended majority-rule consensus (splits in more than half of the trees, then the
  remaining splits from most to least frequent while compatible), resolved by one ML search on the reference MSA
  constrained to it (IQ-TREE -g), which also estimates the branch lengths;
* alignment support of every branch: the percentage of the k alternative-MSA trees that contain it;
* the same alignment support for a user-supplied tree.

Trees are written in Newick with alignment support (0-100) as internal-node labels and the original sequence names.
"""
import io
import os
import random
import shutil
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

from Bio import Phylo, SeqIO

from guidance3.trees.tree_viewer import write_tree_html

DEFAULT_N_ALT_TREES = 20
_SEED = 1
_NUC = set("ACGTUN-?.")


# ------------------------------------------------------------------------------------------------ splits

def leaf_names(tree):
    return sorted(t.name for t in tree.get_terminals())


def tree_splits(tree, leaves):
    """Non-trivial bipartitions of an (unrooted) tree, each as the frozenset of leaves on the side WITHOUT the
    reference leaf leaves[0], so the same split has one representation whatever the rooting."""
    ref, all_leaves, n = leaves[0], frozenset(leaves), len(leaves)
    out = set()
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal() or clade is tree.root:
            continue
        side = frozenset(t.name for t in clade.get_terminals())
        if ref in side:
            side = all_leaves - side
        if 2 <= len(side) <= n - 2:
            out.add(side)
    return out


def _compatible(a, b):
    return not (a & b) or a <= b or b <= a


def extended_majority_splits(split_sets):
    """Extended majority-rule consensus (= IQ-TREE -con default): majority splits, then the remaining splits from
    most to least frequent while compatible with all accepted ones. Returns (accepted splits, Counter)."""
    counts = Counter()
    for s in split_sets:
        counts.update(s)
    half = len(split_sets) / 2
    keep = [s for s, c in counts.items() if c > half]
    for s, _c in sorted(counts.items(), key=lambda x: (-x[1], sorted(x[0]))):
        if s not in keep and all(_compatible(s, k) for k in keep):
            keep.append(s)
    return keep, counts


def _quote(name):
    return f"'{name}'" if any(ch in name for ch in " ():;,[]'\t") else name


def splits_to_newick(splits, leaves):
    """Newick topology (no branch lengths) of compatible splits on `leaves`."""
    root = {"leaves": set(leaves), "children": []}
    for split in sorted(splits, key=len, reverse=True):  # largest first: nest under the smallest containing node
        node = root
        while True:
            nxt = next((c for c in node["children"] if split <= c["leaves"]), None)
            if nxt is None:
                break
            node = nxt
        moved = [c for c in node["children"] if c["leaves"] <= split]
        node["children"] = [c for c in node["children"] if c not in moved]
        node["children"].append({"leaves": set(split), "children": moved})

    def render(node):
        placed = set().union(*(c["leaves"] for c in node["children"])) if node["children"] else set()
        parts = [render(c) for c in node["children"]] + [_quote(l) for l in sorted(node["leaves"] - placed)]
        return "(" + ",".join(parts) + ")"

    return render(root) + ";"


# ------------------------------------------------------------------------------------------------ helpers

def _read_codes(codes_file):
    """Seqs.Codes: '<original name>\\t<code>' per line -> {code: name}."""
    codes = {}
    with open(codes_file) as fh:
        for line in fh:
            if line.strip():
                name, code = line.rstrip("\n").split("\t")
                codes[code] = name
    return codes


def _is_nucleotide(fasta):
    chars = Counter()
    for rec in SeqIO.parse(fasta, "fasta"):
        chars.update(str(rec.seq).upper())
    total = sum(chars.values())
    return total > 0 and sum(c for ch, c in chars.items() if ch in _NUC) / total > 0.9


def _model(fasta):
    return "GTR+G" if _is_nucleotide(fasta) else "LG+G"


def _run(cmd, log_file):
    with open(log_file, "a") as log:
        log.write("alignment_trees: " + " ".join(cmd) + "\n")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"IQ-TREE failed ({res.returncode}): {' '.join(cmd)}\n{res.stderr[-2000:]}")


def _infer(iqtree, msa, prefix, model, log_file, extra=()):
    _run([iqtree, "-s", msa, "-m", model, "-fast", "-T", "1", "-seed", str(_SEED), "-pre", prefix,
          "-redo", "-quiet", *extra], log_file)
    return Phylo.read(f"{prefix}.treefile", "newick")


def _clean_copy(tree):
    """Copy of a tree without internal labels (IQ-TREE -te input)."""
    import copy
    t = copy.deepcopy(tree)
    for c in t.find_clades():
        if not c.is_terminal():
            c.confidence = None
            c.name = None
    return t


def standard_supports(iqtree, msa, tree, model, prefix, log_file):
    """aBayes and SH-aLRT (1000 replicates) of each branch of `tree` (leaves = MSA names), computed on `msa` with
    the topology fixed. Returns {split: (sh_alrt, abayes)}; {} if IQ-TREE cannot compute them (e.g. polytomies)."""
    tree_file = f"{prefix}.input.nwk"
    _write_tree(_clean_copy(tree), tree_file, plain=True)
    try:
        _run([iqtree, "-s", msa, "-m", model, "-te", tree_file, "-abayes", "-alrt", "1000", "-T", "1",
              "-seed", str(_SEED), "-pre", prefix, "-redo", "-quiet"], log_file)
        labelled = Phylo.read(f"{prefix}.treefile", "newick")
    except Exception as e:
        with open(log_file, "a") as log:
            log.write(f"alignment_trees: aBayes / SH-aLRT not computed: {e}\n")
        return {}
    leaves = leaf_names(labelled)
    ref, all_leaves, n = leaves[0], frozenset(leaves), len(leaves)
    out = {}
    for clade in labelled.find_clades():
        if clade.is_terminal() or clade is labelled.root:
            continue
        label = clade.name if clade.name is not None else (str(clade.confidence) if clade.confidence is not None else "")
        parts = label.split("/")
        if len(parts) != 2:
            continue
        side = frozenset(t.name for t in clade.get_terminals())
        split = all_leaves - side if ref in side else side
        try:
            out[split] = (float(parts[0]), float(parts[1]))
        except ValueError:
            continue
    return out


def _annotate(tree, counts, k, codes):
    """Set each internal node's label to its alignment support (%) and rename leaves to original names.
    Returns rows (dicts) for the branch table."""
    leaves = leaf_names(tree)
    ref, all_leaves, n = leaves[0], frozenset(leaves), len(leaves)
    rows, seen = [], set()
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        clade.confidence = None
        clade.name = None
        if clade is tree.root:
            continue
        side = frozenset(t.name for t in clade.get_terminals())
        split = all_leaves - side if ref in side else side
        if not 2 <= len(split) <= n - 2:
            continue
        support = round(100 * counts.get(split, 0) / k)
        clade.confidence = support
        if split in seen:  # the two root edges of a rooted tree are one unrooted branch
            continue
        seen.add(split)
        smaller = split if len(split) <= n / 2 else all_leaves - split
        rows.append({"branch": len(rows) + 1, "split": split, "alignment_support": support,
                     "n_smaller": len(smaller), "taxa": ",".join(sorted(codes.get(x, x) for x in smaller)),
                     "key": _viewer_key(split, all_leaves, codes)})
    for t in tree.get_terminals():
        t.name = codes.get(t.name, t.name)
    return rows


def _viewer_key(split, all_leaves, codes):
    """Branch id shared with the tree viewer: original names of the side without the alphabetically first name."""
    names = {x: codes.get(x, x) for x in all_leaves}
    first = min(names.values())
    side = {names[x] for x in split}
    if first in side:
        side = set(names.values()) - side
    return "|".join(sorted(side))


def _write_tree(tree, path, plain=False):
    buf = io.StringIO()
    if plain:
        Phylo.write(tree, buf, "newick")
    else:
        Phylo.write(tree, buf, "newick", format_confidence="%d")
    with open(path, "w") as fh:
        fh.write(buf.getvalue())


def _fmt(x, digits):
    return "NA" if x is None else f"{x:.{digits}f}"


def _write_table(rows, path, k, std):
    """Branch table: alignment support, aBayes and SH-aLRT on the base MSA, and aBayes x alignment support.
    Returns {viewer key: {...}} for the tree viewer."""
    info = {}
    with open(path, "w") as fh:
        fh.write(f"# alignment_support: % of the {k} trees inferred from alternative MSAs that contain the branch\n")
        fh.write("# abayes, sh_alrt: aBayes and SH-aLRT (1000 replicates) on the base MSA, topology fixed (IQ-TREE)\n")
        fh.write("# combined: abayes x alignment_support / 100 (in our simulations it identified incorrect branches "
                 "better than either measure alone)\n")
        fh.write("branch\talignment_support\tabayes\tsh_alrt\tcombined\tn_taxa_smaller_side\ttaxa_smaller_side\n")
        for r in rows:
            alrt, abayes = std.get(r["split"], (None, None))
            combined = None if abayes is None else abayes * r["alignment_support"] / 100
            fh.write("\t".join([str(r["branch"]), str(r["alignment_support"]), _fmt(abayes, 3), _fmt(alrt, 1),
                                 _fmt(combined, 3), str(r["n_smaller"]), r["taxa"]]) + "\n")
            info[r["key"]] = {"abayes": abayes, "alrt": alrt, "combined": combined}
    return info


# ------------------------------------------------------------------------------------------------ main step

def run_alignment_trees(config):
    """Pipeline step (GUIDANCE3 only; call before create_tar_archives, which removes the alternative MSAs).
    Writes into config.WorkingDir:
      <prefix>_GUIDANCE3_tree.nwk                    resolved extended consensus, alignment support as labels
      <prefix>_GUIDANCE3_tree_branch_support.tsv     one row per branch
      <prefix>_GUIDANCE3_tree.html                   interactive tree, branches coloured by alignment support
      <prefix>_AltMSA_trees.nwk                      the k per-alternative trees
      <prefix>_UserTree_alignment_support.nwk/.tsv   with --userTree
    Returns True on success."""
    work = config.WorkingDir
    prefix = os.path.join(work, config.Output_Prefix)
    log_file = config.OutLogFile
    want_consensus = not getattr(config, "no_consensus_tree", False)
    user_tree_file = getattr(config, "user_tree", "") or ""
    if not want_consensus and not user_tree_file:
        return True

    alt_dir = config.GUIDANCE3_MSAs_Dir
    alts = sorted(f for f in os.listdir(alt_dir) if f.endswith(".fasta")) if os.path.isdir(alt_dir) else []
    if len(alts) < 3:
        print(f"[alignment_trees] only {len(alts)} alternative MSAs found in {alt_dir}, skipping")
        return False
    k = min(int(getattr(config, "n_alt_trees", DEFAULT_N_ALT_TREES) or DEFAULT_N_ALT_TREES), len(alts))
    chosen = sorted(random.Random(_SEED).sample(alts, k))

    tmp = os.path.join(work, "alignment_trees")
    os.makedirs(tmp, exist_ok=True)
    reference = os.path.join(work, config.Alignment_File)
    model = _model(os.path.join(alt_dir, chosen[0]))
    if _is_nucleotide(reference) != (model == "GTR+G"):
        # codon runs: alternatives are amino-acid alignments; translate a copy of the codon reference to match
        from guidance3.sequences.fasta import _translate_codon_msa_inplace
        ref_copy = os.path.join(tmp, "reference.fasta")
        shutil.copy(reference, ref_copy)
        _translate_codon_msa_inplace(ref_copy)
        reference = ref_copy
    codes = _read_codes(os.path.join(work, config.code_fileName))

    print(f"[alignment_trees] inferring {k} trees from alternative MSAs (IQ-TREE {model}, -fast)...")
    iqtree = config.iqtree_prog

    def one(i_name):
        i, name = i_name
        return _infer(iqtree, os.path.join(alt_dir, name), os.path.join(tmp, f"alt{i:03d}"), model, log_file)

    with ThreadPoolExecutor(max_workers=max(1, int(config.proc_num))) as pool:
        trees = list(pool.map(one, enumerate(chosen)))

    leaves = leaf_names(trees[0])
    split_sets = [tree_splits(t, leaves) for t in trees]
    keep, counts = extended_majority_splits(split_sets)

    # the k trees, with original names
    with open(f"{prefix}_AltMSA_trees.nwk", "w") as fh:
        for t in trees:
            for leaf in t.get_terminals():
                leaf.name = codes.get(leaf.name, leaf.name)
            buf = io.StringIO()
            Phylo.write(t, buf, "newick")
            fh.write(buf.getvalue())

    ok = True
    if want_consensus:
        constraint = os.path.join(tmp, "consensus.constraint.nwk")
        with open(constraint, "w") as fh:
            fh.write(splits_to_newick(keep, leaves) + "\n")
        resolved = _infer(iqtree, reference, os.path.join(tmp, "guidance3_tree"), model, log_file,
                          extra=("-g", constraint))
        std = standard_supports(iqtree, reference, resolved, model, os.path.join(tmp, "guidance3_tree_support"),
                                log_file)
        rows = _annotate(resolved, counts, k, codes)
        _write_tree(resolved, f"{prefix}_GUIDANCE3_tree.nwk")
        info = _write_table(rows, f"{prefix}_GUIDANCE3_tree_branch_support.tsv", k, std)
        write_tree_html(f"{prefix}_GUIDANCE3_tree.nwk", f"{prefix}_GUIDANCE3_tree.html", "GUIDANCE3 tree", k, info)
        print(f"[alignment_trees] GUIDANCE3 tree -> {prefix}_GUIDANCE3_tree.nwk")

    if user_tree_file:
        ok = evaluate_user_tree(user_tree_file, counts, k, codes, leaves, prefix, log_file,
                                support_args=(iqtree, reference, model, os.path.join(tmp, "user_tree_support"))) and ok

    shutil.rmtree(tmp, ignore_errors=True)
    return ok


def evaluate_user_tree(user_tree_file, counts, k, codes, leaves, prefix, log_file, support_args=None):
    """Alignment support for each branch of a user tree whose leaves are the original sequence names.
    support_args = (iqtree, base MSA, model, work prefix): also aBayes and SH-aLRT on the base MSA."""
    name_to_code = {name: code for code, name in codes.items()}
    try:
        tree = Phylo.read(user_tree_file, "newick")
    except Exception as e:
        msg = f"Could not read the user tree as Newick: {e}"
        _report_user_tree_error(msg, prefix, log_file)
        return False
    names = [t.name for t in tree.get_terminals()]
    unknown = sorted(n for n in names if n not in name_to_code)
    missing = sorted(set(codes.values()) - set(names))
    if unknown or missing or len(names) != len(set(names)):
        msg = ("The user tree must contain each input sequence exactly once, with the same names."
               + (f" Not among the sequences: {', '.join(unknown[:10])}." if unknown else "")
               + (f" Missing from the tree: {', '.join(missing[:10])}." if missing else ""))
        _report_user_tree_error(msg, prefix, log_file)
        return False
    for t in tree.get_terminals():
        t.name = name_to_code[t.name]
    std = {}
    if support_args:
        iqtree, msa, model, work_prefix = support_args
        std = standard_supports(iqtree, msa, tree, model, work_prefix, log_file)
    rows = _annotate(tree, counts, k, codes)
    _write_tree(tree, f"{prefix}_UserTree_alignment_support.nwk")
    info = _write_table(rows, f"{prefix}_UserTree_alignment_support.tsv", k, std)
    write_tree_html(f"{prefix}_UserTree_alignment_support.nwk", f"{prefix}_UserTree_alignment_support.html",
                    "Your tree: alignment support", k, info)
    print(f"[alignment_trees] user tree with alignment support -> {prefix}_UserTree_alignment_support.nwk")
    return True


def _report_user_tree_error(msg, prefix, log_file):
    print(f"[alignment_trees] {msg}")
    with open(log_file, "a") as log:
        log.write(f"alignment_trees: user tree skipped: {msg}\n")
    with open(f"{prefix}_UserTree_ERROR.txt", "w") as fh:
        fh.write(msg + "\n")
