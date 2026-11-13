"""guidance3-tree: the GUIDANCE3 tree and alignment support from any set of alternative MSAs.

The same computation as the tree step of the guidance3 pipeline, for MSAs made elsewhere (other aligners,
other settings, the alternative MSAs of an earlier run):

* one ML tree (IQ-TREE 3, -fast; LG+G for proteins, GTR+G for nucleotides) from each of --n-trees randomly
  chosen MSAs of the set;
* the GUIDANCE3 tree: their extended majority-rule consensus, resolved by an ML search on the base MSA
  constrained to it;
* for every branch: alignment support (% of the trees that contain it), and aBayes and SH-aLRT on the base MSA;
* the same table for a tree of your own (--user-tree).

--alt-trees reuses trees already inferred (e.g. <prefix>_AltMSA_trees.nwk of a guidance3 run) instead of
inferring them, so evaluating another tree takes seconds.
"""
import argparse
import os
import shutil
import sys
import tempfile

from Bio import Phylo

from guidance3.trees import alignment_trees as at


def _write_coded(src, dst, name_to_code):
    from guidance3.sequences.msa_set import read_fasta
    with open(dst, "w") as fh:
        for name, seq in read_fasta(src):
            fh.write(f">{name_to_code[name]}\n{seq}\n")


def _fail(msg):
    sys.exit(f"guidance3-tree: error: {msg}")


def build_parser():
    p = argparse.ArgumentParser(prog="guidance3-tree", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--msas-dir", help="Folder, or .zip / .tar.gz archive, of aligned FASTA files of the same "
                                        "sequences (at least 3).")
    src.add_argument("--alt-trees", help="Newick file with one tree per line, inferred from alternative MSAs "
                                         "(e.g. <prefix>_AltMSA_trees.nwk of a guidance3 run).")
    p.add_argument("--base-msa", help="MSA used to resolve the consensus and to compute aBayes and SH-aLRT. "
                                      "Default with --msas-dir: the MSA of the set that the GUIDANCE3 DL model "
                                      "ranks best. Required with --alt-trees, except with --no-consensus.")
    p.add_argument("--user-tree", help="Newick tree of the same sequence names: report its branch supports.")
    p.add_argument("--out-dir", required=True, help="Output directory.")
    p.add_argument("--prefix", default="MSA", help="Prefix of the output files. Default: MSA.")
    p.add_argument("--n-trees", type=int, default=at.DEFAULT_N_ALT_TREES,
                   help=f"Number of alternative MSAs (chosen at random, seeded) to infer trees from (>= 3). "
                        f"Default: {at.DEFAULT_N_ALT_TREES}.")
    p.add_argument("--no-consensus", action="store_true",
                   help="Do not build the GUIDANCE3 tree (only evaluate --user-tree).")
    p.add_argument("--proc-num", type=int, default=2, help="IQ-TREE runs in parallel. Default: 2.")
    p.add_argument("--keep-temp", action="store_true", help="Keep the intermediate IQ-TREE files.")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    from guidance3.constants import IQTREE
    from guidance3.sequences.msa_set import (MSASetError, collect_msa_files, read_fasta, validate_msa_set,
                                             rank_msas)

    if args.no_consensus and not args.user_tree:
        _fail("--no-consensus needs --user-tree (otherwise there is nothing to compute).")
    if args.n_trees < 3:
        _fail("--n-trees must be at least 3.")
    if args.alt_trees and not args.base_msa and not args.no_consensus:
        _fail("--alt-trees needs --base-msa to build the GUIDANCE3 tree (or use --no-consensus with --user-tree).")
    for f in (args.alt_trees, args.base_msa, args.user_tree):
        if f and not os.path.isfile(f):
            _fail(f"file not found: {f}")
    if args.msas_dir and not os.path.exists(args.msas_dir):
        _fail(f"not found: {args.msas_dir}")

    os.makedirs(args.out_dir, exist_ok=True)
    prefix = os.path.join(args.out_dir, args.prefix)
    log_file = f"{prefix}_guidance3_tree.log"
    open(log_file, "w").close()
    tmp = tempfile.mkdtemp(prefix="guidance3_tree_", dir=args.out_dir)
    try:
        alt_files, alt_trees = [], None
        if args.msas_dir:
            try:
                files = collect_msa_files(args.msas_dir, os.path.join(tmp, "msa_set"))
                validate_msa_set(files + ([args.base_msa] if args.base_msa else []))
            except MSASetError as e:
                _fail(str(e))
            if len(files) < 3:
                _fail(f"at least 3 MSAs are needed for the trees (found {len(files)}).")
            base = args.base_msa
            if not base:
                ranking = rank_msas(files, at._is_nucleotide(files[0]), os.path.join(tmp, "ranking"))
                ranking.insert(0, "rank", range(1, len(ranking) + 1))
                ranking.to_csv(f"{prefix}_DL_ranking.csv", index=False)
                base = os.path.join(tmp, "msa_set", ranking.loc[0, "code"])
                shutil.copyfile(base, f"{prefix}_base_MSA.fasta")
                print(f"[guidance3-tree] base MSA chosen by the DL model: {ranking.loc[0, 'code']} "
                      f"(ranking: {prefix}_DL_ranking.csv)")
            names = [n for n, _s in read_fasta(files[0])]
        else:
            try:
                alt_trees = list(Phylo.parse(args.alt_trees, "newick"))
            except Exception as e:
                _fail(f"could not read --alt-trees as Newick: {e}")
            if len(alt_trees) < 3:
                _fail(f"--alt-trees must contain at least 3 trees (found {len(alt_trees)}).")
            names = sorted(t.name for t in alt_trees[0].get_terminals())
            for i, t in enumerate(alt_trees[1:], 2):
                if sorted(x.name for x in t.get_terminals()) != names:
                    _fail(f"tree {i} of --alt-trees does not have the same leaves as tree 1.")
            base = args.base_msa
            if base:
                base_names = sorted(n for n, _s in read_fasta(base))
                if base_names != names:
                    diff = sorted(set(base_names) ^ set(names))[:5]
                    _fail(f"--base-msa and --alt-trees do not have the same sequence names (e.g. {', '.join(diff)}).")

        # Work on short codes, as the pipeline does: IQ-TREE rewrites names with special characters.
        name_to_code = {n: f"S{i}" for i, n in enumerate(names, 1)}
        codes = {c: n for n, c in name_to_code.items()}
        coded_dir = os.path.join(tmp, "coded")
        os.makedirs(coded_dir)
        if alt_trees is None:
            for f in files:
                dst = os.path.join(coded_dir, os.path.basename(f))
                _write_coded(f, dst, name_to_code)
                alt_files.append(dst)
        else:
            for t in alt_trees:
                for leaf in t.get_terminals():
                    leaf.name = name_to_code[leaf.name]
        coded_base = None
        if base:
            coded_base = os.path.join(tmp, "base_msa.coded.fasta")
            _write_coded(base, coded_base, name_to_code)

        ok = at.build_alignment_trees(alt_files, coded_base, codes, prefix, os.path.join(tmp, "iqtree"), IQTREE,
                                      args.proc_num, log_file, n_trees=args.n_trees,
                                      want_consensus=not args.no_consensus, user_tree_file=args.user_tree or "",
                                      alt_trees=alt_trees)
    finally:
        if args.keep_temp:
            print(f"[guidance3-tree] intermediate files kept in {tmp}")
        else:
            shutil.rmtree(tmp, ignore_errors=True)

    outputs = sorted(f for f in os.listdir(args.out_dir) if f.startswith(args.prefix + "_"))
    print("[guidance3-tree] outputs in " + args.out_dir + ":\n  " + "\n  ".join(outputs))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
