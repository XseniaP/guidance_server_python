"""
guidance3-rank: rank a set of MSAs of the same sequences with the GUIDANCE3 DL model and return the best one.

Usage:
    guidance3-rank --msas-dir <folder | .zip | .tar.gz> --out-dir <dir> [--seq-type {aa,nuc,codon}]

The model predicts each MSA's distance from the (unknown) true MSA; the MSA with the LOWEST predicted score is
ranked first. Writes to --out-dir:
  ranking.csv      rank, MSA file, predicted_score (best first)
  best_msa.fasta   the top-ranked MSA (unchanged)
  features.csv     the alignment features the model used
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import tempfile

TRUE_MSA_FILENAME = "MSA_TRUE.fas"

AA_MODELS = [
    {"gap_open_cost": -10, "gap_extend_cost": -0.5, "matrix_file_name": "BLOSUM62"},
    {"gap_open_cost": -6,  "gap_extend_cost": -0.5, "matrix_file_name": "BLOSUM62"},
    {"gap_open_cost": -10, "gap_extend_cost": -1,   "matrix_file_name": "BLOSUM62"},
    {"gap_open_cost": -6,  "gap_extend_cost": -1,   "matrix_file_name": "BLOSUM62"},
    {"gap_open_cost": -10, "gap_extend_cost": -0.2, "matrix_file_name": "BLOSUM62"},
    {"gap_open_cost": -6,  "gap_extend_cost": -0.2, "matrix_file_name": "BLOSUM62"},
    {"gap_open_cost": -10, "gap_extend_cost": -0.5, "matrix_file_name": "PAM250"},
    {"gap_open_cost": -6,  "gap_extend_cost": -0.5, "matrix_file_name": "PAM250"},
    {"gap_open_cost": -10, "gap_extend_cost": -1,   "matrix_file_name": "PAM250"},
    {"gap_open_cost": -6,  "gap_extend_cost": -1,   "matrix_file_name": "PAM250"},
    {"gap_open_cost": -10, "gap_extend_cost": -0.2, "matrix_file_name": "PAM250"},
    {"gap_open_cost": -6,  "gap_extend_cost": -0.2, "matrix_file_name": "PAM250"},
]
NUC_MODELS = [
    {"gap_open_cost": -1,   "gap_extend_cost": 0, "matrix_file_name": "NucleotidesPAM250"},
    {"gap_open_cost": -1.5, "gap_extend_cost": 0, "matrix_file_name": "NucleotidesPAM250"},
    {"gap_open_cost": -3,   "gap_extend_cost": 0, "matrix_file_name": "NucleotidesPAM250"},
    {"gap_open_cost": -6,   "gap_extend_cost": 0, "matrix_file_name": "NucleotidesPAM250"},
]


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="guidance3-rank",
        description="Rank a set of MSAs of the same sequences with the GUIDANCE3 DL model (predicted distance "
                    "from the true MSA; lower = better) and return the best one.")
    p.add_argument("--msas-dir", required=True,
                   help="Folder, or .zip / .tar.gz archive, of aligned FASTA files of the same sequences "
                        "(extensions: .fasta .fas .fa .faa .fna .aln .afa .mfa).")
    p.add_argument("--out-dir", required=True,
                   help="Output directory: ranking.csv, best_msa.fasta, features.csv.")
    p.add_argument("--seq-type", choices=["aa", "nuc", "codon"], default=None,
                   help="Sequence type. Default: detected (aa or nuc). Use codon for codon alignments, which are "
                        "ranked on their amino-acid translation.")
    p.add_argument("--verbose", action="store_true", default=False,
                   help="Print detailed progress to stdout.")
    return p


def main(argv=None) -> int:
    from guidance3.sequences.msa_set import MSASetError, select_best_from_msa_set
    from guidance3.trees.alignment_trees import _is_nucleotide

    args = _build_parser().parse_args(argv)
    src = os.path.abspath(args.msas_dir)
    out_dir = os.path.abspath(args.out_dir)
    if not os.path.exists(src):
        sys.exit(f"guidance3-rank: error: not found: {src}")
    os.makedirs(out_dir, exist_ok=True)

    work = tempfile.mkdtemp(prefix="guidance3_rank_", dir=out_dir)
    try:
        seq_type = args.seq_type
        try:
            if seq_type is None:
                from guidance3.sequences.msa_set import collect_msa_files
                probe = collect_msa_files(src, os.path.join(work, "probe"))
                if not probe:
                    raise MSASetError("no alignment files found (extensions: .fasta .fas .fa .faa .fna .aln "
                                      ".afa .mfa).")
                seq_type = "nuc" if _is_nucleotide(probe[0]) else "aa"
                shutil.rmtree(os.path.join(work, "probe"), ignore_errors=True)
                print(f"Sequence type: {seq_type} (detected; use --seq-type codon for codon alignments)")
            best = select_best_from_msa_set(src, work, seq_type, verbose=args.verbose)
        except MSASetError as e:
            sys.exit(f"guidance3-rank: error: {e}")
        except RuntimeError as e:
            sys.exit(f"guidance3-rank: error: {e}")

        ranking_csv = os.path.join(out_dir, "ranking.csv")
        shutil.copyfile(os.path.join(work, "MSA_set_ranking.csv"), ranking_csv)
        shutil.copyfile(best, os.path.join(out_dir, "best_msa.fasta"))
        features = os.path.join(work, "MSA_set_ranking", "features.csv")
        if os.path.exists(features):
            shutil.copyfile(features, os.path.join(out_dir, "features.csv"))
    finally:
        shutil.rmtree(work, ignore_errors=True)

    print(f"Ranking (best first): {ranking_csv}")
    print(f"Best MSA            : {os.path.join(out_dir, 'best_msa.fasta')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
