"""
guidance3-predict: select the best MSA from a folder of alternatives.

Usage:
    guidance3-predict --msas-dir <dir> --seq-type {aa,nuc} --out-dir <dir>

Steps:
  1. Copies all .fasta files from --msas-dir into a temporary workspace.
  2. Picks one MSA at random as the _TRUE.fas reference required by
     features_for_msas (only a structural placeholder — not used for scoring).
  3. Runs the bundled features_for_msas binary to extract alignment features.
  4. Runs the bundled DL model to predict MSA quality scores.
  5. Copies the MSA with the LOWEST predicted score (predicted distance from the true MSA) to
     <out-dir>/best_msa.fasta and prints its name.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

from guidance3.constants import (
    FEATURES_EXTRACTION_PROG,
    FEATURES_EXTRACTION_MATRIX_DIR,
    DL_MODEL_PYTHON,
    DL_MODEL_PREDICT_SCRIPT,
    DL_MODEL_PATH,
    DL_MODEL_SCALER_PATH,
    DL_MODEL_NUC_PATH,
    DL_MODEL_NUC_SCALER_PATH,
)

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
        description="Select the best MSA from a folder of alternatives using the GUIDANCE3 DL model."
    )
    p.add_argument("--msas-dir", required=True,
                   help="Directory containing alternative MSA .fasta files.")
    p.add_argument("--seq-type", required=True, choices=["aa", "nuc"],
                   help="Sequence type: aa (amino acid) or nuc (nucleotide).")
    p.add_argument("--out-dir", required=True,
                   help="Output directory. Best MSA is written as best_msa.fasta.")
    p.add_argument("--verbose", action="store_true", default=False,
                   help="Print detailed progress to stdout.")
    return p


def _copy_stripped(src: str, dst: str) -> None:
    """Copy a FASTA stripping blank lines — features_for_msas crashes on empty lines."""
    if os.path.isdir(dst):
        dst = os.path.join(dst, os.path.basename(src))
    with open(src) as fin, open(dst, "w") as fout:
        for line in fin:
            if line.strip():
                fout.write(line)


def main(argv=None) -> int:
    from guidance3.sequences.msa_set import rank_msas

    args = _build_parser().parse_args(argv)

    msas_dir = os.path.abspath(args.msas_dir)
    out_dir  = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    fasta_files = sorted(glob.glob(os.path.join(msas_dir, "*.fasta")))
    if not fasta_files:
        sys.exit(f"No .fasta files found in: {msas_dir}")
    print(f"Found {len(fasta_files)} alternative MSAs in {msas_dir}")

    # The model predicts each MSA's distance from the true MSA: the LOWEST predicted_score is the best MSA
    # (as in guidance3.dl_model.evaluation.pick_best and the pipeline's select_best_msa).
    try:
        preds = rank_msas(fasta_files, args.seq_type == "nuc", out_dir, verbose=args.verbose)
    except RuntimeError as e:
        sys.exit(str(e))

    best_code, best_score = preds.loc[0, "code"], preds.loc[0, "predicted_score"]
    best_dst = os.path.join(out_dir, "best_msa.fasta")
    shutil.copy(os.path.join(msas_dir, best_code), best_dst)
    preds.to_csv(os.path.join(out_dir, "predictions.csv"), index=False)  # best first

    print(f"Best MSA : {best_code}  (predicted distance from the true MSA: {best_score:.6f}; lower = better)")
    print(f"Saved to : {best_dst}")
    print(f"Features : {os.path.join(out_dir, 'features.csv')}")
    print(f"All predictions (best first) : {os.path.join(out_dir, 'predictions.csv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
