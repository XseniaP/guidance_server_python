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
  5. Copies the top-scoring MSA to <out-dir>/best_msa.fasta and prints its name.
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
    args = _build_parser().parse_args(argv)

    msas_dir = os.path.abspath(args.msas_dir)
    out_dir  = os.path.abspath(args.out_dir)
    os.makedirs(out_dir, exist_ok=True)

    fasta_files = sorted(glob.glob(os.path.join(msas_dir, "*.fasta")))
    if not fasta_files:
        sys.exit(f"No .fasta files found in: {msas_dir}")
    print(f"Found {len(fasta_files)} alternative MSAs in {msas_dir}")

    is_nuc = (args.seq_type == "nuc")
    models_list  = NUC_MODELS  if is_nuc else AA_MODELS
    k_values     = [8, 16, 32] if is_nuc else [5, 10, 20]
    stats_output = ["ALL"]     if is_nuc else ["ALL_NO_SUBS_MATRIX"]
    model_path   = DL_MODEL_NUC_PATH    if is_nuc else DL_MODEL_PATH
    scaler_path  = DL_MODEL_NUC_SCALER_PATH if is_nuc else DL_MODEL_SCALER_PATH

    work_dir = tempfile.mkdtemp(prefix="guidance3_predict_")
    try:
        # features_for_msas expects: input_files_dir_path/<subfolder>/*.fasta
        # The subfolder name is arbitrary; we mirror what select_best_msa uses.
        features_input_dir = os.path.join(work_dir, "features_input")
        all_msas_dir = os.path.join(features_input_dir, "all_msas")
        os.makedirs(all_msas_dir)

        for f in fasta_files:
            _copy_stripped(f, all_msas_dir)

        # Pick any MSA as the reference placeholder required by features_for_msas.
        # The TRUE label is not used during inference — only needed structurally.
        ref_msa = random.choice(fasta_files)
        _copy_stripped(ref_msa, os.path.join(all_msas_dir, TRUE_MSA_FILENAME))
        if args.verbose:
            print(f"Reference placeholder: {os.path.basename(ref_msa)} → {TRUE_MSA_FILENAME}")

        features_cfg = {
            "models_list": models_list,
            "sop_calc_type": 1,
            "additional_weights": [
                "HENIKOFF_WG", "HENIKOFF_WOG",
                "CLUSTAL_MID_ROOT", "CLUSTAL_DIFFERENTIAL_SUM",
            ],
            "k_values": k_values,
            "stats_output": stats_output,
            "input_files_dir_path": features_input_dir,  # parent of all_msas/
            "output_file_dir_path": work_dir,
            "is_unified_file": True,
            "matrix_dir_path": FEATURES_EXTRACTION_MATRIX_DIR,
        }
        config_path = os.path.join(work_dir, "config_features.json")
        with open(config_path, "w") as f:
            json.dump(features_cfg, f, indent=4)

        # Step 1: feature extraction
        print("Extracting MSA features...")
        features_cmd = [FEATURES_EXTRACTION_PROG, config_path]
        if args.verbose:
            print(f"  {' '.join(features_cmd)}")
        result = subprocess.run(features_cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            sys.exit(f"Feature extraction failed (rc={result.returncode}):\n{result.stderr}")
        if args.verbose and result.stdout:
            print(result.stdout)

        features_file = os.path.join(work_dir, "unified_stats_all_msas.csv")
        if not os.path.exists(features_file):
            sys.exit(f"Features file not produced: {features_file}")

        df = pd.read_csv(features_file)
        run_id = os.path.basename(msas_dir)
        df.insert(1, "code1", run_id)
        df.dropna(axis=1, how="all", inplace=True)
        df.to_csv(features_file, index=False)

        # Save features to out_dir: drop label columns and non-candidate rows.
        label_cols = [c for c in df.columns
                      if c.startswith("dseq") or c.startswith("dpos") or c == "ssp_from_true"]
        features_out = df.drop(columns=label_cols, errors="ignore")
        features_out = features_out[
            ~features_out["code"].isin([TRUE_MSA_FILENAME, "all_msas"])
        ]
        features_out.to_csv(os.path.join(out_dir, "features.csv"), index=False)

        # Step 2: DL model prediction
        print("Running DL model...")
        pred_out_dir = os.path.join(work_dir, "pred")
        os.makedirs(pred_out_dir, exist_ok=True)
        predict_cmd = DL_MODEL_PYTHON + [
            DL_MODEL_PREDICT_SCRIPT,
            "--features-file", features_file,
            "--model-path",    model_path,
            "--scaler-path",   scaler_path,
            "--scaler-type-features", "rank",
            "--scaler-type-labels",   "rank",
            "--out-dir",  pred_out_dir,
            "--no-metrics",
        ]
        if args.verbose:
            print(f"  {' '.join(predict_cmd)}")
        result = subprocess.run(predict_cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            sys.exit(f"DL model failed (rc={result.returncode}):\n{result.stderr}")
        if args.verbose:
            print(result.stdout)

        pred_file = os.path.join(pred_out_dir, "prediction_pretrained_0_mode1_dseq_from_true.csv")
        if not os.path.exists(pred_file):
            sys.exit(f"Prediction file not produced: {pred_file}")

        preds = pd.read_csv(pred_file)
        if "code" not in preds.columns or "predicted_score" not in preds.columns:
            sys.exit("Prediction CSV is missing expected columns (code, predicted_score)")

        candidates = preds[preds["code"] != TRUE_MSA_FILENAME].copy()
        if candidates.empty:
            sys.exit("No candidate predictions found (all rows were the reference placeholder)")

        best_row   = candidates.loc[candidates["predicted_score"].idxmax()]
        best_code  = best_row["code"]
        best_score = best_row["predicted_score"]

        best_src = os.path.join(all_msas_dir, best_code)
        if not os.path.exists(best_src):
            sys.exit(f"Best MSA file not found in input dir: {best_code}")

        best_dst = os.path.join(out_dir, "best_msa.fasta")
        shutil.copy(best_src, best_dst)

        # Save predictions (all candidates, sorted best-first) to out_dir.
        preds_out = preds[preds["code"] != TRUE_MSA_FILENAME].copy()
        preds_out = preds_out.sort_values("predicted_score", ascending=False)
        preds_out.to_csv(os.path.join(out_dir, "predictions.csv"), index=False)

        print(f"Best MSA : {best_code}  (predicted score: {best_score:.6f})")
        print(f"Saved to : {best_dst}")
        print(f"Features : {os.path.join(out_dir, 'features.csv')}")
        print(f"All predictions : {os.path.join(out_dir, 'predictions.csv')}")

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
