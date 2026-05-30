"""
subsample_predict_main.py

For each sample size N in SAMPLE_SIZES, samples N rows per code1 from the full
test set, runs predict_pretrained_main, then runs pick_best analysis.
Outputs are organised under --out-dir/n{N}/.

A cross-level summary CSV is written to --out-dir/subsampling_summary.csv.

Example
-------
python scripts/subsample_predict_main.py \
    --test-file out/clean_data/orthomam_test_set.csv \
    --model-path  out/OrthoMaM12/new_model1_1_340/regressor_model_0_mode1_dseq_from_true.keras \
    --scaler-path out/OrthoMaM12/new_model1_1_340/scaler_0_mode1_dseq_from_true.pkl \
    --out-dir     out/subsampled_predictions
"""
from __future__ import annotations

import argparse
import sys

import pandas as pd

from guidance3.dl_model.pipeline.pretrained_predictor import PretrainedPredictConfig, PretrainedPredictor
from guidance3.dl_model.evaluation.pick_best import PickBest

DEFAULT_SAMPLE_SIZES = [20, 50, 100, 200, 400, 800, 1600]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _is_default(codes: pd.Series, code1: str) -> pd.Series:
    """Boolean mask for the four default-aligner rows within a code1 group."""
    mafft   = codes == "MSA.MAFFT.aln.With_Names"
    prank   = codes.isin(["MSA.PRANK.aln.With_Names", "MSA.PRANK.aln.best.fas"])
    muscle  = codes.str.match(r"^MSA\.MUSCLE\.aln\.best\.\S+\.fas$", na=False)
    baliphy = codes.str.match(r"^MSA\.BALIPHY\.aln\.best\.\S+\.fas$", na=False)
    return mafft | prank | muscle | baliphy


def _sample_group(g: pd.DataFrame, n: int, seed: int) -> pd.DataFrame:
    """
    Always keep the four default-aligner rows; fill the remaining quota
    (up to n total) with a random sample from the rest.
    """
    code1 = g["code1"].iloc[0]
    default_mask = _is_default(g["code"], code1)
    defaults = g[default_mask]
    rest     = g[~default_mask]

    remaining_quota = max(0, n - len(defaults))
    if remaining_quota > 0 and len(rest) > 0:
        sampled_rest = rest.sample(min(remaining_quota, len(rest)), random_state=seed)
    else:
        sampled_rest = rest.iloc[0:0]  # empty with same columns

    return pd.concat([defaults, sampled_rest])


def sample_per_code1(df: pd.DataFrame, n: int, seed: int) -> pd.DataFrame:
    """
    Return up to *n* rows per code1, always keeping the four default-aligner
    rows (MAFFT, PRANK, MUSCLE, BALIPHY) and filling the remainder randomly.
    """
    return (
        df.groupby("code1", group_keys=False)
        .apply(lambda g: _sample_group(g, n, seed))
        .reset_index(drop=True)
    )


def run_predict(subset_csv: str, n: int, args: argparse.Namespace, out_dir: Path) -> str:
    """Run PretrainedPredictor on *subset_csv*, return path to predictions CSV."""
    cfg = PretrainedPredictConfig(
        features_file=subset_csv,
        true_score_name=args.true_score_name,
        mode=args.mode,
        remove_correlated_features=args.remove_correlated_features,
        corr_threshold=args.corr_threshold,
        scaler_type_features=args.scaler_type_features,
        scaler_type_labels=args.scaler_type_labels,
        model_path=args.model_path,
        scaler_path=args.scaler_path,
        out_dir=str(out_dir),
        run_id=f"n{n}",
        verbose=(not args.quiet),
    )
    result = PretrainedPredictor(cfg).run(custom_objects=None)
    return result["predictions_csv"]


def run_pickbest(subset_csv: str, pred_csv: str, args: argparse.Namespace, out_dir: Path) -> PickBest:
    """Run PickBest on *subset_csv* features + *pred_csv* predictions."""
    pickme = PickBest(
        features_file=subset_csv,
        prediction_file=pred_csv,
        true_score_name=args.true_score_name,
        error=args.pickbest_error,
        subset=None,
        output_dir=str(out_dir),
    )
    pickme.run(0)
    pickme.summarize()
    pickme.save_to_csv(0)
    pickme.plot_results(0)
    pickme.plot_overall_results(0)
    return pickme


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Predict with a pretrained model at multiple subsampling levels, then run pick_best."
    )
    p.add_argument("--test-file", required=True,
                   help="Full test-set CSV (e.g. out/clean_data/orthomam_test_set.csv).")
    p.add_argument("--model-path",  required=True, help="Path to pretrained .keras model.")
    p.add_argument("--scaler-path", required=True, help="Path to pretrained scaler (.pkl or dir).")

    p.add_argument("--sample-sizes", type=int, nargs="+", default=DEFAULT_SAMPLE_SIZES,
                   metavar="N",
                   help=f"Number of codes to sample per code1. Default: {DEFAULT_SAMPLE_SIZES}")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for reproducible sampling (default: 42).")

    p.add_argument("--out-dir", default="out/subsampled_predictions",
                   help="Base output directory. A subdirectory n{N} is created for each sample size.")
    p.add_argument("--true-score-name", default="dseq_from_true",
                   choices=["ssp_from_true", "dseq_from_true", "dpos_from_true"])
    p.add_argument("--mode", type=int, default=1, choices=[1, 3])
    p.add_argument("--remove-correlated-features", action="store_true")
    p.add_argument("--corr-threshold", type=float, default=0.90)
    p.add_argument("--scaler-type-features", default="standard",
                   choices=["standard", "rank", "zscore"])
    p.add_argument("--scaler-type-labels", default="standard",
                   choices=["standard", "rank", "zscore"])

    p.add_argument("--no-pickbest", action="store_true",
                   help="Skip PickBest analysis; only run predictions.")
    p.add_argument("--pickbest-error", type=float, default=0.0,
                   help="Error tolerance passed to PickBest (default: 0.0).")

    p.add_argument("--quiet", action="store_true", help="Suppress verbose output.")
    return p


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    args = build_parser().parse_args(argv)

    test_file = Path(args.test_file)
    if not test_file.exists():
        print(f"ERROR: test file not found: {test_file}", file=sys.stderr)
        return 1

    base_out = Path(args.out_dir)
    base_out.mkdir(parents=True, exist_ok=True)

    if not args.quiet:
        print(f"Loading full test set: {test_file}")
    full_df = pd.read_csv(test_file)
    full_df["code1"] = full_df["code1"].astype(str)
    if not args.quiet:
        print(f"  {len(full_df):,} rows | {full_df['code1'].nunique()} unique code1 values")

    summary_rows = []

    for n in sorted(args.sample_sizes):
        print(f"\n{'='*60}")
        print(f"  Sample size: {n} codes per code1")
        print(f"{'='*60}")

        n_out = base_out / f"n{n}"
        n_out.mkdir(parents=True, exist_ok=True)

        # --- subsample --------------------------------------------------
        subset = sample_per_code1(full_df, n, args.seed)
        subset_csv = str(n_out / "subset_features.csv")
        subset.to_csv(subset_csv, index=False)
        if not args.quiet:
            actual_per_code1 = subset.groupby("code1").size()
            print(f"  Subset: {len(subset):,} rows "
                  f"(min {actual_per_code1.min()}, max {actual_per_code1.max()} per code1)")
            print(f"  Saved:  {subset_csv}")

        # --- predict ----------------------------------------------------
        pred_csv = run_predict(subset_csv, n, args, n_out)
        if not args.quiet:
            print(f"  Predictions: {pred_csv}")

        # --- aligner counts in the subset --------------------------------
        aligner_totals = (
            subset["aligner"].value_counts().to_dict()
            if "aligner" in subset.columns
            else {}
        )
        aligner_mean_per_code1 = (
            subset.groupby("code1")["aligner"].value_counts()
            .unstack(fill_value=0)
            .mean()
            .to_dict()
            if "aligner" in subset.columns
            else {}
        )

        row: dict = {
            "n_per_code1": n,
            "n_rows_total": len(subset),
            "n_code1": full_df["code1"].nunique(),
        }
        for aligner, total in sorted(aligner_totals.items()):
            row[f"n_{aligner}_total"] = int(total)
            row[f"n_{aligner}_mean_per_code1"] = round(aligner_mean_per_code1.get(aligner, 0), 2)
        row["predictions_csv"] = pred_csv

        # --- pick best --------------------------------------------------
        if not args.no_pickbest:
            try:
                pickme = run_pickbest(subset_csv, pred_csv, args, n_out)

                overall_counts = (
                    pickme.overall_winners_df["overall_winner"]
                    .value_counts()
                    .to_dict()
                )
                row.update({f"overall_{k}": v for k, v in overall_counts.items()})

                aligner_counts = {
                    col: pickme.winners_df[col].value_counts().to_dict()
                    for col in pickme.winners_df.columns
                }
                for aligner, counts in aligner_counts.items():
                    for winner, cnt in counts.items():
                        row[f"{aligner}_{winner}"] = cnt

                print(f"  PickBest overall winners: {overall_counts}")
            except Exception as exc:
                print(f"  WARNING: PickBest failed for n={n}: {exc}")
                row["pickbest_error"] = str(exc)

        summary_rows.append(row)

    # --- cross-level summary --------------------------------------------
    summary_df = pd.DataFrame(summary_rows).set_index("n_per_code1")
    summary_path = base_out / "subsampling_summary.csv"
    summary_df.to_csv(summary_path)
    print(f"\nDone. Summary: {summary_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
