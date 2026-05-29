import logging

from dl_model.config.config import DataConfig, FeatureConfig, TrainConfig, OutputConfig, ShapExplainConfig
from dl_model.pipeline.experiment import RegressionExperiment

logging.basicConfig(level=logging.INFO)

data_cfg = DataConfig(
    features_file="/Users/kpolonsky/Documents/sp_alternative/dl_model/out/clean_data/nucleotide_features_200426.csv", #replace with your features file
    true_score_name="dseq_from_true",
    test_size=0.2,
    deduplicated=False,
    empirical=False, #change for empirical
    min_rows_per_code_after_dedup=1,
    number_of_msas_threshold_simulated=400,
    number_of_msas_threshold_empirical=400,
    random_state=12345
)

feat_cfg = FeatureConfig(
    mode=1,
    remove_correlated_features=False,
    scaler_type_features="rank",
    scaler_type_labels="rank",
)

train_cfg = TrainConfig(
    epochs=50,
    batch_size=32,
    learning_rate=0.0022,
    neurons=(64, 128, 64, 512),
    dropout_rate=0.24,
    regularizer_name="l2",
    l1=0,
    l2=1.65e-5,

    loss_fn="custom_mse",
    alpha=0,
    eps=0,
    top_k=8,
    ranking_weight=1.33,
    margin=0.0,

    batch_generation="custom",   # "custom" (for model2) or "standard" (for model1)
    repeats=1,
    mixed_portion=0.0,
    per_aligner=False,
)

out_cfg = OutputConfig(
    out_dir="../out",
    run_id="0",
    verbose=True,
    save_model=True,
    save_scaled_csv=True,
    save_predictions_csv=True,
    save_plots=True,
    save_scaler=True,
)

explain_cfg = ShapExplainConfig(
    enabled=True,
    sample_n=500
)

exp = RegressionExperiment(data_cfg, feat_cfg, train_cfg, out_cfg, explain_cfg)
out = exp.run()
print(out["metrics"])
