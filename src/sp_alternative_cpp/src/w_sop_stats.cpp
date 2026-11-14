#include "w_sop_stats.h"
#include <numeric>
#include <iostream>
#include <chrono>
#include "utils.h"

WSopStats::WSopStats(const std::string& code, int taxa_num, int msa_length)
    : BasicStats(code, taxa_num, msa_length,
                  {"code", "sp_HENIKOFF_with_gaps", "sp_HENIKOFF_without_gaps",
                   "sp_CLUSTAL_WEIGHTS_mid_root", "sp_CLUSTAL_WEIGHTS_diff_sum"}),
      sp_HENIKOFF_with_gaps(-1), sp_HENIKOFF_without_gaps(-1),
      sp_CLUSTAL_WEIGHTS_mid_root(-1), sp_CLUSTAL_WEIGHTS_diff_sum(-1) {}

std::pair<std::vector<double>, std::vector<double>> WSopStats::compute_seq_w_henikoff_vars(
    const std::vector<std::string>& sequences) const {
    std::vector<double> seq_weights_with_gap(taxa_num, 0.0);
    std::vector<double> seq_weights_no_gap(taxa_num, 0.0);

    for (int k = 0; k < msa_length; k++) {
        std::unordered_map<char, std::vector<int>> seq_dict;
        for (int i = 0; i < taxa_num; i++) {
            char ch = sequences[i][k];
            seq_dict[ch].push_back(i);
        }
        for (const auto& [cluster_key, indices] : seq_dict) {
            double w = 1.0 / indices.size();
            for (int seq_inx : indices) {
                seq_weights_with_gap[seq_inx] += w;
                if (cluster_key != '-') {
                    seq_weights_no_gap[seq_inx] += w;
                }
            }
        }
    }

    double sum_wg = std::accumulate(seq_weights_with_gap.begin(), seq_weights_with_gap.end(), 0.0);
    double sum_ng = std::accumulate(seq_weights_no_gap.begin(), seq_weights_no_gap.end(), 0.0);
    for (int i = 0; i < taxa_num; i++) {
        seq_weights_with_gap[i] /= sum_wg;
        seq_weights_no_gap[i] /= sum_ng;
    }
	normalize(seq_weights_with_gap);
	normalize(seq_weights_no_gap);
    return {seq_weights_with_gap, seq_weights_no_gap};
}

std::vector<double> WSopStats::get_weight_list(const UnrootedTree& tree, RootingMethods rooting_method,
                                                 const std::vector<std::string>& seq_names) {
    std::string key = rooting_method_to_string(rooting_method);
    RootingPoint rp = get_rooting_point(rooting_method, tree);
    auto rt = std::make_unique<RootedTree>(RootedTree(tree, rp));
    string x = rt->print_newick();
    rt->calc_seq_w();
    std::vector<double> weights;
    for (const auto& s_name : seq_names) {
        weights.push_back(rt->seq_weight_dict[s_name]);
    }

    rooted_trees[key] = std::move(rt);
 //   for (size_t i = 0; i < weights.size(); i++) {
 //       std::cout << "Weight for sequence " << seq_names[i] << " with method " << key << ": " << weights[i] << std::endl;
	//}

    normalize(weights);
    return weights;
}

void normalize(std::vector<double>& v) {
    double sum = std::accumulate(v.begin(), v.end(), 0.0);

    if (sum == 0.0) {
        // Handle edge case: avoid division by zero
        return; // or throw, or set uniform values
    }

    for (double& x : v) {
        x /= sum;
    }
}

void WSopStats::calc_seq_weights(const std::set<WeightMethods>& additional_weights,
                                   const std::vector<std::string>& sequences,
                                   const std::vector<std::string>& seq_names,
                                   const UnrootedTree& tree) {
    if (additional_weights.empty()) return;

    if (additional_weights.count(WeightMethods::HENIKOFF_WG) ||
        additional_weights.count(WeightMethods::HENIKOFF_WOG)) {
        auto [wg, wog] = compute_seq_w_henikoff_vars(sequences);
        if (additional_weights.count(WeightMethods::HENIKOFF_WG)) {
            seq_weights_options.push_back(wg);
            weight_names.push_back(weight_method_to_string(WeightMethods::HENIKOFF_WG));
        }
        if (additional_weights.count(WeightMethods::HENIKOFF_WOG)) {
            seq_weights_options.push_back(wog);
            weight_names.push_back(weight_method_to_string(WeightMethods::HENIKOFF_WOG));
        }
    }

    if (additional_weights.count(WeightMethods::CLUSTAL_MID_ROOT)) {
        seq_weights_options.push_back(get_weight_list(tree, RootingMethods::LONGEST_PATH_MID, seq_names));
        weight_names.push_back(weight_method_to_string(WeightMethods::CLUSTAL_MID_ROOT));
        // cout << "End weights calc for CLUSTAL_MID_ROOT " << endl;
    }

    if (additional_weights.count(WeightMethods::CLUSTAL_DIFFERENTIAL_SUM)) {
        seq_weights_options.push_back(get_weight_list(tree, RootingMethods::MIN_DIFFERENTIAL_SUM, seq_names));
        weight_names.push_back(weight_method_to_string(WeightMethods::CLUSTAL_DIFFERENTIAL_SUM));
        // cout << "End weights calc for CLUSTAL_DIFFERENTIAL_SUM " << endl;
    }
}

void WSopStats::calc_w_sp(const std::vector<std::string>& sequences, const SPScore& sp) {
    std::unordered_map<std::string, double> sop_w_options_dict;
    for (size_t index = 0; index < weight_names.size(); index++) {
        sop_w_options_dict[weight_names[index]] = sp.compute_efficient_w_sp(sequences, seq_weights_options[index]);
    }

    auto get_or_zero = [&](WeightMethods wm) -> double {
        std::string key = weight_method_to_string(wm);
        auto it = sop_w_options_dict.find(key);
        return (it != sop_w_options_dict.end()) ? it->second : 0.0;
    };

    sp_HENIKOFF_with_gaps = get_or_zero(WeightMethods::HENIKOFF_WG);
    sp_HENIKOFF_without_gaps = get_or_zero(WeightMethods::HENIKOFF_WOG);
    sp_CLUSTAL_WEIGHTS_mid_root = get_or_zero(WeightMethods::CLUSTAL_MID_ROOT);
    sp_CLUSTAL_WEIGHTS_diff_sum = get_or_zero(WeightMethods::CLUSTAL_DIFFERENTIAL_SUM);
}

std::vector<std::string> WSopStats::get_ordered_col_names_with_model(const std::string& model_name,
                                                                       double go_val, double ge_val) const {
    std::string suffix = get_model_name_suffix(model_name, go_val, ge_val);
    std::vector<std::string> result;
    for (const auto& col_name : ordered_col_names) {
        result.push_back(col_name + suffix);
    }
    return result;
}

std::vector<StatValue> WSopStats::get_my_features_as_list() const {
    return {
        StatValue(code),
        StatValue(sp_HENIKOFF_with_gaps),
        StatValue(sp_HENIKOFF_without_gaps),
        StatValue(sp_CLUSTAL_WEIGHTS_mid_root),
        StatValue(sp_CLUSTAL_WEIGHTS_diff_sum)
    };
}
