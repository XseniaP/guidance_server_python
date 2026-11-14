#pragma once
#include "msa_basic_stats.h"
#include "sp_score.h"
#include "rooted_tree.h"
#include "unrooted_tree.h"
#include "enums.h"
#include <vector>
#include <string>
#include <set>
#include <unordered_map>
#include <memory>

class WSopStats : public BasicStats {
public:
    double sp_HENIKOFF_with_gaps;
    double sp_HENIKOFF_without_gaps;
    double sp_CLUSTAL_WEIGHTS_mid_root;
    double sp_CLUSTAL_WEIGHTS_diff_sum;

    std::vector<std::string> weight_names;
    std::vector<std::vector<double>> seq_weights_options;
    std::unordered_map<std::string, std::unique_ptr<RootedTree>> rooted_trees;

    WSopStats(const std::string& code, int taxa_num, int msa_length);

    std::pair<std::vector<double>, std::vector<double>> compute_seq_w_henikoff_vars(
        const std::vector<std::string>& sequences) const;

    std::vector<double> get_weight_list(const UnrootedTree& tree, RootingMethods rooting_method,
                                         const std::vector<std::string>& seq_names);

    void calc_seq_weights(const std::set<WeightMethods>& additional_weights,
                           const std::vector<std::string>& sequences,
                           const std::vector<std::string>& seq_names,
                           const UnrootedTree& tree);

    void calc_w_sp(const std::vector<std::string>& sequences, const SPScore& sp);

    std::vector<std::string> get_ordered_col_names_with_model(const std::string& model_name,
                                                               double go_val, double ge_val) const;

    std::vector<StatValue> get_my_features_as_list() const override;
};

void normalize(std::vector<double>& v);
