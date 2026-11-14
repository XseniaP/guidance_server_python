#pragma once
#include "msa_basic_stats.h"
#include "sp_score.h"
#include "enums.h"
#include <vector>
#include <string>
#include <set>

class SopStats : public BasicStats {
public:
    int sp_match_count;
    double sp_match_count_norm;
    int sp_mismatch_count;
    double sp_mismatch_count_norm;
    int sp_go;
    double sp_go_norm;
    int sp_ge;
    double sp_ge_norm;

    double sp_match;
    double sp_match_norm;
    double sp_mismatch;
    double sp_mismatch_norm;
    double sp;
    double sp_norm;

    std::set<std::string> model_agnostic_col_names;
    std::set<std::string> gaps_agnostic_col_names;

    SopStats(const std::string& code, int taxa_num, int msa_length);

    void set_my_sop_score_parts(const SPScore& sp_score, const std::vector<std::string>& sequences,
                                vector<vector<int>>& subs_matrix_counts, const set<StatsOutput>& statsOutput,
                                const int iteration);
    std::vector<std::string> get_ordered_col_names_with_model(const std::string& model_name,
                                                              double go_val, double ge_val) const;
    std::vector<StatValue> get_my_features_as_list() const override;
};
