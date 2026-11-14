#pragma once
#include "evo_model.h"
#include "gap_interval.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <tuple>
#include <filesystem>
using namespace std;

namespace fs = std::filesystem;

class SPScore {
public:
    std::vector<std::vector<double>> w_matrix;
    std::unordered_map<std::string, int> code_to_index_dict;
    double go_cost;
    double ge_cost;
    std::string model_name;

    SPScore(const EvoModel& evo_model, const fs::path& matrix_path);

    std::vector<double> compute_naive_sp_score(const std::vector<std::string>& profile,
                                                const std::vector<std::vector<double>>* seq_w_options = nullptr) const;

    struct SpPerCol {
        std::vector<double> subs;
        std::vector<double> gap_o;
        std::vector<double> gap_e;
    };
    SpPerCol compute_naive_sp_score_per_col(const std::vector<std::string>& profile) const;

    static std::vector<GapInterval> compute_gap_intervals(const std::vector<char>& seq_i, double seq_w = 1.0);

    struct SpSAndGe {
        double sp_match_score;
        double sp_mismatch_score;
        double ge_score;
        int sp_match_count;
        int sp_mismatch_count;
        int ge_count;
    };
    SpSAndGe compute_sp_s_and_sp_ge(const std::vector<std::string>& profile, vector<vector<int>>& substitutions_matrix, bool is_filling_substitutions_matrix) const;

    int subst(char a, char b) const;

    struct SpGapOpen {
        double sp_gp_open;
        int sp_gpo_count;
    };
    SpGapOpen compute_sp_gap_open(const std::vector<std::string>& profile) const;

    double compute_efficient_sp(const std::vector<std::string>& profile, bool is_filling_substitutions_matrix) const;

    struct EfficientSpParts {
        double sp_match_score;
        double sp_mismatch_score;
        double go_score;
        double sp_score_gap_e;
        int sp_match_count;
        int sp_mismatch_count;
        int go_count;
        int ge_count;
        vector<vector<int>> subs_matrix_counts;
    };
    EfficientSpParts compute_efficient_sp_parts(const std::vector<std::string>& profile, vector<vector<int>>& subs_matrix_counts, bool is_filling_substitutions_matrix) const;

    double get_pair_score(int i, int j) const;

    std::vector<double> compute_w_sp_s_and_sp_ge(const std::vector<std::string>& alignment,
                                                   const std::vector<double>& seq_w) const;
    double compute_w_sp_gap_open(const std::vector<std::string>& alignment,
                                  const std::vector<double>& w) const;
    double compute_efficient_w_sp(const std::vector<std::string>& alignment,
                                   const std::vector<double>& w) const;
};
