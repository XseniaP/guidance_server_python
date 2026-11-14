#include "sop_stats.h"
#include "utils.h"

SopStats::SopStats(const std::string& code, int taxa_num, int msa_length)
    : BasicStats(code, taxa_num, msa_length,
                  {"code", "sp_match_count", "sp_match_count_norm", "sp_mismatch_count", "sp_mismatch_count_norm",
                   "sp_go", "sp_go_norm", "sp_ge", "sp_ge_norm",
                   "sp_match", "sp_match_norm", "sp_mismatch", "sp_mismatch_norm", "sp", "sp_norm"}),
      sp_match_count(-1), sp_match_count_norm(-1),
      sp_mismatch_count(-1), sp_mismatch_count_norm(-1),
      sp_go(-1), sp_go_norm(-1),
      sp_ge(-1), sp_ge_norm(-1),
      sp_match(-1), sp_match_norm(-1),
      sp_mismatch(-1), sp_mismatch_norm(-1),
      sp(-1), sp_norm(-1),
      
      model_agnostic_col_names{ "sp_match_count", "sp_match_count_norm", "sp_mismatch_count","sp_mismatch_count_norm",
                                "sp_go", "sp_go_norm", "sp_ge", "sp_ge_norm" },
      gaps_agnostic_col_names{ "sp_match", "sp_match_norm", "sp_mismatch", "sp_mismatch_norm" } {}

void SopStats::set_my_sop_score_parts(const SPScore& sp_score, const std::vector<std::string>& sequences, vector<vector<int>>& subs_matrix_counts, 
                                      const set<StatsOutput>& statsOutput, const int iteration) {
    bool is_filling_substitutions_matrix = false;
    if (iteration == 0 && has_any(statsOutput, {
        StatsOutput::ALL,
        StatsOutput::SUBS_MATRIX,
        }))
    {
        is_filling_substitutions_matrix = true;
    }
    auto parts = sp_score.compute_efficient_sp_parts(sequences, subs_matrix_counts, is_filling_substitutions_matrix);
    double number_of_pairs_with_msa_length = taxa_num * (taxa_num - 1) / 2.0 * msa_length;

    sp_match_count = parts.sp_match_count;
    sp_match_count_norm = parts.sp_match_count / number_of_pairs_with_msa_length;
    sp_mismatch_count = parts.sp_mismatch_count;
    sp_mismatch_count_norm = parts.sp_mismatch_count / number_of_pairs_with_msa_length;
    sp_go = parts.go_count;
    sp_go_norm = parts.go_count / number_of_pairs_with_msa_length;
    sp_ge = parts.ge_count;
    sp_ge_norm = parts.ge_count / number_of_pairs_with_msa_length;
    sp_match = parts.sp_match_score;
    sp_match_norm = sp_match / number_of_pairs_with_msa_length;
    sp_mismatch = parts.sp_mismatch_score;
    sp_mismatch_norm = sp_mismatch / number_of_pairs_with_msa_length;
    

    sp = parts.sp_match_score + parts.sp_mismatch_score + parts.go_score + parts.sp_score_gap_e;
    sp_norm = sp / number_of_pairs_with_msa_length;
}

std::vector<std::string> SopStats::get_ordered_col_names_with_model(const std::string& model_name,
                                                                      double go_val, double ge_val) const {
    std::vector<std::string> col_names;
    std::string suffix = get_model_name_suffix(model_name, go_val, ge_val);
    for (const auto& col_name : ordered_col_names) {
        if (model_agnostic_col_names.count(col_name)) {
                col_names.push_back(col_name);
        } else if (gaps_agnostic_col_names.count(col_name)) {
            col_names.push_back(col_name + "_" + model_name);
        } else {
            col_names.push_back(col_name + suffix);
        }
    }
    return col_names;
}

std::vector<StatValue> SopStats::get_my_features_as_list() const {
    return {
        StatValue(code),
        StatValue(sp_match_count), StatValue(sp_match_count_norm),
        StatValue(sp_mismatch_count), StatValue(sp_mismatch_count_norm),
        StatValue(sp_go), StatValue(sp_go_norm),
        StatValue(sp_ge), StatValue(sp_ge_norm),
        StatValue(sp_match), StatValue(sp_match_norm),
        StatValue(sp_mismatch), StatValue(sp_mismatch_norm),
        StatValue(sp), StatValue(sp_norm)
    };
}
