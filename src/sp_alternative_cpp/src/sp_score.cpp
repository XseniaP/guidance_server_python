#include "sp_score.h"
#include "utils.h"
#include <numeric>
#include <algorithm>
#include <filesystem>

using namespace std;

namespace fs = std::filesystem;



SPScore::SPScore(const EvoModel& evo_model, const fs::path& matrix_path) {
    auto [wm, cdi] = read_matching_matrix(matrix_path);
    w_matrix = wm;
    code_to_index_dict = cdi;
    go_cost = evo_model.go_cost;
    ge_cost = evo_model.ge_cost;
    model_name = evo_model.name;
}

std::vector<double> SPScore::compute_naive_sp_score(const std::vector<std::string>& profile,
                                                      const std::vector<std::vector<double>>* seq_w_options_ptr) const {
    std::vector<std::vector<double>> default_w;
    const std::vector<std::vector<double>>* seq_w_options;
    if (seq_w_options_ptr == nullptr) {
        default_w.push_back(std::vector<double>(profile.size(), 1.0));
        seq_w_options = &default_w;
    } else {
        seq_w_options = seq_w_options_ptr;
    }

    int weight_options_count = static_cast<int>(seq_w_options->size());
    if (profile.empty()) return std::vector<double>(weight_options_count, 0.0);

    int seq_len = static_cast<int>(profile[0].size());
    std::vector<double> sp_score_subs(weight_options_count, 0.0);
    std::vector<double> sp_score_gaps(weight_options_count, 0.0);

    for (int i = 0; i < static_cast<int>(profile.size()); i++) {
        for (int j = i + 1; j < static_cast<int>(profile.size()); j++) {
            std::vector<char> clean_seq_i, clean_seq_j;
            std::vector<double> seq_weights_mult(weight_options_count);
            for (int w = 0; w < weight_options_count; w++) {
                seq_weights_mult[w] = (*seq_w_options)[w][i] * (*seq_w_options)[w][j];
            }
            for (int k = 0; k < seq_len; k++) {
                if (!(profile[i][k] == '-' && profile[j][k] == '-')) {
                    clean_seq_i.push_back(profile[i][k]);
                    clean_seq_j.push_back(profile[j][k]);
                }
                if (profile[i][k] != '-' && profile[j][k] != '-') {
                    int s = subst(profile[i][k], profile[j][k]);
                    for (int w = 0; w < weight_options_count; w++) {
                        sp_score_subs[w] += s * seq_weights_mult[w];
                    }
                }
            }
            auto gaps_i = compute_gap_intervals(clean_seq_i);
            auto gaps_j = compute_gap_intervals(clean_seq_j);
            std::vector<GapInterval> all_gaps;
            all_gaps.insert(all_gaps.end(), gaps_i.begin(), gaps_i.end());
            all_gaps.insert(all_gaps.end(), gaps_j.begin(), gaps_j.end());
            for (const auto& gi : all_gaps) {
                for (int w = 0; w < weight_options_count; w++) {
                    sp_score_gaps[w] += gi.g_cost(go_cost, ge_cost) * seq_weights_mult[w];
                }
            }
        }
    }

    std::vector<double> result(weight_options_count);
    for (int w = 0; w < weight_options_count; w++) {
        result[w] = sp_score_subs[w] + sp_score_gaps[w];
    }
    return result;
}

std::vector<GapInterval> SPScore::compute_gap_intervals(const std::vector<char>& seq_i, double seq_w) {
    int seq_len = static_cast<int>(seq_i.size());
    std::vector<GapInterval> gap_intervals_list;
    GapInterval gap_interval(seq_w);

    for (int k = 0; k < seq_len; k++) {
        if (seq_i[k] == '-' && gap_interval.is_empty()) {
            gap_interval.set_start(k);
        }
        if (seq_i[k] != '-' && !gap_interval.is_empty()) {
            gap_interval.set_end(k - 1);
            gap_intervals_list.push_back(gap_interval.copy_me());
            gap_interval = GapInterval(seq_w);
        }
    }
    if (!gap_interval.is_empty()) {
        gap_interval.set_end(seq_len - 1);
        gap_intervals_list.push_back(gap_interval.copy_me());
    }
    return gap_intervals_list;
}

SPScore::SpSAndGe SPScore::compute_sp_s_and_sp_ge(const std::vector<std::string>& profile, vector<vector<int>> & substitutions_matrix, bool is_filling_substitutions_matrix) const {
    int options_count = static_cast<int>(w_matrix[0].size());
    int seq_len = static_cast<int>(profile[0].size());
    double sp_match_score = 0;
    double sp_mismatch_score = 0;
    int ge_count = 0;
    int sp_match_count = 0;
    int sp_mismatch_count = 0;

    if (is_filling_substitutions_matrix) {
        substitutions_matrix.resize(options_count);
        for (int i=0; i < options_count; i++) {
            substitutions_matrix[i].resize(options_count, 0);
		}
	}

    for (int k = 0; k < seq_len; k++) {
        std::vector<int> histo_count(options_count + 1, 0);
        for (int i = 0; i < static_cast<int>(profile.size()); i++) {
            char ch = profile[i][k];
            if (ch == '-') {
                histo_count[options_count]++;
            } else {
                int char_index = translate_to_matrix_index(ch, code_to_index_dict);
                histo_count[char_index]++;
            }
        }

        for (int i = 0; i < options_count; i++) {
            if (histo_count[i] != 0) {
                sp_match_score += static_cast<double>(w_matrix[i][i]) *
                                  histo_count[i] * (histo_count[i] - 1) / 2.0;
                sp_match_count += histo_count[i] * (histo_count[i] - 1) / 2;
                for (int j = i + 1; j < options_count; j++) {
                    if (histo_count[j] != 0) {
                        sp_mismatch_score += static_cast<double>(w_matrix[i][j]) *
                            histo_count[i] * histo_count[j];
                        sp_mismatch_count += histo_count[i] * histo_count[j];
                        if (is_filling_substitutions_matrix) {
							substitutions_matrix[i][j] += histo_count[i] * histo_count[j];
                        }
                    }
                }
            }
        }
        if (histo_count[options_count] > 0) {
            ge_count += (static_cast<int>(profile.size()) - histo_count[options_count]) *
                        histo_count[options_count];
        }
    }

    return {sp_match_score, sp_mismatch_score, ge_count * ge_cost,
            sp_match_count, sp_mismatch_count, ge_count};
}

int SPScore::subst(char a, char b) const {
    int ai = translate_to_matrix_index(a, code_to_index_dict);
    int bi = translate_to_matrix_index(b, code_to_index_dict);
    return w_matrix[ai][bi];
}

SPScore::SpGapOpen SPScore::compute_sp_gap_open(const std::vector<std::string>& profile) const {
    if (profile.empty()) return {0, 0};

    int seq_len = static_cast<int>(profile[0].size());
    int n = static_cast<int>(profile.size());
    std::vector<int> nb_open_gap(seq_len, 0);
    std::vector<std::vector<GapInterval>> gap_closing(seq_len);

    for (const auto& seq_i : profile) {
        std::vector<char> seq_chars(seq_i.begin(), seq_i.end());
        auto gap_intervals_list = compute_gap_intervals(seq_chars);
        for (const auto& gi : gap_intervals_list) {
            for (int k = gi.start; k <= gi.end; k++) {
                nb_open_gap[k]++;
            }
            gap_closing[gi.end].push_back(gi);
        }
    }

    double sp_gp_open = 0;
    int sp_gpo_count = 0;

    for (int i = 0; i < seq_len; i++) {
        for (const auto& gi : gap_closing[i]) {
            int gpo_count = n - nb_open_gap[gi.start];
            sp_gp_open += gpo_count * go_cost;
            sp_gpo_count += gpo_count;
        }
        for (const auto& gi : gap_closing[i]) {
            for (int k = gi.start; k <= gi.end; k++) {
                nb_open_gap[k]--;
            }
        }
    }

    return {sp_gp_open, sp_gpo_count};
}

double SPScore::compute_efficient_sp(const std::vector<std::string>& profile, bool is_filling_substitutions_matrix) const {
    vector<vector<int>> subs_matrix_counts;
    auto [sp_match, sp_mismatch, sp_ge, mc, mmc, gc] = compute_sp_s_and_sp_ge(profile, subs_matrix_counts, is_filling_substitutions_matrix);
    auto [go_score, gpo_count] = compute_sp_gap_open(profile);
    return sp_match + sp_mismatch + sp_ge + go_score;
}

SPScore::EfficientSpParts SPScore::compute_efficient_sp_parts(const std::vector<std::string>& profile, vector<vector<int>>& subs_matrix_counts, 
                                                              bool is_filling_substitutions_matrix) const {
    auto [sp_match, sp_mismatch, sp_ge, mc, mmc, gc] = compute_sp_s_and_sp_ge(profile, subs_matrix_counts, is_filling_substitutions_matrix);
    auto [go_score, go_count] = compute_sp_gap_open(profile);
    return {sp_match, sp_mismatch, go_score, sp_ge, mc, mmc, go_count, gc, subs_matrix_counts};
}

double SPScore::get_pair_score(int i, int j) const {
    int min_i = std::min(i, j);
    int max_j = std::max(i, j);
    if (max_j < static_cast<int>(w_matrix.size())) {
        return w_matrix[min_i][max_j];
    }
    return ge_cost;
}

std::vector<double> SPScore::compute_w_sp_s_and_sp_ge(const std::vector<std::string>& alignment,
                                                        const std::vector<double>& seq_w) const {
    int options_count = static_cast<int>(w_matrix[0].size()) + 1;
    int seq_len = static_cast<int>(alignment[0].size());
    std::vector<double> sp_subs_and_ge_score(seq_len, 0.0);

    for (int col = 0; col < seq_len; col++) {
        struct Hist { int count; double sq_sum; double sum; };
        std::vector<Hist> histo(options_count, {0, 0.0, 0.0});

        for (int i = 0; i < static_cast<int>(alignment.size()); i++) {
            char ch = alignment[i][col];
            int char_index = (ch == '-') ? (options_count - 1)
                              : translate_to_matrix_index(ch, code_to_index_dict);
            histo[char_index].count++;
            histo[char_index].sum += seq_w[i];
            histo[char_index].sq_sum += seq_w[i] * seq_w[i];
        }

        for (int ci = 0; ci < options_count - 1; ci++) {
            if (histo[ci].count != 0) {
                sp_subs_and_ge_score[col] +=
                    static_cast<double>(get_pair_score(ci, ci) *
                    (histo[ci].sum * histo[ci].sum - histo[ci].sq_sum) / 2.0);
                for (int j = ci + 1; j < options_count; j++) {
                    if (histo[j].count != 0) {
                        sp_subs_and_ge_score[col] +=
                            static_cast<double>(get_pair_score(ci, j) *
                            histo[ci].sum * histo[j].sum);
                    }
                }
            }
        }
    }

    return sp_subs_and_ge_score;
}

double SPScore::compute_w_sp_gap_open(const std::vector<std::string>& alignment,
                                        const std::vector<double>& w) const {
    if (alignment.empty()) return 0;

    int seq_len = static_cast<int>(alignment[0].size());
    int n = static_cast<int>(alignment.size());

    struct GapInfo { int count; double sum; };
    std::vector<GapInfo> nb_open_gap(seq_len, {0, 0.0});
    std::vector<std::vector<GapInterval>> gap_closing(seq_len);
    double sum_w = std::accumulate(w.begin(), w.end(), 0.0);

    for (int index = 0; index < static_cast<int>(alignment.size()); index++) {
        std::vector<char> seq_chars(alignment[index].begin(), alignment[index].end());
        auto gap_intervals_list = compute_gap_intervals(seq_chars, w[index]);
        for (const auto& gi : gap_intervals_list) {
            for (int k = gi.start; k <= gi.end; k++) {
                nb_open_gap[k].count++;
                nb_open_gap[k].sum += gi.w;
            }
            gap_closing[gi.end].push_back(gi);
        }
    }

    double sp_gp_open = 0;
    for (int col = 0; col < seq_len; col++) {
        for (const auto& gi : gap_closing[col]) {
            double gpo_w = sum_w - nb_open_gap[gi.start].sum;
            sp_gp_open += gpo_w * gi.w * go_cost;
        }
        for (const auto& gi : gap_closing[col]) {
            for (int k = gi.start; k <= gi.end; k++) {
                nb_open_gap[k].count--;
                nb_open_gap[k].sum -= gi.w;
            }
        }
    }

    return sp_gp_open;
}

double SPScore::compute_efficient_w_sp(const std::vector<std::string>& alignment,
                                         const std::vector<double>& w) const {
    auto sp_s_and_ge_score = compute_w_sp_s_and_sp_ge(alignment, w);
    double sp_gp_open = compute_w_sp_gap_open(alignment, w);
    double sum = std::accumulate(sp_s_and_ge_score.begin(), sp_s_and_ge_score.end(), 0.0);
    return sum + sp_gp_open;
}
