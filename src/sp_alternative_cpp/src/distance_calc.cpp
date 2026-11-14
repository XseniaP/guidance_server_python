#include "distance_calc.h"
#include <algorithm>
#include <numeric>
#include <iterator>
#include <unordered_map>
#include <iostream>

std::vector<std::string> translate_seq_h(const std::string& sequence, int seq_index, DistanceType distance_type) {
    std::vector<std::string> res;
    int s_count = 1;
    for (char c : sequence) {
        if (c != '-') {
            res.push_back("S^" + std::to_string(seq_index) + "_" + std::to_string(s_count));
            s_count++;
        } else {
            if (distance_type == DistanceType::D_POS) {
                res.push_back("G^" + std::to_string(seq_index) + "_" + std::to_string(s_count - 1));
            } else if (distance_type == DistanceType::D_SEQ) {
                res.push_back("G^" + std::to_string(seq_index));
            } else { // D_SSP
                res.push_back("-");
            }
        }
    }
    return res;
}

std::vector<std::vector<std::string>> translate_profile_naming(const std::vector<std::string>& profile, DistanceType distance_type) {
    std::vector<std::vector<std::string>> res;
    for (int index = 0; index < static_cast<int>(profile.size()); index++) {
        res.push_back(translate_seq_h(profile[index], index + 1, distance_type));
    }
    return res;
}

std::vector<std::string> get_column(const std::vector<std::vector<std::string>>& profile, int col_index) {
    std::vector<std::string> res;
    for (const auto& s : profile) {
        res.push_back(s[col_index]);
    }
    return res;
}

optional<std::set<std::string>> get_place_h(const std::vector<std::string>& column, int seq_index) {
    std::vector<std::string> col = column;
    std::string self_item = col[seq_index];
    col.erase(col.begin() + seq_index);

    if (!self_item.empty() && (self_item[0] == 'G' || self_item[0] == '-')) {
        return std::nullopt;
    }
    return std::set<std::string>(col.begin(), col.end());
}

double get_place_d(const std::set<std::string>& set_a, const std::set<std::string>& set_b, DistanceType distance_type) {
    if (distance_type == DistanceType::D_SSP) {
        std::set<std::string> intersection;
        std::set_intersection(set_a.begin(), set_a.end(), set_b.begin(), set_b.end(),
                              std::inserter(intersection, intersection.begin()));
        std::set<std::string> union_set;
        std::set_union(set_a.begin(), set_a.end(), set_b.begin(), set_b.end(),
                       std::inserter(union_set, union_set.begin()));
        return static_cast<double>(intersection.size()) / union_set.size();
    } else {
        // Symmetric difference
        std::set<std::string> diff_ab, diff_ba;
        std::set_difference(set_a.begin(), set_a.end(), set_b.begin(), set_b.end(),
                            std::inserter(diff_ab, diff_ab.begin()));
        std::set_difference(set_b.begin(), set_b.end(), set_a.begin(), set_a.end(),
                            std::inserter(diff_ba, diff_ba.begin()));
        return static_cast<double>(diff_ab.size() + diff_ba.size()) / (set_a.size() + set_b.size());
    }
}

std::vector<std::vector<std::set<std::string>>> create_h_table(const std::vector<std::vector<std::string>>& profile_naming, DistanceType distance_type) {
    int seq_len = static_cast<int>(profile_naming[0].size());
    int seq_count = static_cast<int>(profile_naming.size());
    std::vector<std::vector<std::set<std::string>>> h_table(seq_count);

    for (int col_index = 0; col_index < seq_len; col_index++) {
        auto col = get_column(profile_naming, col_index);
        for (int i = 0; i < seq_count; i++) {
            auto h = get_place_h(col, i);
            if (h.has_value()) {
                if (distance_type == DistanceType::D_SSP) {
                    // Remove '-' from set
                    std::set<std::string> filtered = h.value();
                    filtered.erase("-");
                    h_table[i].push_back(filtered);
                } else {
                    h_table[i].push_back(h.value());
                }
            }
        }
    }
    return h_table;
}

double compute_distance(const std::vector<std::string>& profile_a, const std::vector<std::string>& profile_b, DistanceType distance_type) {
    auto profile_a_naming = translate_profile_naming(profile_a, distance_type);
    auto profile_b_naming = translate_profile_naming(profile_b, distance_type);
    int seq_count = static_cast<int>(profile_a.size());
    auto profile_a_h = create_h_table(profile_a_naming, distance_type);
    auto profile_b_h = create_h_table(profile_b_naming, distance_type);

    std::vector<double> d_list;
    double numerator_sum = 0;
    double denominator_sum = 0;

    for (int i = 0; i < seq_count; i++) {
        for (int j = 0; j < static_cast<int>(profile_a_h[i].size()); j++) {
            const auto& h_a_i = profile_a_h[i][j];
            const auto& h_b_i = profile_b_h[i][j];
            if (!h_a_i.empty() || !h_b_i.empty()) {
                if (distance_type == DistanceType::D_SSP) {
                    std::set<std::string> intersection, union_set;
                    std::set_intersection(h_a_i.begin(), h_a_i.end(), h_b_i.begin(), h_b_i.end(),
                                          std::inserter(intersection, intersection.begin()));
                    std::set_union(h_a_i.begin(), h_a_i.end(), h_b_i.begin(), h_b_i.end(),
                                   std::inserter(union_set, union_set.begin()));
                    numerator_sum += intersection.size();
                    denominator_sum += union_set.size();
                } else {
                    double d_i_j = get_place_d(h_a_i, h_b_i, distance_type);
                    d_list.push_back(d_i_j);
                }
            } else {
                d_list.push_back(0);
            }
        }
    }

    if (distance_type == DistanceType::D_SSP) {
        if (denominator_sum == 0) return -1;
        return 1.0 - (numerator_sum / denominator_sum);
    } else if (!d_list.empty()) {
        return std::accumulate(d_list.begin(), d_list.end(), 0.0) / d_list.size();
    }
    return -1;
}

double compute_eff_d_seq(const std::vector<std::string>& profile_a, const std::vector<std::string>& profile_b) {
    const int rows_num = profile_a.size();
	const int cols_a_num = profile_a[0].size();
	const int cols_b_num = profile_b[0].size();
    vector<vector<int>> a_vectors(cols_a_num, std::vector<int>(rows_num));
    vector<vector<int>> b_vectors(cols_b_num, std::vector<int>(rows_num));

    vector<long> map_a;
    vector<long> map_b;

    fill_d_seq_vectors(profile_a, a_vectors, map_a, rows_num, cols_a_num);
    fill_d_seq_vectors(profile_b, b_vectors, map_b, rows_num, cols_b_num);

    std::vector<int> counts(cols_a_num * cols_b_num, 0);

	int total_count = map_a.size();
    for (int i = 0; i < total_count; ++i) {
        counts[map_a[i] * cols_b_num + map_b[i]] += 1;
	}

    double total_distance = 0;
    for (int i = 0; i < cols_a_num; ++i) {
        for (int j = 0; j < cols_b_num; ++j) {
            if (counts[i * cols_b_num + j] > 0) {
                double distance = vectors_distance(a_vectors[i], b_vectors[j]) ;
                total_distance += distance * counts[i * cols_b_num + j];
            }
        }
    }
	return total_distance / total_count;
}

void fill_d_seq_vectors(const vector<string>& msa_a, vector<vector<int>>& vectors_list, vector<long>& v_map, int rows_num, int cols_num) {
    for (int i = 0; i < rows_num; ++i) {
        int char_count = 0;
        for (int j = 0; j < cols_num; ++j) {
            if (msa_a[i][j] == '-') {
                vectors_list[j][i] = -1;
            }
            else {
                char_count += 1;
                vectors_list[j][i] = char_count;
                v_map.push_back(j);
            }
        }
    }
}

double vectors_distance(const vector<int>& a, const vector<int>& b) {
    int intersection_count = 0;
    for (int i = 0; i < a.size(); ++i) {
        if (a[i] != b[i]) {
            intersection_count += 1;
        }
    }
    return double(intersection_count) / (a.size() - 1);
}

string key_from_char_counts(int a, int b) {
	return to_string(min(a, b)) + "_" + to_string(max(a, b));
}

double compute_distance_from_known_msa(const std::vector<std::string>& profile_a, const std::vector<std::vector<std::set<std::string>>> &profile_b_h, DistanceType distance_type) {
    auto profile_a_naming = translate_profile_naming(profile_a, distance_type);
    int seq_count = static_cast<int>(profile_a.size());
    auto profile_a_h = create_h_table(profile_a_naming, distance_type);
     
    std::vector<double> d_list;
    double numerator_sum = 0;
    double denominator_sum = 0;

    for (int i = 0; i < seq_count; i++) {
        for (int j = 0; j < static_cast<int>(profile_a_h[i].size()); j++) {
            const auto& h_a_i = profile_a_h[i][j];
            const auto& h_b_i = profile_b_h[i][j];
            if (!h_a_i.empty() || !h_b_i.empty()) {
                if (distance_type == DistanceType::D_SSP) {
                    std::set<std::string> intersection, union_set;
                    std::set_intersection(h_a_i.begin(), h_a_i.end(), h_b_i.begin(), h_b_i.end(),
                        std::inserter(intersection, intersection.begin()));
                    std::set_union(h_a_i.begin(), h_a_i.end(), h_b_i.begin(), h_b_i.end(),
                        std::inserter(union_set, union_set.begin()));
                    numerator_sum += intersection.size();
                    denominator_sum += union_set.size();
                }
                else {
                    double d_i_j = get_place_d(h_a_i, h_b_i, distance_type);
                    d_list.push_back(d_i_j);
                }
            }
            else {
                d_list.push_back(0);
            }
        }
    }

    if (distance_type == DistanceType::D_SSP) {
        if (denominator_sum == 0) return -1;
        return 1.0 - (numerator_sum / denominator_sum);
    }
    else if (!d_list.empty()) {
        return std::accumulate(d_list.begin(), d_list.end(), 0.0) / d_list.size();
    }
    return -1;
}

std::vector<std::vector<std::set<std::string>>> compute_msa_dist_h(const std::vector<std::string>& msa, DistanceType distance_type) {
    auto profile_naming = translate_profile_naming(msa, distance_type);
	return create_h_table(profile_naming, distance_type);
}

double compute_eff_d_seq_from_true(const std::vector<std::string>& msa, const vector<vector<int>>& true_msa_vectors, const vector<long>& true_map) {
    const int rows_num = msa.size();
    const int cols_num = msa[0].size();
    vector<vector<int>> msa_vectors(cols_num, std::vector<int>(rows_num));
    vector<long> msa_map;

    cout << "computing d_seq..." << endl;
    fill_d_seq_vectors(msa, msa_vectors, msa_map, rows_num, cols_num);
    //cout << "End filling msa vectors for efficient d_seq calculation..." << endl;
	int true_msa_cols_num = true_msa_vectors.size();
    std::vector<long> counts(cols_num * true_msa_cols_num, 0);

    //cout << "Filling total count..." << endl;
    int total_count = msa_map.size();
    for (int i = 0; i < total_count; ++i) {
        counts[msa_map[i] * true_msa_cols_num + true_map[i]] += 1;
    }
    //cout << "End filling total count..." << endl;

    double total_distance = 0;
    for (int i = 0; i < cols_num; ++i) {
        //cout << "processing column..." << i<< endl;
        for (int j = 0; j < true_msa_cols_num; ++j) {
            if (counts[i * true_msa_cols_num + j] > 0) {
                double distance = vectors_distance(msa_vectors[i], true_msa_vectors[j]);
                total_distance += distance * counts[i * true_msa_cols_num + j];
            }
        }
    }
    return total_distance / total_count;
}

