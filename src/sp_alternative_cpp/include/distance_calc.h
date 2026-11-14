#pragma once
#include "enums.h"
#include <vector>
#include <string>
#include <set>
#include <optional>

std::vector<std::string> translate_seq_h(const std::string& sequence, int seq_index, DistanceType distance_type);
std::vector<std::vector<std::string>> translate_profile_naming(const std::vector<std::string>& profile, DistanceType distance_type);
std::vector<std::string> get_column(const std::vector<std::vector<std::string>>& profile, int col_index);
optional<std::set<std::string>> get_place_h(const std::vector<std::string>& column, int seq_index);
double get_place_d(const std::set<std::string>& set_a, const std::set<std::string>& set_b, DistanceType distance_type);
std::vector<std::vector<std::set<std::string>>> create_h_table(const std::vector<std::vector<std::string>>& profile_naming, DistanceType distance_type);
double compute_distance(const std::vector<std::string>& profile_a, const std::vector<std::string>& profile_b, DistanceType distance_type);
double compute_distance_from_known_msa(const std::vector<std::string>& profile_a, const std::vector<std::vector<std::set<std::string>>> & profile_b_h, DistanceType distance_type);
std::vector<std::vector<std::set<std::string>>> compute_msa_dist_h(const std::vector<std::string>& msa, DistanceType distance_type);
double compute_eff_d_seq(const std::vector<std::string>& profile_a, const std::vector<std::string>& profile_b);
void fill_d_seq_vectors(const std::vector<std::string>& msa_a, vector<vector<int>>& vectors_list, vector<long>& v_map, int rows_num, int cols_num);
double vectors_distance(const vector<int>& a, const vector<int>& b);
double compute_eff_d_seq_from_true(const std::vector<std::string>& msa, const vector<vector<int>>& true_msa_vectors, const vector<long>& true_map);

