#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <utility>
#include <filesystem>
#include "node.h"
#include "enums.h"


namespace fs = std::filesystem;

std::pair<std::vector<std::vector<double>>, std::unordered_map<std::string, int>>
read_matching_matrix(const fs::path& file_path);

int translate_to_matrix_index(char letter, const std::unordered_map<std::string, int>& code_to_index_dict);

double calc_p_distance_from_other(const std::string& aligned_seq, const std::string& other_aligned_seq);

double calc_kimura_distance_from_other(const std::string& aligned_seq, const std::string& other_aligned_seq);

double calc_percentile(std::vector<double>& values, int percentile);
double calc_percentile_int(std::vector<int>& values, int percentile);
std::vector<Node*> get_raw_pointers_from_unique(const std::vector<std::unique_ptr<Node>>& all_nodes);

bool has_any(const auto& set, std::initializer_list<StatsOutput> values) {
    for (auto v : values) {
        if (set.count(v)) return true;
    }
    return false;
}

string get_model_name_suffix(const std::string& model_name, double go_val, double ge_val);
string format_double(double x);