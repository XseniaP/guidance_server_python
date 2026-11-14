#include "utils.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <filesystem>
#include "node.h"

namespace fs = std::filesystem;

std::pair<std::vector<std::vector<double>>, std::unordered_map<std::string, int>>
read_matching_matrix(const fs::path& file_path) {
    std::unordered_map<std::string, int> codes_dict_to_inx;
    std::vector<std::vector<double>> match_matrix;

    std::ifstream infile(file_path);
    
    if (!infile.is_open()) {
        std::cerr << "Cannot open file: [" << file_path << "]" << std::endl;
        throw std::runtime_error("Cannot open matrix file: " + file_path.string());
    }

    std::string line;
    int l_inx = 0;
    while (std::getline(infile, line)) {
        // Trim trailing whitespace/carriage return
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' ')) {
            line.pop_back();
        }
        if (line.empty()) continue;

        if (l_inx == 0) {
            std::istringstream iss(line);
            std::string token;
            while (iss >> token) {
                codes_dict_to_inx[token] = static_cast<int>(codes_dict_to_inx.size());
            }
        } else {
            match_matrix.emplace_back();
            std::istringstream iss(line);
            std::string token;
            bool is_first_col = true;
            while (iss >> token) {
                if (!is_first_col) {
                    match_matrix[l_inx - 1].push_back(std::stod(token));
                }
                is_first_col = false;
            }
        }
        l_inx++;
    }

    return {match_matrix, codes_dict_to_inx};
}

int translate_to_matrix_index(char letter, const std::unordered_map<std::string, int>& code_to_index_dict) {
    std::string key(1, letter);
    auto it = code_to_index_dict.find(key);
    if (it != code_to_index_dict.end()) {
        return it->second;
    }
    auto star_it = code_to_index_dict.find("*");
    if (star_it != code_to_index_dict.end()) {
        return star_it->second;
    }
    throw std::runtime_error("Character not found in matrix index: " + key);
}

double calc_p_distance_from_other(const std::string& aligned_seq, const std::string& other_aligned_seq) {
    int changes_count = 0;
    for (size_t i = 0; i < aligned_seq.size(); i++) {
        if (aligned_seq[i] != other_aligned_seq[i]) {
            changes_count++;
        }
    }
    return static_cast<double>(changes_count) / aligned_seq.size();
}

double calc_kimura_distance_from_other(const std::string& aligned_seq, const std::string& other_aligned_seq) {
    double fractional_identity = calc_p_distance_from_other(aligned_seq, other_aligned_seq);
    double kimura_exponent = 1.0 - fractional_identity - 0.2 * fractional_identity * fractional_identity;
    if (kimura_exponent < 0) {
        return 2.0;
    }
    return -std::log(kimura_exponent);
}

double calc_percentile(std::vector<double>& values, int percentile) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    double index = percentile / 100.0 * (static_cast<double>(values.size()) - 1.0);
    int lower = static_cast<int>(std::floor(index));
    int upper = static_cast<int>(std::ceil(index));
    if (lower == upper) {
        return values[lower];
    }
    double frac = index - lower;
    return values[lower] * (1.0 - frac) + values[upper] * frac;
}

double calc_percentile_int(std::vector<int>& values, int percentile) {
    if (values.empty()) return 0.0;
    std::sort(values.begin(), values.end());
    double index = percentile / 100.0 * (static_cast<double>(values.size()) - 1.0);
    int lower = static_cast<int>(std::floor(index));
    int upper = static_cast<int>(std::ceil(index));
    if (lower == upper) {
        return static_cast<double>(values[lower]);
    }
    double frac = index - lower;
    return values[lower] * (1.0 - frac) + values[upper] * frac;
}

std::vector<Node*> get_raw_pointers_from_unique(const std::vector<std::unique_ptr<Node>>& all_nodes) {
    std::vector<Node*> raw_nodes;
    raw_nodes.reserve(all_nodes.size());
    for (const auto& ptr : all_nodes) {
        raw_nodes.push_back(ptr.get());
    }
    return raw_nodes;
}

string get_model_name_suffix(const string& model_name, double go_val, double ge_val) {
    return "_" + model_name +
        "_GO_" + format_double(go_val) +
        "_GE_" + format_double(ge_val);
}

string format_double(double x) {
    ostringstream oss;
    oss << setprecision(15) << x;  // enough precision, no forced zeros
    return oss.str();
}
