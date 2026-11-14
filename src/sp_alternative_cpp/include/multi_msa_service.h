#pragma once

#include <string>
#include <vector>
#include "config.h"

struct FileNames {
    std::string true_file_name;
    std::string true_tree_file_name;
    std::vector<std::string> other_file_names;
};

FileNames get_file_names_ordered(const std::vector<std::string>& file_names);

void multiple_msa_calc_features_and_labels(const Configuration& config);
