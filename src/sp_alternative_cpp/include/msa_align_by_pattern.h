#pragma once
#include "msa.h"
#include <string>
#include <filesystem>

class MsaAlignByPattern : public MSA {
public:
    MSA pattern_msa;

    MsaAlignByPattern(const std::string& dataset_name,
                       const std::string& original_file_path,
                       const std::string& pattern_file_path,
                       const std::string& output_dir_path);

    void read_data(const std::string& original_file_path, const std::string& pattern_file_path);
    void align_me();
};
