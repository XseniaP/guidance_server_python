#include "msa.h"
#include "config.h"
#include "sp_score.h"
#include "unrooted_tree.h"
#include "enums.h"
#include <iostream>
#include <filesystem>
#include <algorithm>
#include <chrono>
#include <vector>
#include <string>
#include "distance_calc.h"

namespace fs = std::filesystem;

struct FileNames {
    std::string true_file_name;
    std::string true_tree_file_name;
    std::vector<std::string> other_file_names;
};

FileNames get_file_names_ordered(const std::vector<std::string>& file_names) {
    FileNames result;
    for (const auto& fn : file_names) {
        if (fn.find("_TRUE.") != std::string::npos) {
            result.true_file_name = fn;
        } else if (fn.find("_true_tree.") != std::string::npos) {
            result.true_tree_file_name = fn;
        } else {
            result.other_file_names.push_back(fn);
        }
    }
    return result;
}

void multiple_msa_calc_features_and_labels(const Configuration& config) {
    auto start = std::chrono::steady_clock::now();

    fs::path comparison_dir(config.input_files_dir_path);

    std::vector<SPScore> sp_models;
    sp_models.reserve(config.models.size());

    for (const auto& m : config.models) {
        fs::path matrix_path = fs::path(config.matrix_dir_path) / m.matrix_file_name;
        matrix_path.replace_extension(".txt");
        sp_models.emplace_back(m, matrix_path);
    }
    cout << "End subs matrices: " << endl;

    fs::path output_dir_path = fs::path(config.output_file_dir_path);
    fs::create_directories(output_dir_path);

    for (const auto& dir_entry : fs::directory_iterator(comparison_dir)) {
        if (!dir_entry.is_directory()) continue;
        std::string dir_name = dir_entry.path().filename().string();
        fs::path dir_path = dir_entry.path();

        // Collect file names
        std::vector<std::string> files;
        for (const auto& f : fs::directory_iterator(dir_path)) {
            if (f.is_regular_file()) {
                files.push_back(f.path().filename().string());
            }
        }

        auto ordered = get_file_names_ordered(files);
        cout << "Starting to work on an msa: " << endl;
        MSA true_msa(dir_name);
        if (!ordered.true_file_name.empty()) {
            true_msa.read_from_fasta(dir_path / ordered.true_file_name);
        }
        //std::vector<std::vector<std::set<std::string>>> true_msa_dist_h = compute_msa_dist_h(true_msa.sequences, DistanceType::D_SEQ);


        bool is_init_files = true;
        unordered_map<string, int> col_names_dict;
        for (const auto& inferred_file_name : ordered.other_file_names) {
            std::string msa_name = inferred_file_name;
            std::cout << "=================================="<< msa_name << "===================================" << std::endl;

            MSA inferred_msa(msa_name);
            inferred_msa.read_from_fasta(dir_path / inferred_file_name);
            inferred_msa.order_sequences(true_msa.seq_names);
            inferred_msa.calc_and_print_stats(true_msa, config, sp_models, output_dir_path, dir_name,
                                               true_msa.tree.get(), is_init_files, col_names_dict);
            //string newick_str = inferred_msa.tree->print_newick();
            is_init_files = false;
        }

        // Also process the true MSA
        true_msa.calc_and_print_stats(true_msa, config, sp_models, output_dir_path, dir_name,
                                       true_msa.tree.get(), is_init_files, col_names_dict);
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Done - Total duration: "
              << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
}
