#include "msa.h"
#include "neighbor_joining.h"
#include "node.h"
#include "sop_stats.h"
#include "w_sop_stats.h"
#include "gaps_stats.h"
#include "entropy_stats.h"
#include "kmer_stats.h"
#include "tree_stats.h"
#include "dist_labels_stats.h"
#include "subs_matrix_counter.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <chrono>
#include <sstream>
#include <algorithm>
#include <unordered_map>

MSA::MSA(const std::string& dataset_name)
    : dataset_name(dataset_name) {}

void MSA::add_sequence(const std::string& sequence, const std::string& seq_name) {
    sequences.push_back(sequence);
    seq_names.push_back(seq_name);
}

void MSA::set_sequences(const std::vector<std::string>& seqs, const std::vector<std::string>& names) {
    sequences = seqs;
    seq_names = names;
}

void MSA::read_from_fasta(const std::filesystem::path& file_path) {
    std::ifstream in_file(file_path);
    if (!in_file.is_open()) {
        std::cerr << "Cannot open FASTA file: " << file_path << std::endl;
        return;
    }

    std::string seq;
    std::string seq_name;
    std::string line;

    while (std::getline(in_file, line)) {
        // Trim
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) {
            line.pop_back();
        }
        if (line.empty()) {
            if (!seq.empty()) {
                seq.erase(seq.find_last_not_of(" ") + 1);
                add_sequence(seq, seq_name);
            }
            return;
        }
        if (line[0] == '>') {
            if (!seq.empty()) {
                seq.erase(seq.find_last_not_of(" ") + 1);
                add_sequence(seq, seq_name);
                seq.clear();
            }
            seq_name = line.substr(1);
            seq_name.erase(seq_name.find_last_not_of(" ") + 1);
        } else {
            seq += line;
        }
    }
    if (!seq.empty()) {
        seq.erase(seq.find_last_not_of(" ") + 1);
        add_sequence(seq, seq_name);
    }
}

void MSA::order_sequences(const std::vector<std::string>& ordered_seq_names) {
    std::unordered_map<std::string, int> names_dict;
    for (int i = 0; i < static_cast<int>(seq_names.size()); i++) {
        names_dict[seq_names[i]] = i;
    }
    std::vector<std::string> ordered_seq;
    for (const auto& name : ordered_seq_names) {
        auto it = names_dict.find(name);
        if (it != names_dict.end()) {
            ordered_seq.push_back(sequences[it->second]);
        }
    }
    if (ordered_seq.size() != sequences.size()) { // TODO: continue from here
		throw std::runtime_error("Sequence count mismatch during ordering");
    } else {
        sequences = ordered_seq;
        seq_names = ordered_seq_names;
    }
}

void MSA::build_nj_tree() {
    int n = static_cast<int>(sequences.size());
    std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n, 0.0));
    std::vector<unique_ptr<Node>>nodes;

    std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(nodes);

	cout << "Calculating distance matrix for NJ tree construction..." << endl;
    for (int i = 0; i < n; i++) {
        auto node = std::make_unique<Node>(i, std::set<std::string>{seq_names[i]}, std::vector<int>{}, 0.0);
        node->fill_newick(raw_nodes);
        nodes.push_back(move(node));
        for (int j = i; j < n; j++) {
            double kd = calc_kimura_distance_from_other(sequences[i], sequences[j]);
            distance_matrix[i][j] = kd;
            distance_matrix[j][i] = kd;
        }
    }
    cout << "End calculating distance matrix for NJ tree construction..." << endl;

    NeighborJoining nj(distance_matrix, std::move(nodes));

   
    // Access the calculated tree (tree_res is std::optional)
    if (nj.tree_res.has_value()) {
        tree = std::make_unique<UnrootedTree>(std::move(nj.tree_res.value()));
        tree->print_newick();
    }
}

void MSA::set_tree(UnrootedTree&& t) {
    tree = std::make_unique<UnrootedTree>(std::move(t));
}

int MSA::get_msa_len() const {
    return sequences.empty() ? 0 : static_cast<int>(sequences[0].size());
}

int MSA::get_taxa_num() const {
    return static_cast<int>(sequences.size());
}

void MSA::calc_and_print_stats(const MSA& true_msa, const Configuration& config,
                                 const std::vector<SPScore>& sp_models,
                                 const std::filesystem::path& output_dir_path,
                                 const std::string& batch_name,
                                 const UnrootedTree* true_tree,
                                 bool is_init_file,
                                 unordered_map<string, int>& col_names_dict) {
    vector<vector<StatValue>>data;
    vector<vector<string>>data_names;

    // Basic stats
    {
        std::vector<std::string> basic_cols = {"code", "taxa_num", "msa_length"};
        std::vector<StatValue> basic_vals = {StatValue(dataset_name), StatValue(get_taxa_num()), StatValue(get_msa_len())};
        if (config.is_unified_file) {
            data.push_back(basic_vals);
            data_names.push_back(basic_cols);
        }
        else {
            print_stats_files(basic_vals, batch_name, output_dir_path, "basic_stats", is_init_file, basic_cols);
        }
    }

    DistanceLabelsStats dist_labels_stats(dataset_name, get_taxa_num(), get_msa_len());

    // Distance labels
    if (has_any(config.stats_output, {
        StatsOutput::ALL,
        StatsOutput::ALL_NO_SUBS_MATRIX,
        StatsOutput::DISTANCE_LABELS
        }))
    {
        auto start = std::chrono::steady_clock::now();
        dist_labels_stats.set_my_distance_from_true(sequences, true_msa.sequences);
        if (config.is_unified_file) {
            data.push_back(dist_labels_stats.get_my_features_as_list());
            data_names.push_back(dist_labels_stats.get_ordered_col_names());
        }
        else {
            print_stats_files(dist_labels_stats.get_my_features_as_list(), batch_name, output_dir_path,
                stats_output_to_string(StatsOutput::DISTANCE_LABELS), is_init_file,
                dist_labels_stats.get_ordered_col_names());
        }
        auto end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time for Labels: "
                  << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

    // Entropy
    if (has_any(config.stats_output, {
        StatsOutput::ALL,
        StatsOutput::ALL_NO_SUBS_MATRIX,
        StatsOutput::ENTROPY
        }))
    {
        auto start = std::chrono::steady_clock::now();
        EntropyStats entropy_stats(dataset_name, get_taxa_num(), get_msa_len());
        entropy_stats.calc_entropy(sequences);
        if (config.is_unified_file) {
            data.push_back(entropy_stats.get_my_features_as_list());
            data_names.push_back(entropy_stats.get_ordered_col_names());
        }
        else {
            print_stats_files(entropy_stats.get_my_features_as_list(), batch_name, output_dir_path,
                stats_output_to_string(StatsOutput::ENTROPY), is_init_file,
                entropy_stats.get_ordered_col_names());
        }
        auto end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time for Entropy: "
                  << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

    // Gaps
    if (has_any(config.stats_output, {
        StatsOutput::ALL,
        StatsOutput::ALL_NO_SUBS_MATRIX,
        StatsOutput::GAPS
        }))
    {
        auto start = std::chrono::steady_clock::now();
        GapStats gaps_stats(dataset_name, get_taxa_num(), get_msa_len());
        gaps_stats.calc_gaps_values(sequences);
        if (config.is_unified_file) {
            data.push_back(gaps_stats.get_my_features_as_list());
            data_names.push_back(gaps_stats.get_ordered_col_names());
        }
        else {
            print_stats_files(gaps_stats.get_my_features_as_list(), batch_name, output_dir_path,
                stats_output_to_string(StatsOutput::GAPS), is_init_file,
                gaps_stats.get_ordered_col_names());
        }
        auto end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time for Gaps: "
                  << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

    // K-mer
    if (has_any(config.stats_output, {
        StatsOutput::ALL,
        StatsOutput::ALL_NO_SUBS_MATRIX,
        StatsOutput::K_MER
        }))
    {
        auto start = std::chrono::steady_clock::now();
        for (int k_value : config.k_values) {
            KMerStats kmer_stats(dataset_name, get_taxa_num(), get_msa_len(), k_value);
            kmer_stats.set_k_mer_features(sequences);
            if (config.is_unified_file) {
                data.push_back(kmer_stats.get_my_features_as_list());
                data_names.push_back(kmer_stats.get_ordered_col_names_with_k_value());
            }
            else {
                print_stats_files(kmer_stats.get_my_features_as_list(), batch_name, output_dir_path,
                    stats_output_to_string(StatsOutput::K_MER), is_init_file,
                    kmer_stats.get_ordered_col_names_with_k_value(), "", 0, 0,
                    std::to_string(k_value));
            }
        }
        auto end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time for Kmer: "
                  << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

    // Tree
    if (has_any(config.stats_output, {
        StatsOutput::ALL,
        StatsOutput::ALL_NO_SUBS_MATRIX,
        StatsOutput::TREE
        }))
    {
        auto start = std::chrono::steady_clock::now();
        build_nj_tree();
        TreeStats tree_stats(dataset_name, get_taxa_num(), get_msa_len());
        tree_stats.set_tree_stats(tree->get_branches_lengths_list(), *tree, sequences, seq_names);
        if (config.is_unified_file) {
            data.push_back(tree_stats.get_my_features_as_list());
            data_names.push_back(tree_stats.get_ordered_col_names());
        }
        else {
            print_stats_files(tree_stats.get_my_features_as_list(), batch_name, output_dir_path,
                stats_output_to_string(StatsOutput::TREE), is_init_file,
                tree_stats.get_ordered_col_names());
        }

        auto end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time for Tree: "
                  << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

    // SP/SOP
    vector<vector<int>>subs_matrix_counts;

    if (has_any(config.stats_output, {
            StatsOutput::ALL,
            StatsOutput::ALL_NO_SUBS_MATRIX,
            StatsOutput::SUBS_MATRIX,
            StatsOutput::SP
        }))
    {
        auto start = std::chrono::steady_clock::now();
        for (int i = 0; i < sp_models.size(); ++i) {
            SPScore sp = sp_models[i];
            SopStats sop_stats(dataset_name, get_taxa_num(), get_msa_len());
            sop_stats.set_my_sop_score_parts(sp, sequences, subs_matrix_counts, config.stats_output, i);
            if (config.is_unified_file) {
                data.push_back(sop_stats.get_my_features_as_list());
                data_names.push_back(sop_stats.get_ordered_col_names_with_model(sp.model_name, sp.go_cost, sp.ge_cost));
            }
            else {
                print_stats_files(sop_stats.get_my_features_as_list(), batch_name, output_dir_path,
                    stats_output_to_string(StatsOutput::SP), is_init_file,
                    sop_stats.get_ordered_col_names_with_model(sp.model_name, sp.go_cost, sp.ge_cost),
                    sp.model_name, sp.go_cost, sp.ge_cost);
            }
        }
        auto end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time for regular Sop: "
                  << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }

	// Weighted SP 
    // // TODO: change this, don't root for each model, root once and calculate the weights for each model, then calculate wSOP for each model.
    if (has_any(config.stats_output, {
        StatsOutput::ALL,
        StatsOutput::ALL_NO_SUBS_MATRIX,
        StatsOutput::W_SP
        }))
    {
        auto start = std::chrono::steady_clock::now();
        if (config.stats_output.count(StatsOutput::W_SP)) {
            build_nj_tree();
        }
        for (const auto& sp : sp_models) {
            WSopStats w_sop_stats(dataset_name, get_taxa_num(), get_msa_len());
            w_sop_stats.calc_seq_weights(config.additional_weights, sequences, seq_names, *tree);
            w_sop_stats.calc_w_sp(sequences, sp);
            if (config.is_unified_file) {
                data.push_back(w_sop_stats.get_my_features_as_list());
                data_names.push_back(w_sop_stats.get_ordered_col_names_with_model(sp.model_name, sp.go_cost, sp.ge_cost));
            } else {
                print_stats_files(w_sop_stats.get_my_features_as_list(), batch_name, output_dir_path,
                    stats_output_to_string(StatsOutput::W_SP), is_init_file,
                    w_sop_stats.get_ordered_col_names_with_model(sp.model_name, sp.go_cost, sp.ge_cost),
                    sp.model_name, sp.go_cost, sp.ge_cost);
            }
        }
        auto end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time for wSop: "
                  << std::chrono::duration<double>(end - start).count() << " seconds" << std::endl;
    }
    
    // Subs Matrix Counter
    if (has_any(config.stats_output, {
        StatsOutput::ALL,
        StatsOutput::SUBS_MATRIX,
        }))
    {
        unordered_map<std::string, int> code_to_index_dict = sp_models[0].code_to_index_dict;
        vector<string> labels(subs_matrix_counts[0].size());
        for (const auto& [code, index] : code_to_index_dict) {
            labels[index] = code;
        }
        SubsMatrixCounterStats subs_matrix_counter_stats(dataset_name, get_taxa_num(), get_msa_len(), subs_matrix_counts, labels);
        if (config.is_unified_file) {
            data.push_back(subs_matrix_counter_stats.get_my_features_as_list());
            data_names.push_back(subs_matrix_counter_stats.get_ordered_col_names_with_labels_value(labels));
        }
        else {
            print_stats_files(subs_matrix_counter_stats.get_my_features_as_list(), batch_name, output_dir_path,
                stats_output_to_string(StatsOutput::SUBS_MATRIX), is_init_file, subs_matrix_counter_stats.get_ordered_col_names_with_labels_value(labels), "", 0, 0);
        }
    }
    if (config.is_unified_file) {
        print_unified_stats_files(output_dir_path, batch_name, is_init_file, col_names_dict, data_names, data);
    }
}

void MSA::print_stats_files(const std::vector<StatValue>& stats_data,
                             const string& batch_name,   
                             const std::filesystem::path& output_dir_path,
                             const std::string& file_name,
                             bool is_init_file,
                             const std::vector<std::string>& col_names,
                             const std::string& model_name,
                             double go_val, double ge_val,
                             const std::string& k_value) {
    std::string model_str;
    if (!model_name.empty()) {
        model_str = get_model_name_suffix(model_name, go_val, ge_val);
    }
    std::string k_value_str;
    if (!k_value.empty()) {
        k_value_str = "_K" + k_value;
    }

    std::filesystem::path output_file = output_dir_path / (batch_name + "_" + file_name + model_str + k_value_str + ".csv");

    if (is_init_file) {
        std::ofstream outfile(output_file);
        // Write header
        for (size_t i = 0; i < col_names.size(); i++) {
            if (i > 0) outfile << ",";
            outfile << col_names[i];
        }
        outfile << "\n";
        // Write data
        for (size_t i = 0; i < stats_data.size(); i++) {
            if (i > 0) outfile << ", ";
            outfile << BasicStats::stat_value_to_string(stats_data[i]);
        }
        outfile << "\n";
    } else {
        std::ofstream outfile(output_file, std::ios::app);
        for (size_t i = 0; i < stats_data.size(); i++) {
            if (i > 0) outfile << ", ";
            outfile << BasicStats::stat_value_to_string(stats_data[i]);
        }
        outfile << "\n";
    }
}

string MSA::calc_unified_stats_header(unordered_map<string, int>& col_names_dict, const vector<vector<string>>& col_names) {
    string header = "code,";
    int count = 1;
    col_names_dict["code"] = 0;
    for (size_t j = 0; j < col_names.size(); j++) {
        for (size_t i = 1; i < col_names[j].size(); ++i) {
            if (!col_names_dict.contains(col_names[j][i])) {
                col_names_dict[col_names[j][i]] = count;
                header += (col_names[j][i] + ",");
                count++;
            }
        }
    }
    return header;
}

void MSA::print_unified_stats_files(const std::filesystem::path& output_dir_path,
                                    const string& batch_name,
                                    bool is_init_file,
                                    unordered_map<string, int>& col_names_dict,
                                    const vector<vector<string>>& col_names,
                                    const vector<vector<StatValue>>& data) {
    
    std::filesystem::path output_file = output_dir_path / ("unified_stats_" + batch_name + ".csv");
    std::ofstream outfile(
        output_file,
        is_init_file ? std::ios::out : std::ios::app
    );

    // Helper lambda to write a CSV row
    auto write_row = [&](const auto& vec, auto formatter) {
        outfile << formatter(vec[0][0]) << ",";
        int count = 1;
        for (size_t j = 0; j < vec.size(); j++) {
            for (size_t i = 1; i < vec[j].size(); ++i) {
                if (col_names_dict.contains(col_names[j][i])) {
                    int index = col_names_dict[col_names[j][i]];
                    if (index == count) {
                        outfile << formatter(vec[j][i]) << ",";
                        count++;
                    }
                }
            }
        }
        outfile << "\n";
    };

    if (is_init_file) {
        string header = calc_unified_stats_header(col_names_dict, col_names);
        outfile << header << endl;
    }

    write_row(data, [](const auto& v) {
        return BasicStats::stat_value_to_string(v);
    });
}

void MSA::print_to_fasta_file(const std::filesystem::path& dir_path) const {
    std::filesystem::path output_file = dir_path / (dataset_name + ".fas");
    std::ofstream outfile(output_file);
    for (size_t i = 0; i < sequences.size(); i++) {
        outfile << ">" << seq_names[i] << "\n";
        outfile << sequences[i] << "\n";
    }
}

