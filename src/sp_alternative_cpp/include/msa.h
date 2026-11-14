#pragma once
#include "config.h"
#include "sp_score.h"
#include "unrooted_tree.h"
#include "msa_basic_stats.h"
#include <vector>
#include <string>
#include <filesystem>
#include <memory>

class MSA {
public:
    std::string dataset_name;
    std::vector<std::string> sequences;
    std::vector<std::string> seq_names;
    std::unique_ptr<UnrootedTree> tree;

    explicit MSA(const std::string& dataset_name);

    void add_sequence(const std::string& sequence, const std::string& seq_name);
    void set_sequences(const std::vector<std::string>& seqs, const std::vector<std::string>& names);
    void read_from_fasta(const std::filesystem::path& file_path);
    void order_sequences(const std::vector<std::string>& ordered_seq_names);
    void build_nj_tree();
    void set_tree(UnrootedTree&& t);

    int get_msa_len() const;
    int get_taxa_num() const;

    void calc_and_print_stats(const MSA& true_msa, const Configuration& config,
                               const std::vector<SPScore>& sp_models,
                               const std::filesystem::path& output_dir_path,
                               const std::string& batch_name,
                               const UnrootedTree* true_tree,
                               bool is_init_file,
                               unordered_map<string, int>& col_names_dict);

    static void print_stats_files(const std::vector<StatValue>& stats_data,
                                  const string& batch_name,
                                  const std::filesystem::path& output_dir_path,
                                  const std::string& file_name,
                                  bool is_init_file,
                                  const std::vector<std::string>& col_names,
                                  const std::string& model_name = "",
                                  double go_val = 0, double ge_val = 0,
                                  const std::string& k_value = "");

    static void print_unified_stats_files(const std::filesystem::path& output_dir_path,
                                   const std::string& batch_name,
                                bool is_init_file,
                                unordered_map<string, int>& col_names_dict,
                                const vector<vector<string>>& col_names,
                                const vector<vector<StatValue>>& data);

    static string calc_unified_stats_header(unordered_map<string, int>& col_names_dict,
                                            const vector<vector<string>>& col_names);

    void print_to_fasta_file(const std::filesystem::path& dir_path) const;
};