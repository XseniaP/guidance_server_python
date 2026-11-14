#pragma once
#include "msa_basic_stats.h"
#include <vector>
#include <string>
#include <unordered_map>

class KMerStats : public BasicStats {
public:
    int k_value;
    int k_mer_max;
    double k_mer_average;
    int k_mer_95_pct;
    int k_mer_90_pct;
    int kmer_sum_of_top_10;

    KMerStats(const std::string& code, int taxa_num, int msa_length, int k_value);

    void set_k_mer_features(const std::vector<std::string>& aln);

    std::vector<std::string> get_ordered_col_names_with_k_value() const;

    std::vector<StatValue> get_my_features_as_list() const override;
};

std::vector<int> calc_kmer_histo(const std::vector<std::string>& aln, int k);
