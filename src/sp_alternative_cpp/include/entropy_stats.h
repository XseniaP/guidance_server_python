#pragma once
#include "msa_basic_stats.h"
#include <vector>
#include <string>

class EntropyStats : public BasicStats {
public:
    double constant_sites_pct;
    int n_unique_sites;
    double entropy_mean;
    double entropy_25_pct;
    double entropy_75_pct;
    double entropy_sum;
    double entropy_max;

    EntropyStats(const std::string& code, int taxa_num, int msa_length);

    void calc_entropy(const std::vector<std::string>& aln);

    std::vector<StatValue> get_my_features_as_list() const override;
};
