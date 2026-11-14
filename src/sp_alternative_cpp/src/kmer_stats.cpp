#include "kmer_stats.h"
#include "utils.h"
#include <algorithm>
#include <numeric>

KMerStats::KMerStats(const std::string& code, int taxa_num, int msa_length, int k_value)
    : BasicStats(code, taxa_num, msa_length,
                  {"code", "k_mer_max", "k_mer_average", "k_mer_95_pct", "k_mer_90_pct", "kmer_sum_of_top_10"}),
      k_value(k_value), k_mer_max(0), k_mer_average(0),
      k_mer_95_pct(0), k_mer_90_pct(0), kmer_sum_of_top_10(0) {}

void KMerStats::set_k_mer_features(const std::vector<std::string>& aln) {
    int k = std::min(k_value, static_cast<int>(aln[0].size()) - 1);
    auto histo = calc_kmer_histo(aln, k);
    if (!histo.empty()) {
        k_mer_max = *std::max_element(histo.begin(), histo.end());
        k_mer_average = std::accumulate(histo.begin(), histo.end(), 0.0) / histo.size();
        k_mer_95_pct = static_cast<int>(calc_percentile_int(histo, 95));
        k_mer_90_pct = static_cast<int>(calc_percentile_int(histo, 90));
        std::sort(histo.begin(), histo.end(), std::greater<int>());
        int top = std::min(10, static_cast<int>(histo.size()));
        kmer_sum_of_top_10 = std::accumulate(histo.begin(), histo.begin() + top, 0);
        k_mer_max = histo[0];
    }
}

std::vector<std::string> KMerStats::get_ordered_col_names_with_k_value() const {
    std::string suffix = "_K" + std::to_string(k_value);
    std::vector<std::string> result;
    for (const auto& col_name : ordered_col_names) {
        result.push_back(col_name + suffix);
    }
    return result;
}

std::vector<int> calc_kmer_histo(const std::vector<std::string>& aln, int k) {
    std::unordered_map<std::string, int> histo;
    for (const auto& seq : aln) {
        for (int i = 0; i <= static_cast<int>(seq.size()) - k; i++) {
            std::string kmer = seq.substr(i, k);
            histo[kmer]++;
        }
    }
    std::vector<int> values;
    for (const auto& [key, val] : histo) {
        values.push_back(val);
    }
    return values;
}

std::vector<StatValue> KMerStats::get_my_features_as_list() const {
    return {
        StatValue(code),
        StatValue(k_mer_max), StatValue(k_mer_average),
        StatValue(k_mer_95_pct), StatValue(k_mer_90_pct), StatValue(kmer_sum_of_top_10)
    };
}
