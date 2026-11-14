#include "entropy_stats.h"
#include "utils.h"
#include <cmath>
#include <map>
#include <set>
#include <algorithm>
#include <numeric>

EntropyStats::EntropyStats(const std::string& code, int taxa_num, int msa_length)
    : BasicStats(code, taxa_num, msa_length,
                  {"code", "constant_sites_pct", "n_unique_sites", "entropy_mean",
                   "entropy_25_pct", "entropy_75_pct", "entropy_sum", "entropy_max"}),
      constant_sites_pct(0), n_unique_sites(0), entropy_mean(0),
      entropy_25_pct(0), entropy_75_pct(0), entropy_sum(0), entropy_max(0) {}

void EntropyStats::calc_entropy(const std::vector<std::string>& aln) {
    int seq_count = static_cast<int>(aln.size());
    int seq_len = static_cast<int>(aln[0].size());

    std::vector<double> entropy_per_col;
    std::set<std::string> unique_columns;

    for (int col = 0; col < seq_len; col++) {
        // Count non-gap characters
        std::map<char, int> counts;
        int total = 0;
        std::string col_str;
        for (int i = 0; i < seq_count; i++) {
            char ch = aln[i][col];
            col_str += ch;
            if (ch != '-') {
                counts[ch]++;
                total++;
            }
        }
        unique_columns.insert(col_str);

        // Calculate entropy
        double col_entropy = 0;
        if (total > 0) {
            for (const auto& [ch, cnt] : counts) {
                double p = static_cast<double>(cnt) / total;
                col_entropy += -p * std::log(p);
            }
        }
        entropy_per_col.push_back(col_entropy);
    }

    // Constant sites: columns with 0 entropy
    int constant_count = 0;
    for (double e : entropy_per_col) {
        if (e == 0.0) constant_count++;
    }
    constant_sites_pct = static_cast<double>(constant_count) / entropy_per_col.size();
    n_unique_sites = static_cast<int>(unique_columns.size());

    entropy_max = *std::max_element(entropy_per_col.begin(), entropy_per_col.end());
    entropy_sum = std::accumulate(entropy_per_col.begin(), entropy_per_col.end(), 0.0);
    entropy_mean = entropy_sum / entropy_per_col.size();
    entropy_25_pct = calc_percentile(entropy_per_col, 25);
    entropy_75_pct = calc_percentile(entropy_per_col, 75);
}

std::vector<StatValue> EntropyStats::get_my_features_as_list() const {
    return {
        StatValue(code),
        StatValue(constant_sites_pct), StatValue(n_unique_sites), StatValue(entropy_mean),
        StatValue(entropy_25_pct), StatValue(entropy_75_pct), StatValue(entropy_sum), StatValue(entropy_max)
    };
}
