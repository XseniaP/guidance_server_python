#pragma once
#include <vector>
#include <string>

#include "msa_basic_stats.h"


using namespace std;

class SubsMatrixCounterStats : public BasicStats {
public:
    std::vector<int> values;
    std::vector<string> col_names;


    static std::vector<std::string> build_col_names(
        const std::vector<std::string>& labels
    ) {
        std::vector<std::string> cols;
        size_t n = labels.size();
        cols.reserve(n * n / 2);
        cols.push_back("code_subs_matrix");

        for (size_t i = 0; i < n - 1; ++i) {
            for (size_t j = i + 1; j < n - 1; ++j) {
                cols.push_back(labels[i] + "<->" + labels[j]);
            }
        }
        return cols;
    }

    SubsMatrixCounterStats(
        const std::string& code,
        int taxa_num,
        int msa_length,
        const std::vector<std::vector<int>>& matrix,
        const std::vector<std::string>& labels
    )
        : BasicStats(code, taxa_num, msa_length, build_col_names(labels))
    {
        size_t n = labels.size();

        values.reserve(n * n / 2);

        for (size_t i = 0; i < n - 1; ++i) {
            for (size_t j = i + 1; j < n - 1; ++j) {
                values.push_back(matrix[i][j]);
            }
        }
    }

    vector<string>get_ordered_col_names_with_labels_value(const std::vector<std::string>& labels) const {
        return build_col_names(labels);
    }

    vector<StatValue> get_my_features_as_list() const override {
        vector<StatValue> res;
        res.reserve(values.size() + 1);

        res.push_back(StatValue(code));
        res.insert(res.end(), values.begin(), values.end());
        return res;
    }
};