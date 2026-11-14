#pragma once
#include "msa_basic_stats.h"
#include "unrooted_tree.h"
#include "enums.h"
#include <vector>
#include <string>

class DistanceLabelsStats : public BasicStats {
public:
    double ssp_from_true;
    double dseq_from_true;
    double dpos_from_true;
    int rf_from_true;

    DistanceLabelsStats(const std::string& code, int taxa_num, int msa_length);

    void set_my_distance_from_true(const std::vector<std::string>& inferred_msa,
                                   const std::vector<std::string>& true_msa);
    void set_rf_from_true(const UnrootedTree& my_tree, const UnrootedTree& true_tree);

    std::pair<std::vector<StatValue>, std::vector<std::string>> get_print_rf() const;

    std::vector<StatValue> get_my_features_as_list() const override;
};
