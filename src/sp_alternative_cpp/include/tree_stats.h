#pragma once
#include "msa_basic_stats.h"
#include "node.h"
#include "unrooted_tree.h"
#include <vector>
#include <string>

class TreeStats : public BasicStats {
public:
    double bl_mean;
    double bl_sum;
    double bl_max;
    double bl_min;
    double bl_25_pct;
    double bl_75_pct;

    double parsimony_mean;
    double parsimony_sum;
    double parsimony_max;
    double parsimony_min;
    double parsimony_25_pct;
    double parsimony_75_pct;

    TreeStats(const std::string& code, int taxa_num, int msa_length);

    void set_tree_stats(const std::vector<double>& bl_list, const UnrootedTree& tree,
                         const std::vector<std::string>& aln, const std::vector<std::string>& names);

    std::vector<StatValue> get_my_features_as_list() const override;
};

std::vector<int> calc_parsimony(const UnrootedTree& unrooted_tree, const std::vector<std::string>& aln,
                                 const std::vector<std::string>& names);
