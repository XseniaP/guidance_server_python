#pragma once
#include "node.h"
#include "unrooted_tree.h"
#include <vector>
#include <memory>
#include <optional>


class NeighborJoining {
public:
    std::vector<std::vector<double>> distance_matrix;
    std::vector<Node*> working_nodes;
    std::vector<std::vector<double>> q_matrix;
    std::vector<std::unique_ptr<Node>> all_nodes;
    std::optional<UnrootedTree> tree_res;

    NeighborJoining(
                const std::vector<std::vector<double>>& distanceMatrix,
                std::vector<std::unique_ptr<Node>> initial_nodes
            );

    std::vector<std::vector<double>> calc_q_matrix() const;

    struct ClosestPair {
        int f_inx;
        int s_inx;
        double min_dist;
    };
    ClosestPair find_closest_pair() const;

    void merge_two_clusters();
    std::pair<double, double> find_delta(int f_inx, int s_inx) const;
    Node* merge_last_three();
    UnrootedTree build_tree();
};