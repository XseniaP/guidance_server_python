#pragma once
#include "node.h"
#include <string>
#include <vector>
#include <set>
#include <memory>
#include <tuple>
#include <filesystem>
using namespace std;

struct PathEntry {
    int start_id;
    int end_id;
    double dist;
};

class UnrootedTree {
public:
    Node* anchor;
    std::vector<std::unique_ptr<Node>> all_nodes;
    std::set<std::string> keys;
    std::string differentiator_key;

    UnrootedTree()
        : anchor(nullptr)
    {
    }

    UnrootedTree(Node* _anchor, const std::vector<std::unique_ptr<Node>>& _all_nodes);
    UnrootedTree(const std::string& newick);

    std::set<std::string> get_internal_edges_set() const;
    int calc_rf(const UnrootedTree& other_tree) const;
    std::vector<double> get_branches_lengths_list() const;

    std::tuple<Node*, std::vector<PathEntry>, double> get_longest_path(Node* u) const;
    std::pair<std::vector<PathEntry>, double> longest_path() const;
    double get_longest_dist_to(Node* dest) const;
    string print_newick() const;
    string get_newick(Node* node) const;

};

// Free functions for Newick parsing
std::string read_newick_from_file(const std::filesystem::path& input_file_path);
Node* create_node_from_children(std::vector<std::vector<Node*>>& open_nodes_per_level,
                                 int level, double branch_length, int node_inx,
                                 std::vector<std::unique_ptr<Node>>& storage);
constexpr double TREE_EPSILON = 0.001;
