#pragma once
#include "node.h"
#include "unrooted_tree.h"
#include "enums.h"
#include <vector>
#include <unordered_map>
#include <string>
#include <deque>

struct RootingPoint {
    int start_id;
    int end_id;
    double dist_from_start;
    double dist_from_end;
    double longest_dist; // used for shallowest tree comparison
};


class RootedTree {
public:
    Node* root;
    std::vector<unique_ptr<Node>> all_nodes;
    std::set<std::string> keys;
    std::unordered_map<std::string, double> seq_weight_dict;

    //RootedTree(Node* root, const std::vector<Node*>& all_nodes, const std::set<std::string>& keys);
    RootedTree(Node* r, std::vector<std::unique_ptr<Node>> nodes, const std::set<std::string>& k);

    RootedTree(const UnrootedTree& unrooted, RootingPoint rp);

    void calc_clustal_w(std::vector<Node*> raw_nodes);
    void calc_seq_w();

    string get_newick(Node* node) const;
    string print_newick() const;
    set<string> get_my_keys_set(Node* node) const;
    int get_my_children_count(Node* node) const;
};



// Helper free functions
RootingPoint calc_mid_point(const UnrootedTree& unrooted);
std::vector<RootingPoint> calc_min_differential_sum(const UnrootedTree& unrooted);
void recalc_tree_down(Node* node, Node* father, int broke_id,
                       std::vector<std::tuple<Node*, Node*, int>>& nodes_to_recalc,
                       std::vector<Node*>& all_nodes);
void fill_nodes_w(Node* node, std::deque<Node*>& nodes_to_recalc, std::vector<Node*>& all_nodes);
double sum_bl_up_to_node_id(Node* origin, int dest_id, const std::vector<Node*>& all_nodes);
RootingPoint calc_potential_root_on_branch(double bl, double w_to_orig, double w_to_dest,
                                            double delta, int origin, int dest, double min_bl);

RootingPoint find_shallowest_tree(
    const UnrootedTree& ut,
    const std::vector<RootingPoint>& rooting_points);

RootingPoint get_rooting_point(RootingMethods rooting_method, const UnrootedTree& unrooted);


double longest_dist_from_virtual_root(
    const UnrootedTree& tree,
    Node* start,
    Node* end,
    double d_start,
    double d_end);