#pragma once
#include <set>
#include <string>
#include <vector>
#include <algorithm>

struct AdjEntry {
    int node_id;
    double dist;
};

class Node {
public:
    int id;
    std::set<std::string> keys;
    int father_id;
    double branch_length;
    std::vector<int> children_ids;
    std::string newick_part;
    std::set<std::string> parsimony_set;
    double children_bl_sum;
    double bl_sum_on_differentiator;
    double bl_sum_on_non_differentiator;
    int rank_from_root;
    std::vector<double> w_from_root;
    double weight;

    Node(int node_id, const std::set<std::string>& keys, const std::vector<int>& children_ids,
         double children_bl_sum, double branch_length = 0.0);

    //static Node* create_from_children(const std::vector<Node*>& children_list, int inx,
    //                                  std::vector<Node*>& storage);
    //// Clone the node without children/father pointers
   

    void set_a_father(int other_node_id);
    void set_branch_length(double bl) { branch_length = bl; };
    std::string get_keys_rooted_string() const;
    std::string get_keys_unrooted_string(const std::set<std::string>& tree_keys,
                                          const std::string& differentiator_key) const;
    void fill_newick(std::vector<Node*>all_nodes);
    void set_parsimony_set(const std::set<std::string>& new_set);
    std::vector<AdjEntry> get_adj(std::vector<Node*>all_nodes) const;
    void update_children_only(const std::vector<int>& children_list);
    void set_rank_from_root(int rank);
    void set_w_from_root(const std::vector<double>& w_list);
    void update_data_from_children(std::vector<Node*>all_nodes);
    void set_weight_from_root();
};
