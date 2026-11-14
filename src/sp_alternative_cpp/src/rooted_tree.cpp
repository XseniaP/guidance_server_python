#include <algorithm>
#include <cmath>
#include <deque>
#include <numeric>
#include <set>
#include <tuple>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <iostream>

#include "utils.h"
#include "rooted_tree.h"

using namespace std;


RootedTree::RootedTree(Node* r,
    std::vector<std::unique_ptr<Node>> nodes,
    const std::set<std::string>& k)
{
    root = r;
    all_nodes = std::move(nodes);
    keys = k;
}


RootedTree::RootedTree(const UnrootedTree& unrooted, RootingPoint rp) {

    all_nodes.reserve(unrooted.all_nodes.size());
    keys = unrooted.keys;

    // copy nodes
    for (const auto& n : unrooted.all_nodes) {
        all_nodes.push_back(std::make_unique<Node>(*n));
    }

    //string newick_str = tree.print_newick();
    Node* start = all_nodes[rp.start_id].get();
    Node* end = all_nodes[rp.end_id].get();

    // Create new root node
    std::set<std::string> root_keys = unrooted.anchor->keys;

    auto new_root_ptr = std::make_unique<Node>(
        static_cast<int>(all_nodes.size()),
        root_keys, std::vector<int>{start->id, end->id}, 0.0, 0.0);
    
    all_nodes.push_back(std::move(new_root_ptr));
    root = all_nodes.back().get();

    Node* lower_node = start;
    Node* higher_node = end;
    double lower_node_bl = rp.dist_from_start;
    double higher_node_bl = rp.dist_from_end;

    if (find(start->children_ids.begin(), start->children_ids.end(), end->id) != start->children_ids.end()) {
        lower_node = end;
        higher_node = start;
        lower_node_bl = rp.dist_from_end;
        higher_node_bl = rp.dist_from_start;
    }

    lower_node->father_id = root->id;
    lower_node->branch_length = lower_node_bl;
    higher_node->children_ids.erase(remove(higher_node->children_ids.begin(), higher_node->children_ids.end(), lower_node->id), higher_node->children_ids.end());
    int prev_father_id = higher_node->father_id;
    double prev_dist_to_father = higher_node->branch_length;
    if (higher_node->father_id != -1) {
        higher_node->children_ids.push_back(higher_node->father_id);
    }
	higher_node->father_id = root->id;
    higher_node->branch_length = higher_node_bl;

    while (prev_father_id != -1) {
        int working_node_id = prev_father_id;
		Node* working_node = all_nodes[working_node_id].get();
        double working_node_to_father_bl = working_node->branch_length;
        prev_father_id = working_node->father_id;
        working_node->father_id = higher_node->id;
        working_node->branch_length = prev_dist_to_father;
        working_node->children_ids.erase(remove(working_node->children_ids.begin(), working_node->children_ids.end(), higher_node->id), working_node->children_ids.end());
        if (prev_father_id != -1) {
            working_node->children_ids.push_back(prev_father_id);
        }
        prev_dist_to_father = working_node_to_father_bl;
        higher_node = working_node;
    }
}

void RootedTree::calc_clustal_w(std::vector<Node*> raw_nodes) {
    std::deque<Node*> nodes_to_recalc;
    nodes_to_recalc.push_back(root);
    while (!nodes_to_recalc.empty()) {
        Node* node = nodes_to_recalc.front();
        nodes_to_recalc.pop_front();
        fill_nodes_w(node, nodes_to_recalc, raw_nodes);
    }
}

void RootedTree::calc_seq_w() {
	std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(all_nodes);
    calc_clustal_w(raw_nodes);
    for (const auto& node : all_nodes) {
        if (node->children_ids.empty()) {
            seq_weight_dict[*node->keys.begin()] = node->weight;
        }
    }
}

RootingPoint calc_mid_point(const UnrootedTree& unrooted) {
    auto [path, max_dist] = unrooted.longest_path();
    double half_length = 0;
    for (const auto& b : path) {
        if (half_length + b.dist > max_dist / 2.0) {
            double dist_from_start = max_dist / 2.0 - half_length;
            return {b.start_id, b.end_id, dist_from_start, b.dist - dist_from_start, 0};
        }
        half_length += b.dist;
    }
    // Fallback: use last edge
    if (!path.empty()) {
        return {path.back().start_id, path.back().end_id, path.back().dist / 2.0, path.back().dist / 2.0, 0};
    }
    return {0, 0, 0, 0, 0};
}

double sum_bl_up_to_node_id(Node* origin, int dest_id, const std::vector<Node*>& all_nodes) {
    double w_to_orig = 0;
    std::set<int> visited_node_ids = {dest_id, origin->id};

    struct QueueEntry { Node* node; double dist; };
    std::vector<QueueEntry> queue;
    for (const auto& adj : origin->get_adj(all_nodes)) {
        if (visited_node_ids.find(adj.node_id) == visited_node_ids.end()) {
            // Find actual node
            Node* adj_node = nullptr;
            for (auto* n : all_nodes) {
                if (n->id == adj.node_id) { adj_node = n; break; }
            }
            if (adj_node) queue.push_back({adj_node, adj.dist});
        }
    }

    while (!queue.empty()) {
        auto next = queue.back();
        queue.pop_back();
        if (visited_node_ids.find(next.node->id) == visited_node_ids.end()) {
            w_to_orig += next.dist;
            visited_node_ids.insert(next.node->id);
            for (const auto& adj : next.node->get_adj(all_nodes)) {
                if (visited_node_ids.find(adj.node_id) == visited_node_ids.end()) {
                    Node* adj_node = nullptr;
                    for (auto* n : all_nodes) {
                        if (n->id == adj.node_id) { adj_node = n; break; }
                    }
                    if (adj_node) queue.push_back({adj_node, adj.dist});
                }
            }
        }
    }
    return w_to_orig;
}

RootingPoint calc_potential_root_on_branch(double bl, double w_to_orig, double w_to_dest,
                                            double delta, int origin, int dest, double min_bl) {
    min_bl /= 10.0;
    double dist_from_start;
    if (delta > 0) {
        dist_from_start = bl / 2.0 - (w_to_orig - w_to_dest) / 2.0;
    } else if (w_to_orig > w_to_dest) {
        dist_from_start = min_bl;
    } else {
        dist_from_start = bl - min_bl;
    }
    return {origin, dest, dist_from_start, bl - dist_from_start, 0};
}

std::vector<RootingPoint> calc_min_differential_sum(const UnrootedTree& unrooted) {
    
    std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(unrooted.all_nodes);
    
    struct BranchInfo {
        int origin, dest;
        double bl, w_to_orig, w_to_dest, delta;
    };

    std::vector<double> all_bl = unrooted.get_branches_lengths_list();
    double total_branch_length = std::accumulate(all_bl.begin(), all_bl.end(), 0.0);
    double min_bl = *std::min_element(all_bl.begin(), all_bl.end());

    std::map<std::string, BranchInfo> all_branches;



    for (const auto& node : unrooted.all_nodes) {
        for (const auto& adj : node->get_adj(raw_nodes)) {
            int min_id = std::min(adj.node_id, node->id);
            int max_id = std::max(adj.node_id, node->id);
            std::string key = "o:" + std::to_string(min_id) + "-d:" + std::to_string(max_id);

            if (all_branches.find(key) == all_branches.end()) {
                double w_to_orig_val = sum_bl_up_to_node_id(raw_nodes[min_id], max_id, raw_nodes);
                double w_to_dest_val = total_branch_length - w_to_orig_val - adj.dist;
                double delta_val = adj.dist - std::abs(w_to_dest_val - w_to_orig_val);

                all_branches[key] = {min_id, max_id, adj.dist, w_to_orig_val, w_to_dest_val, delta_val};
            }
        }
    }

    std::vector<RootingPoint> res;
    for (const auto& [key, b] : all_branches) {
        if (b.delta > 0) {
            res.push_back(calc_potential_root_on_branch(b.bl, b.w_to_orig, b.w_to_dest,
                                                         b.delta, b.origin, b.dest, min_bl));
        }
    }

    if (res.empty()) {
        // Find branch with max delta
        auto best = std::max_element(all_branches.begin(), all_branches.end(),
                                      [](const auto& a, const auto& b) { return a.second.delta < b.second.delta; });
        res.push_back(calc_potential_root_on_branch(best->second.bl, best->second.w_to_orig,
                                                     best->second.w_to_dest, best->second.delta,
                                                     best->second.origin, best->second.dest, min_bl));
    }

    return res;
}

void fill_nodes_w(Node* node, std::deque<Node*>& nodes_to_recalc, std::vector<Node*>& all_nodes) {
    std::vector<double> w_list;
    if (node->father_id != -1) {
		Node* father = all_nodes[node->father_id];
        w_list = father->w_from_root;
    }
    w_list.push_back(node->branch_length);
    node->set_w_from_root(w_list);

    if (!node->children_ids.empty()) {
        nodes_to_recalc.push_back(all_nodes[node->children_ids[0]]);
        if (node->children_ids.size() > 1) {
            nodes_to_recalc.push_back(all_nodes[node->children_ids[1]]);
        }
    } else {
        node->set_weight_from_root();
    }
}

double longest_dist_from_virtual_root(
    const UnrootedTree& tree,
    Node* start,
    Node* end,
    double d_start,
    double d_end)
{
    std::deque<std::pair<Node*, double>> q;
    std::unordered_set<int> visited;

    double max_dist = 0.0;

    q.push_back({ start, d_start });
    q.push_back({ end,   d_end });

    visited.insert(start->id);
    visited.insert(end->id);

    std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(tree.all_nodes);

    while (!q.empty())
    {
        auto [node, dist] = q.front();
        q.pop_front();

        max_dist = std::max(max_dist, dist);

        for (const auto& adj : node->get_adj(raw_nodes))
        {
            if (visited.count(adj.node_id))
                continue;

            Node* next = tree.all_nodes[adj.node_id].get();

            visited.insert(adj.node_id);
            q.push_back({ next, dist + adj.dist });
        }
    }

    return max_dist;
}

string RootedTree::print_newick() const {
    string newick_str = get_newick(root) + ";";
    // cout << newick_str << std::endl;
    return newick_str;


    //std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(all_nodes);
    //anchor->fill_newick(raw_nodes);
    //out << anchor->newick_part << ";" << std::endl;
}


string RootedTree::get_newick(Node* node) const {
    string newick_str;
    if (node->children_ids.size() == 0) {
        newick_str = *node->keys.begin();
        newick_str += ":";
    }
    else {
        newick_str = "(";
        vector<string> strings;
        for (int i = 0; i < node->children_ids.size(); ++i) {
            int child_id = node->children_ids[i];
            Node* c = all_nodes[child_id].get();
            strings.push_back(get_newick(c));
        }
        sort(strings.begin(), strings.end());
        for (int i = 0; i < strings.size(); ++i) {
            newick_str += strings[i];
            if (i < strings.size() - 1) {
                newick_str += ",";
            }
            else {
                newick_str += "):";
            }
        }
    }
    newick_str += to_string(node->branch_length);
    return newick_str;
}

 set<string> RootedTree::get_my_keys_set(Node* node) const {
    set<string> keys = {};
    if (node->children_ids.size() == 0) {
        keys = node->keys;
    }
    else {
        for (int i = 0; i < node->children_ids.size(); ++i) {
            int child_id = node->children_ids[i];
            Node* c = all_nodes[child_id].get();
            set<string> set2 = get_my_keys_set(c);
            std::set_union(
                keys.begin(), keys.end(),
                set2.begin(), set2.end(),
                std::inserter(keys, keys.begin())
            );
        }
    }
    return keys;
 }

 int RootedTree::get_my_children_count(Node* node) const {
     int children_count = 0;
     if (node->children_ids.size() > 0) {
         for (int i = 0; i < node->children_ids.size(); ++i) {
             int child_id = node->children_ids[i];
             Node* c = all_nodes[child_id].get();
             children_count += (get_my_children_count(c) + 1);
         }
     }
     return children_count;
 }

RootingPoint get_rooting_point(RootingMethods rooting_method, const UnrootedTree& unrooted) {
    // Find midpoint rooting location
    if (rooting_method == RootingMethods::LONGEST_PATH_MID) {
        return calc_mid_point(unrooted);
    }
    auto rooting_points = calc_min_differential_sum(unrooted);
    return find_shallowest_tree(unrooted, rooting_points);
    
}

RootingPoint find_shallowest_tree(
    const UnrootedTree& ut,
    const std::vector<RootingPoint>& rooting_points)
{
    RootingPoint best = rooting_points[0];
    double best_depth = std::numeric_limits<double>::max();

    for (const auto& rp : rooting_points)
    {
        Node* start = ut.all_nodes[rp.start_id].get();
        Node* end = ut.all_nodes[rp.end_id].get();

        double depth = longest_dist_from_virtual_root(
            ut,
            start,
            end,
            rp.dist_from_start,
            rp.dist_from_end);

        if (depth < best_depth)
        {
            best_depth = depth;
            best = rp;
            best.longest_dist = depth;
        }
    }

    return best;
}
