#include <fstream>
#include <algorithm>
#include <queue>
#include <stdexcept>
#include <iostream>
#include "utils.h"
#include "unrooted_tree.h"
using namespace std;


UnrootedTree::UnrootedTree(Node* _anchor, const std::vector<std::unique_ptr<Node>>& _all_nodes) {
    // Deep copy unique_ptrs
    for (auto& n : _all_nodes) {
        all_nodes.push_back(std::make_unique<Node>(*n)); // TODO: this is the problem
    }
	anchor = all_nodes[_anchor->id].get();
    keys = anchor->keys;
    std::vector<std::string> sorted(keys.begin(), keys.end());
    std::sort(sorted.begin(), sorted.end());
    differentiator_key = sorted[0];
}

std::set<std::string> UnrootedTree::get_internal_edges_set() const {
    std::set<std::string> edges;
    for (const auto& n : all_nodes) {
        if (!n->children_ids.empty()) {
            std::string edge_str =
                n->get_keys_unrooted_string(keys, differentiator_key);

            if (!edge_str.empty()) {
                edges.insert(edge_str);
            }
        }
    }
    return edges;
}

int UnrootedTree::calc_rf(const UnrootedTree& other_tree) const {
    auto set_a = get_internal_edges_set();
    auto set_b = other_tree.get_internal_edges_set();
    // Symmetric difference
    int count = 0;
    for (const auto& e : set_a) {
        if (set_b.find(e) == set_b.end()) count++;
    }
    for (const auto& e : set_b) {
        if (set_a.find(e) == set_a.end()) count++;
    }
    return count;
}

std::vector<double> UnrootedTree::get_branches_lengths_list() const {
    std::vector<double> bl_list;
    for (const auto& n : all_nodes) {
        if (n->father_id != -1) {
            bl_list.push_back(std::max(n->branch_length, TREE_EPSILON));
        }
    }
    return bl_list;
}

std::tuple<Node*, std::vector<PathEntry>, double> UnrootedTree::get_longest_path(Node* u) const {
    int nodes_count = static_cast<int>(all_nodes.size());
    std::vector<Node*> nodes_by_id(nodes_count + 1, nullptr);
    std::vector<double> distance(nodes_count + 1, -1.0);
    std::vector<std::vector<PathEntry>> path(nodes_count + 1);

    distance[u->id] = 0;
    nodes_by_id[u->id] = u;

    std::deque<Node*> queue_nodes;
    queue_nodes.push_back(u);

    std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(all_nodes);

    while (!queue_nodes.empty()) {
        Node* front = queue_nodes.front();
        queue_nodes.pop_front();

        for (const auto& adj : front->get_adj(raw_nodes)) {
            if (adj.node_id >= 0 && adj.node_id <= nodes_count && nodes_by_id[adj.node_id] == nullptr) {
                // Find the actual node
                Node* adj_node = nullptr;
                for (const auto& n : all_nodes) {
                    if (n->id == adj.node_id) {
                        adj_node = n.get();
                        break;
                    }
                }
                if (!adj_node) continue;

                nodes_by_id[adj.node_id] = adj_node;
                distance[adj.node_id] = distance[front->id] + adj.dist;
                path[adj.node_id] = path[front->id];
                path[adj.node_id].push_back({front->id, adj.node_id, adj.dist});
                queue_nodes.push_back(adj_node);
            }
        }
    }

    double max_dist = 0;
    int node_index = -1;
    for (int i = 0; i < nodes_count; i++) {
        if (distance[i] > max_dist) {
            max_dist = distance[i];
            node_index = i;
        }
    }

    return {nodes_by_id[node_index], path[node_index], max_dist};
}

std::pair<std::vector<PathEntry>, double> UnrootedTree::longest_path() const {
    auto [node1, path1, dist1] = get_longest_path(all_nodes[0].get());
    auto [node2, path2, dist2] = get_longest_path(node1);
    return {path2, dist2};
}

double UnrootedTree::get_longest_dist_to(Node* dest) const {
    int nodes_count = static_cast<int>(all_nodes.size());
    std::vector<Node*> nodes_by_id(nodes_count + 1, nullptr);
    std::vector<double> distance(nodes_count + 1, -1.0);

    distance[dest->id] = 0;
    nodes_by_id[dest->id] = dest;

    std::deque<Node*> queue_nodes;
    queue_nodes.push_back(dest);

    std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(all_nodes);

    while (!queue_nodes.empty()) {
        Node* front = queue_nodes.front();
        queue_nodes.pop_front();

        for (const auto& adj : front->get_adj(raw_nodes)) {
            if (adj.node_id >= 0 && adj.node_id <= nodes_count && nodes_by_id[adj.node_id] == nullptr) {
                Node* adj_node = nullptr;
                for (const auto& n : all_nodes) {
                    if (n->id == adj.node_id) {
                        adj_node = n.get();
                        break;
                    }
                }
                if (!adj_node) continue;

                nodes_by_id[adj.node_id] = adj_node;
                distance[adj.node_id] = distance[front->id] + adj.dist;
                queue_nodes.push_back(adj_node);
            }
        }
    }

    return *std::max_element(distance.begin(), distance.end());
}

// Free functions

std::string read_newick_from_file(const std::filesystem::path& input_file_path) {
    std::ifstream in_file(input_file_path);
    if (!in_file.is_open()) {
        throw std::runtime_error("Cannot open newick file: " + input_file_path.string());
    }
    std::string line;
    if (std::getline(in_file, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n' || line.back() == ' ')) {
            line.pop_back();
        }
        return line;
    }
    return "";
}

Node* create_node_from_children(std::vector<std::vector<Node*>>& open_nodes_per_level,
    int level, double branch_length, int node_inx,
    std::vector<std::unique_ptr<Node>>& storage)
{
    int child_level = level + 1;
    if (child_level >= static_cast<int>(open_nodes_per_level.size())) {
        return nullptr;
    }

    // Collect keys from children
    std::set<std::string> node_keys;
    std::vector<int> children_ids;
    for (auto* child : open_nodes_per_level[child_level]) {
        node_keys.insert(child->keys.begin(), child->keys.end());
		children_ids.push_back(child->id);
    }

    // Children list
    std::vector<Node*> child_list = open_nodes_per_level[child_level];

    // Create node and push into storage
    auto current_node = std::make_unique<Node>(node_inx, node_keys, children_ids, 0.0, branch_length);
    Node* current_node_ptr = current_node.get();
    storage.push_back(std::move(current_node));

    // Link children to this parent
    for (auto* child : child_list) {
        child->set_a_father(current_node_ptr->id);
    }

    return current_node_ptr;
}

UnrootedTree::UnrootedTree(const std::string& newick)
{
    //std::vector<std::unique_ptr<Node>> all_nodes;
    std::vector<std::vector<Node*>> open_nodes_per_level(100);
    int level = 0;
    std::string current_key;
    std::string branch_length_str;
	anchor = nullptr;

    size_t i = 0;
    while (i < newick.size()) {
        char c = newick[i];
        if (c == ' ') { i++; continue; }

        if (c == '(') {
            level++;
            if (level >= static_cast<int>(open_nodes_per_level.size()))
                open_nodes_per_level.resize(level + 10);
            open_nodes_per_level[level].clear();
            i++;
        }
        else if (c == ':') {
            branch_length_str.clear();
            i++;
            while (i < newick.size() && newick[i] != ')' && newick[i] != ',')
                branch_length_str += newick[i++];
        }
        else if (c == ',' || c == ')') {
            double bl = branch_length_str.empty() ? 0.0 : std::stod(branch_length_str);

            if (!current_key.empty()) {
                // Leaf node
                std::set<std::string> keys = { current_key };
                auto node = std::make_unique<Node>(static_cast<int>(all_nodes.size()), keys, std::vector<int>{}, 0.0, bl);
                Node* node_ptr = node.get();
                all_nodes.push_back(std::move(node));
                open_nodes_per_level[level].push_back(node_ptr);
                current_key.clear();
            }
            else {
                // Internal node
                Node* node_ptr = create_node_from_children(open_nodes_per_level, level, bl, static_cast<int>(all_nodes.size()), all_nodes);
                open_nodes_per_level[level + 1].clear();
                open_nodes_per_level[level].push_back(node_ptr);
            }

            if (c == ')') level--;
            i++;
        }
        else if (c == ';') {
            double bl = branch_length_str.empty() ? 0.0 : std::stod(branch_length_str);
            anchor = create_node_from_children(open_nodes_per_level, 0, bl, static_cast<int>(all_nodes.size()), all_nodes);
            keys = anchor->keys;
            return;
        }
        else {
            current_key.clear();
            while (i < newick.size() && newick[i] != ')' && newick[i] != ',' && newick[i] != ':')
                current_key += newick[i++];
        }
    }
}

string UnrootedTree::print_newick() const {
	string newick_str = get_newick(anchor) + ";";
    //cout << newick_str << std::endl;
	return newick_str;

    //std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(all_nodes);
    //anchor->fill_newick(raw_nodes);
    //out << anchor->newick_part << ";" << std::endl;
}

string UnrootedTree::get_newick(Node* node) const {
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
