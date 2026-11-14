#include "neighbor_joining.h"
#include <limits>
#include <numeric>
#include <algorithm>
#include <iostream>

using namespace std;

NeighborJoining::NeighborJoining(
    const std::vector<std::vector<double>>& distanceMatrix,
    std::vector<std::unique_ptr<Node>> initial_nodes)
{
    distance_matrix = distanceMatrix;

    //for (int i = 0; i < distanceMatrix.size(); ++i) {
    //    for (int j = 0; j < distanceMatrix[0].size(); ++j) {
    //        //cout << distanceMatrix[i][j] << " ";
    //    }
    //    //cout << endl;
    //}

    // move ownership of nodes
    all_nodes = std::move(initial_nodes);

    // build working_nodes from the owned nodes
    working_nodes.reserve(all_nodes.size());
    for (auto& node : all_nodes) {
        working_nodes.push_back(node.get());
    }

    // initialize q_matrix size
    q_matrix.resize(distance_matrix.size(),
        std::vector<double>(distance_matrix.size(), 0.0));

    tree_res = build_tree();
}

vector<vector<double>> NeighborJoining::calc_q_matrix() const {
    int n = working_nodes.size();
    vector<vector<double>> qm(distance_matrix.size() - 1, vector<double>(n, 0.0));
    
	vector<double> row_sums(distance_matrix.size(), 0.0);
    for (int i = 0; i < distance_matrix.size(); i++) {
        row_sums[i] = accumulate(distance_matrix[i].begin(), distance_matrix[i].end(), 0.0);
    }

    for (int i = 0; i < distance_matrix.size() - 1; i++) {
        for (int j = i + 1; j < static_cast<int>(distance_matrix.size()); j++) {
            qm[i][j] = (n - 2) * distance_matrix[i][j] - row_sums[i] - row_sums[j];
        }
    }
    return qm;
}

NeighborJoining::ClosestPair NeighborJoining::find_closest_pair() const {
    int f_inx = -1, s_inx = -1;
    double min_dist = std::numeric_limits<double>::max();

    for (int i = 0; i < static_cast<int>(q_matrix.size()); i++) {
        for (int j = i + 1; j < static_cast<int>(q_matrix[i].size()); j++) {
            if (q_matrix[i][j] < min_dist) {
                min_dist = q_matrix[i][j];
                f_inx = i;
                s_inx = j;
            }
        }
    }
    return {f_inx, s_inx, min_dist};
}

std::pair<double, double> NeighborJoining::find_delta(int f_inx, int s_inx) const {
    int n = static_cast<int>(working_nodes.size());
    double sum_f = std::accumulate(distance_matrix[f_inx].begin(), distance_matrix[f_inx].end(), 0.0);
    double sum_s = std::accumulate(distance_matrix[s_inx].begin(), distance_matrix[s_inx].end(), 0.0);

    double delta_f = 0.5 * distance_matrix[f_inx][s_inx] +
                     (sum_f - sum_s) / (2.0 * (n - 2));
    double delta_s = distance_matrix[f_inx][s_inx] - delta_f;
    return {delta_f, delta_s};
}

void NeighborJoining::merge_two_clusters() {
    auto [f_inx, s_inx, min_dist] = find_closest_pair();
    auto [delta_f, delta_s] = find_delta(f_inx, s_inx);

    std::vector<std::vector<double>> matrix;
    for (int r = 0; r < static_cast<int>(distance_matrix.size()); r++) {
        std::vector<double> row = distance_matrix[r];
        if (r == f_inx) {
            row[f_inx] = 0;
        } else {
            row[f_inx] = (row[f_inx] + row[s_inx] - distance_matrix[f_inx][s_inx]) / 2.0;
        }
        row.erase(row.begin() + s_inx);

        if (r != s_inx) {
            matrix.push_back(row);
        } else {
            // Merge into f_inx row
            for (int n_idx = 0; n_idx < static_cast<int>(row.size()); n_idx++) {
                if (n_idx != f_inx) {
                    matrix[f_inx][n_idx] = (matrix[f_inx][n_idx] + row[n_idx] -
                                             distance_matrix[f_inx][s_inx]) / 2.0;
                }
            }
        }
    }

    distance_matrix = matrix;

    working_nodes[f_inx]->set_branch_length(delta_f);
    working_nodes[s_inx]->set_branch_length(delta_s);

    // Create new node
    set<string> keysOfNewNode = working_nodes[f_inx]->keys;
    keysOfNewNode.insert(working_nodes[s_inx]->keys.begin(), working_nodes[s_inx]->keys.end());
	double children_bl_sum = delta_f + delta_s;

    auto new_node = std::make_unique<Node>(static_cast<int>(all_nodes.size()),
                                            keysOfNewNode,
                                            std::vector<int>{working_nodes[f_inx]->id, working_nodes[s_inx]->id},
                                            children_bl_sum,
                                            0);

    Node* new_node_ptr = new_node.get();
    all_nodes.push_back(std::move(new_node));
    working_nodes[f_inx]->set_a_father(new_node_ptr->id);
    working_nodes[s_inx]->set_a_father(new_node_ptr->id);

    
    working_nodes[f_inx] = new_node_ptr;
    working_nodes.erase(working_nodes.begin() + s_inx);
}

Node* NeighborJoining::merge_last_three() {
    auto [delta_f, delta_s] = find_delta(0, 1);
    double delta_t = distance_matrix[0][2] - delta_f;
    double deltas[3] = {delta_f, delta_s, delta_t};

    for (int i = 0; i < 3; i++) {
        working_nodes[i]->set_branch_length(deltas[i]);
    }

    double children_bl_sum = 0;
    std::set<std::string>keys;
    for (int i = 0; i < 3; i++) {
        Node* child = working_nodes[i];
        keys.insert(child->keys.begin(), child->keys.end());
        children_bl_sum += child->children_bl_sum + child->branch_length;
    }
    auto new_node = std::make_unique<Node>(static_cast<int>(all_nodes.size()),
        keys,
        std::vector<int>{working_nodes[0]->id, working_nodes[1]->id, working_nodes[2]->id}, children_bl_sum);
    
    all_nodes.push_back(std::move(new_node));

    Node* new_node_ptr = all_nodes.back().get();

    for (int i = 0; i < 3; i++) {
        working_nodes[i]->set_a_father(new_node_ptr->id);
    }
    return new_node_ptr;
}

UnrootedTree NeighborJoining::build_tree() {
    while (working_nodes.size() > 3) {
		//cout << "Merging clusters, remaining count: " << working_nodes.size() << endl;
        q_matrix = calc_q_matrix();
		//cout << "Finding closest pair..." << endl;
        merge_two_clusters();
    }
    q_matrix = calc_q_matrix();
    Node* root = merge_last_three();
    return UnrootedTree(root, std::move(all_nodes));
}

