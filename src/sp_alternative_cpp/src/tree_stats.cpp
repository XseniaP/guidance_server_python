#include "tree_stats.h"
#include "utils.h"
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <set>
#include <functional>
#include "utils.h"

TreeStats::TreeStats(const std::string& code, int taxa_num, int msa_length)
    : BasicStats(code, taxa_num, msa_length,
                  {"code", "bl_sum", "bl_mean", "bl_25_pct", "bl_75_pct", "bl_max", "bl_min",
                   "parsimony_mean", "parsimony_sum", "parsimony_max", "parsimony_min",
                   "parsimony_25_pct", "parsimony_75_pct"}),
      bl_mean(-1), bl_sum(-1), bl_max(-1), bl_min(-1), bl_25_pct(-1), bl_75_pct(-1),
      parsimony_mean(-1), parsimony_sum(-1), parsimony_max(-1), parsimony_min(-1),
      parsimony_25_pct(-1), parsimony_75_pct(-1) {}

void TreeStats::set_tree_stats(const std::vector<double>& bl_list, const UnrootedTree& tree,
                                 const std::vector<std::string>& aln, const std::vector<std::string>& names) {
    std::vector<double> bl = bl_list;
    bl_sum = std::accumulate(bl.begin(), bl.end(), 0.0);
    bl_mean = bl_sum / bl.size();
    bl_max = *std::max_element(bl.begin(), bl.end());
    bl_min = *std::min_element(bl.begin(), bl.end());
    bl_25_pct = calc_percentile(bl, 25);
    bl_75_pct = calc_percentile(bl, 75);

    auto parsimony_list = calc_parsimony(tree, aln, names);
    std::vector<double> pars_double(parsimony_list.begin(), parsimony_list.end());
    parsimony_sum = std::accumulate(pars_double.begin(), pars_double.end(), 0.0);
    parsimony_mean = parsimony_sum / pars_double.size();
    parsimony_max = *std::max_element(pars_double.begin(), pars_double.end());
    parsimony_min = *std::min_element(pars_double.begin(), pars_double.end());
    parsimony_25_pct = calc_percentile(pars_double, 25);
    parsimony_75_pct = calc_percentile(pars_double, 75);
}


std::vector<StatValue> TreeStats::get_my_features_as_list() const {
    return {
        StatValue(code),
        StatValue(bl_sum), StatValue(bl_mean), StatValue(bl_25_pct), StatValue(bl_75_pct),
        StatValue(bl_max), StatValue(bl_min),
        StatValue(parsimony_mean), StatValue(parsimony_sum), StatValue(parsimony_max), StatValue(parsimony_min),
        StatValue(parsimony_25_pct), StatValue(parsimony_75_pct)
    };
}

std::set<char> intersect_sets(const std::set<char>& a, const std::set<char>& b) {
    std::set<char> result;

    std::set_intersection(
        a.begin(), a.end(),
        b.begin(), b.end(),
        std::inserter(result, result.begin())
    );

    return result;
}


std::vector<int> calc_parsimony(const UnrootedTree& unrooted_tree,
    const std::vector<std::string>& aln,
    const std::vector<std::string>& names) {

	// --- 1. Get raw pointers from unique_ptrs ---
	std::vector<Node*> raw_nodes = get_raw_pointers_from_unique(unrooted_tree.all_nodes);

    // --- 2. Prepare postorder traversal ---
    std::vector<Node*> postorder;
    std::function<void(Node*)> dfs = [&](Node* node) {
        for (int node_id : node->children_ids) {
			Node* c = unrooted_tree.all_nodes[node_id].get();
            if (c->children_ids.size() > 0) {
                dfs(c);
            }
            postorder.push_back(c);
        }   
    };
    dfs(unrooted_tree.anchor);
    postorder.push_back(unrooted_tree.anchor);

    // --- 3. Map sequence names to alignment indices ---
    std::unordered_map<std::string, int> seq_name_to_index;
    for (int i = 0; i < static_cast<int>(names.size()); i++) {
        seq_name_to_index[names[i]] = i;
    }

    // --- 4. Fitch parsimony per column ---
    std::vector<int> parsimony_per_col;
    int n_cols = static_cast<int>(aln[0].size());
    for (int col = 0; col < n_cols; col++) {
        std::unordered_map<Node*, std::set<char>> states;
        int score = 0;

        for (Node* node : postorder) {
            if (node->children_ids.empty()) {
                // leaf

                char c = aln[seq_name_to_index[*node->keys.begin()]][col];
                states[node] = { c };
            }
            else if (node->children_ids.size() == 2) {
                std::set<char> intersection, uni;
                bool first = true;
                for (int child_id : node->children_ids) {
                    if (child_id == node->father_id) continue;
					Node* child = unrooted_tree.all_nodes[child_id].get();
                    if (first) {
                        intersection = states[child];
                        uni = states[child];
                        first = false;
                    }
                    else {
                        std::set<char> tmp;
                        std::set_intersection(intersection.begin(), intersection.end(),
                            states[child].begin(), states[child].end(),
                            std::inserter(tmp, tmp.begin()));
                        intersection = tmp;
                        uni.insert(states[child].begin(), states[child].end());
                    }
                }
                if (!intersection.empty()) states[node] = intersection;
                else {
                    states[node] = uni;
                    score++;
                }
            }
            else { // this is the 3-way root:
                std::set<char> intersection, uni;
                bool first = true;
                Node* first_child = raw_nodes[unrooted_tree.anchor->children_ids[0]];
                Node* second_child = raw_nodes[unrooted_tree.anchor->children_ids[1]];
                Node* third_child = raw_nodes[unrooted_tree.anchor->children_ids[2]];

                std::set<char> state_a = states[first_child];
                std::set<char> state_b = states[second_child];
                std::set<char> state_c = states[third_child];
                std::set<char> state_a_b;

				std::set<char> intersection_a = intersect_sets(state_a, state_b);
                if (intersection_a.size() == 0) {
                    uni.insert(state_a.begin(), state_a.end());
                    uni.insert(state_b.begin(), state_b.end());
                    state_a_b = uni;
                    score++;
                }
                else {
                    state_a_b = intersection_a;
                }
                std::set<char> intersection_b = intersect_sets(state_a_b, state_c);
                if (intersection_b.size() == 0) {
                    uni.insert(state_a_b.begin(), state_a_b.end());
                    uni.insert(state_c.begin(), state_c.end());
                    score++;
                }
            }
        }

        parsimony_per_col.push_back(score);
    }

    return parsimony_per_col;
}
