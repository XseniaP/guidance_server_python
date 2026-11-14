#include "global_alignment.h"
#include "utils.h"
#include <limits>
#include <algorithm>
#include <cmath>
#include <deque>

InternalCell* GlobalAlign::make_cell(int row, int col, const std::vector<InternalCell*>& children, char matrix) {
    auto cell = std::make_unique<InternalCell>(row, col, children, matrix);
    InternalCell* ptr = cell.get();
    cell_storage.push_back(std::move(cell));
    return ptr;
}

GlobalAlign::GlobalAlign(const std::string& seq_a, const std::string& seq_b,
                          const EvoModel& configuration, const std::string& matrix_base_path)
    : seq_a(seq_a), seq_b(seq_b), best_score(0) {
    std::string matrix_file_path = matrix_base_path + "/" + configuration.matrix_file_name + ".txt";
    auto [mm, cdi] = read_matching_matrix(matrix_file_path);
    matching_matrix = mm;
    codes_dict_to_inx = cdi;
    gs = configuration.go_cost;
    ge = configuration.go_cost; // Note: Python code uses go_cost for both gs and ge

    InternalCell* alignment_root = nullptr;
    calculate_alignment_matrix(alignment_root);
    aligned_sequences = get_alignment(alignment_root);
}

int GlobalAlign::compare_codes_score(int row_inx, int col_inx) const {
    if (row_inx == 0 && col_inx == 0) return 0;
    return seq_matrix_matching_score(row_inx, col_inx);
}

int GlobalAlign::seq_matrix_matching_score(int row_inx, int col_inx) const {
    std::string from_key(1, seq_a[row_inx - 1]);
    std::string to_key(1, seq_b[col_inx - 1]);
    int from_inx = codes_dict_to_inx.at(from_key);
    int to_inx = codes_dict_to_inx.at(to_key);
    return matching_matrix[from_inx][to_inx];
}

double GlobalAlign::get_score() const {
    return best_score;
}

void GlobalAlign::calculate_alignment_matrix(InternalCell*& alignment_root) {
    int matrix_row_size = static_cast<int>(seq_a.size()) + 1;
    int matrix_col_size = static_cast<int>(seq_b.size()) + 1;
    double NEG_INF = -std::numeric_limits<double>::max() / 2.0;

    std::vector<std::vector<double>> matrix_m(matrix_row_size, std::vector<double>(matrix_col_size, 0));
    std::vector<std::vector<double>> matrix_x(matrix_row_size, std::vector<double>(matrix_col_size, 0));
    std::vector<std::vector<double>> matrix_y(matrix_row_size, std::vector<double>(matrix_col_size, 0));

    // Only keep current and previous row cells to save memory
    // But we need traceback, so keep last two rows of cells
    // Using a map: row -> col -> {m, x, y} cells
    // Python deletes cells aggressively too, keep prev and current rows
    std::vector<std::vector<InternalCell*>> cells_m(matrix_row_size, std::vector<InternalCell*>(matrix_col_size, nullptr));
    std::vector<std::vector<InternalCell*>> cells_x(matrix_row_size, std::vector<InternalCell*>(matrix_col_size, nullptr));
    std::vector<std::vector<InternalCell*>> cells_y(matrix_row_size, std::vector<InternalCell*>(matrix_col_size, nullptr));

    for (int row_inx = 0; row_inx < matrix_row_size; row_inx++) {
        for (int col_inx = 0; col_inx < matrix_col_size; col_inx++) {
            if (row_inx > 0 && col_inx > 0) {
                // M matrix
                double m_match = matrix_m[row_inx - 1][col_inx - 1];
                double x_match = matrix_x[row_inx - 1][col_inx - 1];
                double y_match = matrix_y[row_inx - 1][col_inx - 1];
                double max_val = std::max({m_match, x_match, y_match});
                matrix_m[row_inx][col_inx] = max_val + compare_codes_score(row_inx, col_inx);

                std::vector<InternalCell*> m_children;
                if (m_match == max_val && cells_m[row_inx - 1][col_inx - 1]) m_children.push_back(cells_m[row_inx - 1][col_inx - 1]);
                if (x_match == max_val && cells_x[row_inx - 1][col_inx - 1]) m_children.push_back(cells_x[row_inx - 1][col_inx - 1]);
                if (y_match == max_val && cells_y[row_inx - 1][col_inx - 1]) m_children.push_back(cells_y[row_inx - 1][col_inx - 1]);
                cells_m[row_inx][col_inx] = make_cell(row_inx, col_inx, m_children, 'm');

                // X matrix
                double m_x_gap = matrix_m[row_inx][col_inx - 1] + gs + ge;
                double x_x_gap = matrix_x[row_inx][col_inx - 1] + ge;
                double y_x_gap = matrix_y[row_inx][col_inx - 1] + gs + ge;
                max_val = std::max({m_x_gap, x_x_gap, y_x_gap});
                matrix_x[row_inx][col_inx] = max_val;

                std::vector<InternalCell*> x_children;
                if (m_x_gap == max_val && cells_m[row_inx][col_inx - 1]) x_children.push_back(cells_m[row_inx][col_inx - 1]);
                if (x_x_gap == max_val && cells_x[row_inx][col_inx - 1]) x_children.push_back(cells_x[row_inx][col_inx - 1]);
                if (y_x_gap == max_val && cells_y[row_inx][col_inx - 1]) x_children.push_back(cells_y[row_inx][col_inx - 1]);
                cells_x[row_inx][col_inx] = make_cell(row_inx, col_inx, x_children, 'x');

                // Y matrix
                double m_y_gap = matrix_m[row_inx - 1][col_inx] + gs + ge;
                double x_y_gap = matrix_x[row_inx - 1][col_inx] + ge + ge;
                double y_y_gap = matrix_y[row_inx - 1][col_inx] + ge;
                max_val = std::max({m_y_gap, x_y_gap, y_y_gap});
                matrix_y[row_inx][col_inx] = max_val;

                std::vector<InternalCell*> y_children;
                if (m_y_gap == max_val && cells_m[row_inx - 1][col_inx]) y_children.push_back(cells_m[row_inx - 1][col_inx]);
                if (x_y_gap == max_val && cells_x[row_inx - 1][col_inx]) y_children.push_back(cells_x[row_inx - 1][col_inx]);
                if (y_y_gap == max_val && cells_y[row_inx - 1][col_inx]) y_children.push_back(cells_y[row_inx - 1][col_inx]);
                cells_y[row_inx][col_inx] = make_cell(row_inx, col_inx, y_children, 'y');

            } else if (col_inx > 0) {
                matrix_m[row_inx][col_inx] = NEG_INF;
                matrix_x[row_inx][col_inx] = col_inx * ge + gs;
                matrix_y[row_inx][col_inx] = NEG_INF;
                cells_x[row_inx][col_inx] = make_cell(row_inx, col_inx,
                    {cells_x[row_inx][col_inx - 1]}, 'x');
            } else if (row_inx > 0) {
                matrix_m[row_inx][col_inx] = NEG_INF;
                matrix_x[row_inx][col_inx] = NEG_INF;
                matrix_y[row_inx][col_inx] = row_inx * ge + gs;
                cells_y[row_inx][col_inx] = make_cell(row_inx, col_inx,
                    {cells_y[row_inx - 1][col_inx]}, 'y');
            } else {
                matrix_m[row_inx][col_inx] = 0;
                matrix_x[row_inx][col_inx] = NEG_INF;
                matrix_y[row_inx][col_inx] = NEG_INF;
                cells_m[row_inx][col_inx] = make_cell(row_inx, col_inx, {}, 'm');
                cells_x[row_inx][col_inx] = make_cell(row_inx, col_inx, {}, 'x');
                cells_y[row_inx][col_inx] = make_cell(row_inx, col_inx, {}, 'y');
            }
        }
    }

    matrix = matrix_m;

    int lr = matrix_row_size - 1;
    int lc = matrix_col_size - 1;
    double m_res = matrix_m[lr][lc];
    double x_res = matrix_x[lr][lc];
    double y_res = matrix_y[lr][lc];
    double max_val = std::max({m_res, x_res, y_res});
    best_score = max_val;

    std::vector<InternalCell*> root_children;
    if (m_res == max_val && cells_m[lr][lc]) root_children.push_back(cells_m[lr][lc]);
    if (x_res == max_val && cells_x[lr][lc]) root_children.push_back(cells_x[lr][lc]);
    if (y_res == max_val && cells_y[lr][lc]) root_children.push_back(cells_y[lr][lc]);
    alignment_root = make_cell(matrix_row_size, matrix_col_size, root_children, 'm');
}

std::vector<AlignmentCandidate> GlobalAlign::get_alignment(InternalCell* root) {
    std::vector<AlignmentCandidate> alignment_results;
    if (!root) return alignment_results;

    struct QueueEntry {
        AlignmentCandidate candidate;
        InternalCell* arrow_from;
        int child_inx;
    };

    std::deque<QueueEntry> alignment_queue;
    for (int ci = 0; ci < static_cast<int>(root->children.size()); ci++) {
        alignment_queue.push_back({AlignmentCandidate(), root, ci});
    }

    while (!alignment_queue.empty()) {
        auto current = alignment_queue.front();
        alignment_queue.pop_front();

        InternalCell* arrow_from = current.arrow_from;
        InternalCell* arrow_to = (current.child_inx < static_cast<int>(arrow_from->children.size()))
                                  ? arrow_from->children[current.child_inx] : nullptr;

        while (arrow_to && arrow_to->row_inx >= 0 && arrow_to->col_inx >= 0) {
            if (arrow_to->row_inx == arrow_from->row_inx - 1 &&
                arrow_to->col_inx == arrow_from->col_inx - 1) {
                if (arrow_to->row_inx < static_cast<int>(seq_a.size())) {
                    current.candidate.add_codes(&seq_a, &seq_b, arrow_to->row_inx, arrow_to->col_inx);
                }
            } else if (arrow_to->row_inx == arrow_from->row_inx - 1) {
                current.candidate.add_codes(&seq_a, nullptr, arrow_to->row_inx, -1);
            } else {
                current.candidate.add_codes(nullptr, &seq_b, -1, arrow_to->col_inx);
            }

            if (arrow_to->children.size() > 1) {
                for (int ci = 1; ci < static_cast<int>(arrow_to->children.size()); ci++) {
                    alignment_queue.push_back({
                        AlignmentCandidate(current.candidate.profile_a, current.candidate.profile_b),
                        arrow_to, ci
                    });
                }
            }

            arrow_from = arrow_to;
            arrow_to = (!arrow_from->children.empty()) ? arrow_from->children[0] : nullptr;
        }
        alignment_results.push_back(current.candidate);
    }
    return alignment_results;
}
