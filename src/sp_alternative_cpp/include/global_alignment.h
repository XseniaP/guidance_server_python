#pragma once
#include "evo_model.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>

struct InternalCell {
    int row_inx;
    int col_inx;
    std::vector<InternalCell*> children;
    char matrix; // 'm', 'x', or 'y'

    InternalCell(int row, int col, const std::vector<InternalCell*>& children, char matrix)
        : row_inx(row), col_inx(col), children(children), matrix(matrix) {}
};

struct AlignmentCandidate {
    std::vector<char> profile_a;
    std::vector<char> profile_b;

    AlignmentCandidate() = default;
    AlignmentCandidate(const std::vector<char>& pa, const std::vector<char>& pb)
        : profile_a(pa), profile_b(pb) {}

    void add_codes(const std::string* seq_a, const std::string* seq_b, int inx_a, int inx_b) {
        char col_a = (seq_a == nullptr) ? '-' : (*seq_a)[inx_a];
        char col_b = (seq_b == nullptr) ? '-' : (*seq_b)[inx_b];
        profile_a.insert(profile_a.begin(), col_a);
        profile_b.insert(profile_b.begin(), col_b);
    }
};

class GlobalAlign {
public:
    std::string seq_a;
    std::string seq_b;
    std::vector<std::vector<double>> matching_matrix;
    std::unordered_map<std::string, int> codes_dict_to_inx;
    double gs;
    double ge;
    std::vector<std::vector<double>> matrix;
    double best_score;
    std::vector<AlignmentCandidate> aligned_sequences;

    // Memory management for internal cells
    std::vector<std::unique_ptr<InternalCell>> cell_storage;

    GlobalAlign(const std::string& seq_a, const std::string& seq_b,
                const EvoModel& configuration, const std::string& matrix_base_path);

    double get_score() const;

private:
    InternalCell* make_cell(int row, int col, const std::vector<InternalCell*>& children, char matrix);
    void calculate_alignment_matrix(InternalCell*& alignment_root);
    int compare_codes_score(int row_inx, int col_inx) const;
    int seq_matrix_matching_score(int row_inx, int col_inx) const;
    std::vector<AlignmentCandidate> get_alignment(InternalCell* root);
};
