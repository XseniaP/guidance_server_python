#pragma once
#include "msa_basic_stats.h"
#include <vector>
#include <string>
#include <map>
#include <utility>

class GapStats : public BasicStats {
public:
    int num_gap_segments;
    double num_gap_segments_norm;
    double av_gap_segment_length;
    double gaps_len_one;
    double gaps_len_two;
    double gaps_len_three;
    double gaps_len_four_plus;
    int num_unique_gaps;
    double num_unique_gaps_norm;
    double avg_unique_gap_length;
    int gaps_1seq_len1, gaps_2seq_len1;
    int gaps_1seq_len2, gaps_2seq_len2;
    int gaps_1seq_len3, gaps_2seq_len3;
    int gaps_1seq_len4plus, gaps_2seq_len4plus;
    int gaps_all_except_1_len1;
    int gaps_all_except_1_len2;
    int gaps_all_except_1_len3;
    int gaps_all_except_1_len4plus;
    int num_cols_no_gaps;
    int num_cols_1_gap;
    int num_cols_2_gaps;
    int num_cols_all_gaps_except1;
    int single_char_count;
    int double_char_count;
    int seq_max_len;
    int seq_min_len;

    GapStats(const std::string& code, int taxa_num, int msa_length);

    void calc_gaps_values(const std::vector<std::string>& aln);
    void record_gap_lengths(const std::string& sequence, int seq_index,
                             std::map<std::pair<int,int>, std::vector<int>>& gap_positions,
                             std::map<int, int>& gaps_length_histogram);
    void calculate_counts(const std::map<std::pair<int,int>, std::vector<int>>& gap_positions);

    std::vector<StatValue> get_my_features_as_list() const override;
};
