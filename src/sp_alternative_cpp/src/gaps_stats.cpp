#include "gaps_stats.h"
#include <algorithm>

GapStats::GapStats(const std::string& code, int taxa_num, int msa_length)
    : BasicStats(code, taxa_num, msa_length,
                  {"code",
                   "num_gap_segments_norm", "av_gap_segment_length", "gaps_len_one", "gaps_len_two",
                   "gaps_len_three", "gaps_len_four_plus", "num_unique_gaps", "num_unique_gaps_norm", "avg_unique_gap_length",
                   "gaps_1seq_len1", "gaps_2seq_len1", "gaps_1seq_len2", "gaps_2seq_len2",
                   "gaps_1seq_len3", "gaps_2seq_len3", "gaps_1seq_len4plus", "gaps_2seq_len4plus",
                   "gaps_all_except_1_len1", "gaps_all_except_1_len2", "gaps_all_except_1_len3", "gaps_all_except_1_len4plus",
                   "num_cols_no_gaps", "num_cols_1_gap", "num_cols_2_gaps", "num_cols_all_gaps_except1",
                   "single_char_count", "double_char_count",
                   "seq_max_len", "seq_min_len"}),
      num_gap_segments(-1), num_gap_segments_norm(0), av_gap_segment_length(0),
      gaps_len_one(0), gaps_len_two(0), gaps_len_three(0), gaps_len_four_plus(0),
      num_unique_gaps(0), num_unique_gaps_norm(0), avg_unique_gap_length(0),
      gaps_1seq_len1(0), gaps_2seq_len1(0), gaps_1seq_len2(0), gaps_2seq_len2(0),
      gaps_1seq_len3(0), gaps_2seq_len3(0), gaps_1seq_len4plus(0), gaps_2seq_len4plus(0),
      gaps_all_except_1_len1(0), gaps_all_except_1_len2(0),
      gaps_all_except_1_len3(0), gaps_all_except_1_len4plus(0),
      num_cols_no_gaps(0), num_cols_1_gap(0), num_cols_2_gaps(0), num_cols_all_gaps_except1(0),
      single_char_count(0), double_char_count(0), seq_max_len(-1), seq_min_len(-1) {}

void GapStats::calc_gaps_values(const std::vector<std::string>& aln) {
    int min_length = 1000000000;
    int max_length = -1;
    std::map<std::pair<int,int>, std::vector<int>> gap_positions;
    std::map<int, int> gaps_length_histogram;

    for (int seq_index = 0; seq_index < static_cast<int>(aln.size()); seq_index++) {
        const std::string& record = aln[seq_index];
        int len_no_gaps = 0;
        for (char c : record) {
            if (c != '-') len_no_gaps++;
        }
        if (len_no_gaps < min_length) min_length = len_no_gaps;
        if (len_no_gaps > max_length) max_length = len_no_gaps;
        record_gap_lengths(record, seq_index, gap_positions, gaps_length_histogram);
    }
    calculate_counts(gap_positions);

    for (int pos = 0; pos < msa_length; pos++) {
        int num_gaps = 0;
        for (const auto& record : aln) {
            if (record[pos] == '-') num_gaps++;
        }
        if (num_gaps == 0) num_cols_no_gaps++;
        else if (num_gaps == 1) num_cols_1_gap++;
        else if (num_gaps == 2) num_cols_2_gaps++;
        else if (num_gaps == taxa_num - 1) num_cols_all_gaps_except1++;
    }

    seq_min_len = min_length;
    seq_max_len = max_length;
    gaps_len_one = gaps_length_histogram[1] / static_cast<double>(taxa_num);
    gaps_len_two = gaps_length_histogram[2] / static_cast<double>(taxa_num);
    gaps_len_three = gaps_length_histogram[3] / static_cast<double>(taxa_num);
    double four_plus = 0;
    for (const auto& [length, count] : gaps_length_histogram) {
        if (length > 3) four_plus += count;
    }
    gaps_len_four_plus = four_plus / taxa_num;
}

void GapStats::record_gap_lengths(const std::string& sequence, int seq_index,
                                    std::map<std::pair<int,int>, std::vector<int>>& gap_positions,
                                    std::map<int, int>& gaps_length_histogram) {
    int start_index = -1;
    int current_length = 0;
    int last_gap_index = -1;
    int sc_count = 0, dc_count = 0;

    for (int i = 0; i < static_cast<int>(sequence.size()); i++) {
        if (sequence[i] == '-') {
            if (current_length == 0) {
                start_index = i;
                if (start_index == last_gap_index + 2) sc_count++;
                else if (start_index == last_gap_index + 3) dc_count++;
            }
            current_length++;
        } else {
            if (current_length > 0) {
                auto key = std::make_pair(current_length, start_index);
                auto it = gap_positions.find(key);
                if (it != gap_positions.end()) {
                    bool found = false;
                    for (int si : it->second) { if (si == seq_index) { found = true; break; } }
                    if (!found) it->second.push_back(seq_index);
                } else {
                    gap_positions[key] = {seq_index};
                }
                gaps_length_histogram[current_length]++;
                last_gap_index = std::max(i - 1, 0);
            }
            current_length = 0;
        }
    }

    if (current_length > 0) {
        auto key = std::make_pair(current_length, start_index);
        auto it = gap_positions.find(key);
        if (it != gap_positions.end()) {
            bool found = false;
            for (int si : it->second) { if (si == seq_index) { found = true; break; } }
            if (!found) it->second.push_back(seq_index);
        } else {
            gap_positions[key] = {seq_index};
        }
        gaps_length_histogram[current_length]++;
    } else {
        int current_index = static_cast<int>(sequence.size());
        if (current_index == last_gap_index + 2) sc_count++;
        else if (current_index == last_gap_index + 3) dc_count++;
    }

    int total = 0;
    for (const auto& [len, cnt] : gaps_length_histogram) total += cnt;
    num_gap_segments = total;
    num_gap_segments_norm = static_cast<double>(num_gap_segments) / taxa_num;
    single_char_count += sc_count;
    double_char_count += dc_count;
}

void GapStats::calculate_counts(const std::map<std::pair<int,int>, std::vector<int>>& gap_positions) {
    std::map<int, std::map<int, int>> length_count = {{1, {}}, {2, {}}, {3, {}}};
    std::map<int, int> length_plus_count;
    int total_length = 0;
    int unique_gaps_cnt = 0;
    int unique_gaps_length = 0;

    for (const auto& [key, seq_set] : gap_positions) {
        int length = key.first;
        int num_seqs = static_cast<int>(seq_set.size());
        if (length >= 1 && length <= 3) {
            length_count[length][num_seqs]++;
        }
        if (length > 3) {
            length_plus_count[num_seqs]++;
        }
        if (num_seqs == 1) {
            unique_gaps_cnt++;
            unique_gaps_length += length;
        }
        total_length += num_seqs * length;
    }

    if (num_gap_segments > 0) {
        av_gap_segment_length = static_cast<double>(total_length) / num_gap_segments;
    }
    num_unique_gaps = unique_gaps_cnt;
    num_unique_gaps_norm = static_cast<double>(unique_gaps_cnt) / taxa_num;
    avg_unique_gap_length = static_cast<double>(unique_gaps_length) / std::max(unique_gaps_cnt, 1);
    gaps_1seq_len1 = length_count[1][1];
    gaps_2seq_len1 = length_count[1][2];
    gaps_all_except_1_len1 = length_count[1][taxa_num - 1];
    gaps_1seq_len2 = length_count[2][1];
    gaps_2seq_len2 = length_count[2][2];
    gaps_all_except_1_len2 = length_count[2][taxa_num - 1];
    gaps_1seq_len3 = length_count[3][1];
    gaps_2seq_len3 = length_count[3][2];
    gaps_all_except_1_len3 = length_count[3][taxa_num - 1];
    gaps_1seq_len4plus = length_plus_count[1];
    gaps_2seq_len4plus = length_plus_count[2];
    gaps_all_except_1_len4plus = length_plus_count[taxa_num - 1];
}

std::vector<StatValue> GapStats::get_my_features_as_list() const {
    return {
        StatValue(code),
        StatValue(num_gap_segments_norm), StatValue(av_gap_segment_length),
        StatValue(gaps_len_one), StatValue(gaps_len_two), StatValue(gaps_len_three), StatValue(gaps_len_four_plus),
        StatValue(num_unique_gaps), StatValue(num_unique_gaps_norm), StatValue(avg_unique_gap_length),
        StatValue(gaps_1seq_len1), StatValue(gaps_2seq_len1),
        StatValue(gaps_1seq_len2), StatValue(gaps_2seq_len2),
        StatValue(gaps_1seq_len3), StatValue(gaps_2seq_len3),
        StatValue(gaps_1seq_len4plus), StatValue(gaps_2seq_len4plus),
        StatValue(gaps_all_except_1_len1), StatValue(gaps_all_except_1_len2),
        StatValue(gaps_all_except_1_len3), StatValue(gaps_all_except_1_len4plus),
        StatValue(num_cols_no_gaps), StatValue(num_cols_1_gap),
        StatValue(num_cols_2_gaps), StatValue(num_cols_all_gaps_except1),
        StatValue(single_char_count), StatValue(double_char_count),
        StatValue(seq_max_len), StatValue(seq_min_len)
    };
}
