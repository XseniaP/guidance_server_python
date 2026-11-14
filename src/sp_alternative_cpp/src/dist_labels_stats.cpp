#include "dist_labels_stats.h"
#include "distance_calc.h"
#include <iostream>

DistanceLabelsStats::DistanceLabelsStats(const std::string& code, int taxa_num, int msa_length)
    : BasicStats(code, taxa_num, msa_length,
                  {"code", "ssp_from_true", "dseq_from_true", "dpos_from_true"}),
      ssp_from_true(-1), dseq_from_true(-1), dpos_from_true(-1), rf_from_true(-1) {}

void DistanceLabelsStats::set_my_distance_from_true(const std::vector<std::string>& inferred_msa,
                                                    const std::vector<std::string>& true_msa){
    //dpos_from_true = compute_distance(true_msa, inferred_msa, DistanceType::D_POS);
    //ssp_from_true = compute_distance(true_msa, inferred_msa, DistanceType::D_SSP);
    dpos_from_true = 0.0;
    ssp_from_true = 0.0;
	const int rows_num = true_msa.size();
	const int cols_num = true_msa[0].size();
    vector<vector<int>> true_msa_vectors(cols_num, std::vector<int>(rows_num));
    vector<long> true_msa_map;
    cout << "Filling true vector and map..." << endl;
    fill_d_seq_vectors(true_msa, true_msa_vectors, true_msa_map, rows_num, cols_num);
    cout << "End filling true vector and map..." << endl;

    dseq_from_true = compute_eff_d_seq_from_true(inferred_msa, true_msa_vectors, true_msa_map);
	//dseq_from_true = compute_distance_from_known_msa(inferred_msa, profile_b_h, DistanceType::D_SEQ);  // older version of dseq_from_true, to compare with the efficient version
 //   dseq_from_true = compute_eff_d_seq(inferred_msa, true_msa); // older version of dseq_from_true, to compare with the efficient version

}

void DistanceLabelsStats::set_rf_from_true(const UnrootedTree& my_tree, const UnrootedTree& true_tree) {
    rf_from_true = my_tree.calc_rf(true_tree);
}

std::pair<std::vector<StatValue>, std::vector<std::string>> DistanceLabelsStats::get_print_rf() const {
    return {{StatValue(code), StatValue(rf_from_true)}, {"code", "rf_from_true"}};
}

std::vector<StatValue> DistanceLabelsStats::get_my_features_as_list() const {
    return {
        StatValue(code),
        StatValue(ssp_from_true), StatValue(dseq_from_true), StatValue(dpos_from_true)
    };
}
