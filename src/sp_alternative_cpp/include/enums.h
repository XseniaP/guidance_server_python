#pragma once
#include <string>
#include <set>

using namespace std;

enum class SopCalcTypes {
    NAIVE = 0,
    EFFICIENT = 1
};

enum class RootingMethods {
    LONGEST_PATH_MID,
    MIN_DIFFERENTIAL_SUM
};

enum class WeightMethods {
    HENIKOFF_WG,
    HENIKOFF_WOG,
    CLUSTAL_MID_ROOT,
    CLUSTAL_DIFFERENTIAL_SUM
};

set<WeightMethods> strToEnumWeightMethods(const set<string>& strings);

enum class AffineGapMatrixTypes {
    M,
    X,
    Y
};

enum class StatsOutput {
    SP,
    GAPS,
    ENTROPY,
    TREE,
    K_MER,
    W_SP,
    SUBS_MATRIX,
    DISTANCE_LABELS,
    ALL,
    ALL_FEATURES,
    ALL_NO_SUBS_MATRIX,
};

set<StatsOutput> strToEnumStatsOutput(const set<string>& strings);

enum class DistanceType {
    D_SSP,
    D_SEQ,
    D_POS
};

inline std::string stats_output_to_string(StatsOutput s) {
    switch (s) {
        case StatsOutput::SP: return "sop";
        case StatsOutput::GAPS: return "gaps";
        case StatsOutput::ENTROPY: return "entropy";
        case StatsOutput::TREE: return "tree";
        case StatsOutput::K_MER: return "k_mer";
        case StatsOutput::W_SP: return "w_sop";
        case StatsOutput::SUBS_MATRIX: return "subs_matrix";
        case StatsOutput::DISTANCE_LABELS: return "distances";
        case StatsOutput::ALL: return "all_stats";
        case StatsOutput::ALL_FEATURES: return "all_features";
        case StatsOutput::ALL_NO_SUBS_MATRIX: return "all_no_subs_matrix";
    }
    return "";
}

inline std::string weight_method_to_string(WeightMethods w) {
    switch (w) {
        case WeightMethods::HENIKOFF_WG: return "henikoff_with_gaps";
        case WeightMethods::HENIKOFF_WOG: return "henikoff_without_gaps";
        case WeightMethods::CLUSTAL_MID_ROOT: return "clustal_mid_root";
        case WeightMethods::CLUSTAL_DIFFERENTIAL_SUM: return "clustal_differential_sum";
    }
    return "";
}

inline std::string rooting_method_to_string(RootingMethods r) {
    switch (r) {
        case RootingMethods::LONGEST_PATH_MID: return "longest_path_mid";
        case RootingMethods::MIN_DIFFERENTIAL_SUM: return "min_differential_sum";
    }
    return "";
}
