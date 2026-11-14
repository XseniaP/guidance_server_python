#include "enums.h"
#include <unordered_map>
#include <set>

using namespace std;

unordered_map<string, WeightMethods> str_to_enum_weight_methods = {
    {"HENIKOFF_WG", WeightMethods::HENIKOFF_WG},
    {"HENIKOFF_WOG", WeightMethods::HENIKOFF_WOG},
    {"CLUSTAL_MID_ROOT", WeightMethods::CLUSTAL_MID_ROOT},
    {"CLUSTAL_DIFFERENTIAL_SUM", WeightMethods::CLUSTAL_DIFFERENTIAL_SUM}
};

set<WeightMethods> strToEnumWeightMethods(const set<string>& weightMethodsSet) {
    set<WeightMethods> resultWeightMethodsSet;

    for (const auto& s : weightMethodsSet) {
        auto it = str_to_enum_weight_methods.find(s);
        if (it != str_to_enum_weight_methods.end()) {
            resultWeightMethodsSet.insert(it->second);
        }
    }
    return resultWeightMethodsSet;
}


unordered_map<string, StatsOutput> str_to_enum_stats_output = {
    {"SP", StatsOutput::SP},
    {"GAPS", StatsOutput::GAPS},
    {"ENTROPY", StatsOutput::ENTROPY},
    {"TREE", StatsOutput::TREE},
    {"K_MER", StatsOutput::K_MER},
    {"W_SP", StatsOutput::W_SP },
    {"SUBS_MATRIX", StatsOutput::SUBS_MATRIX },
    {"DISTANCE_LABELS", StatsOutput::DISTANCE_LABELS },
    {"ALL", StatsOutput::ALL},
    {"ALL_NO_SUBS_MATRIX", StatsOutput::ALL_NO_SUBS_MATRIX},
    {"ALL_FEATURES", StatsOutput::ALL_FEATURES}
};

set<StatsOutput> strToEnumStatsOutput(const set<string>& statsOutputSet) {
    set<StatsOutput> resultStatsOutputSet;

    for (const auto& s : statsOutputSet) {
        auto it = str_to_enum_stats_output.find(s);
        if (it != str_to_enum_stats_output.end()) {
            resultStatsOutputSet.insert(it->second);
        }
    }
    return resultStatsOutputSet;
}