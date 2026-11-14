#pragma once
#include <string>
#include <vector>
#include <variant>
#include <sstream>
#include <cmath>
#include <iomanip>

using StatValue = std::variant<int, double, std::string>;

class BasicStats {
public:
    std::string code;
    int taxa_num;
    int msa_length;
    std::vector<std::string> ordered_col_names;

    BasicStats(const std::string& code, int taxa_num, int msa_length,
               const std::vector<std::string>& ordered_col_names)
        : code(code), taxa_num(taxa_num), msa_length(msa_length),
          ordered_col_names(ordered_col_names) {}

    virtual ~BasicStats() = default;

    virtual std::vector<StatValue> get_my_features_as_list() const = 0;

    const std::vector<std::string>& get_ordered_col_names() const {
        return ordered_col_names;
    }

    static std::string stat_value_to_string(const StatValue& v) {
        return std::visit([](const auto& val) -> std::string {
            using T = std::decay_t<decltype(val)>;
            if constexpr (std::is_same_v<T, int>) {
                return std::to_string(val);
            } else if constexpr (std::is_same_v<T, double>) {
                std::ostringstream oss;
                oss << std::fixed << std::setprecision(3) << val;
                // Remove trailing zeros for clean output but keep at least one decimal
                std::string s = oss.str();
                return s;
            } else {
                return val;
            }
        }, v);
    }
};
