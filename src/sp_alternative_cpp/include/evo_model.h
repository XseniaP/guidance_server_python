#pragma once
#include <string>

struct EvoModel {
    double go_cost;
    double ge_cost;
    std::string matrix_file_name;
    std::string name;

    EvoModel(double go_cost, double ge_cost, const std::string& matrix_file_name,
             const std::string& name = "")
        : go_cost(go_cost), ge_cost(ge_cost), matrix_file_name(matrix_file_name),
          name(name.empty() ? matrix_file_name : name) {}
};
