//#include "config.h"
//#include "evo_model.h"
//#include "enums.h"
//#include <iostream>
//#include <filesystem>
//#include <string>
//
//// Forward declaration from multi_msa_service.cpp
//void multiple_msa_calc_features_and_labels(const Configuration& config,
//    const std::string& matrix_base_path);
//
//int main(int argc, char* argv[]) {
//    // Default configuration - matching the Python project's typical usage
//    std::string input_dir = "input_comparison_files";
//    std::string matrix_base_path = "../input_config_files";
//
//    if (argc > 1) {
//        input_dir = argv[1];
//    }
//    if (argc > 2) {
//        matrix_base_path = argv[2];
//    }
//
//    // Create EvoModel(s) - BLOSUM62 with typical gap costs
//    EvoModel blosum62(-11, -1.0, "BLOSUM62", "BLOSUM62");
//
//    // Build configuration
//    Configuration config(
//        { blosum62 },                          // models
//        SopCalcTypes::EFFICIENT,             // sop_calc_type
//        input_dir,                           // input_files_dir_name
//        { WeightMethods::HENIKOFF_WG,         // additional_weights
//         WeightMethods::HENIKOFF_WOG,
//         WeightMethods::CLUSTAL_MID_ROOT,
//         WeightMethods::CLUSTAL_DIFFERENTIAL_SUM },
//        { 10, 20 },                            // k_values
//        { StatsOutput::ALL }                   // stats_output
//    );
//
//    std::cout << "SP Alternative (C++) - MSA Analysis Tool" << std::endl;
//    std::cout << "Input directory: " << input_dir << std::endl;
//    std::cout << "Matrix path: " << matrix_base_path << std::endl;
//
//    try {
//        multiple_msa_calc_features_and_labels(config, matrix_base_path);
//    }
//    catch (const std::exception& e) {
//        std::cerr << "Error: " << e.what() << std::endl;
//        return 1;
//    }
//
//    return 0;
//}
