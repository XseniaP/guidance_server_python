cmake_minimum_required(VERSION 3.14)
project(sp_alternative_cpp CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES
    src/dist_labels_stats.cpp
    src/distance_calc.cpp
    src/entropy_stats.cpp
    src/enums.cpp
    src/gaps_stats.cpp
    src/global_alignment.cpp
    src/kmer_stats.cpp
    src/msa.cpp
    src/msa_align_by_pattern.cpp
    src/multi_msa_service.cpp
    src/neighbor_joining.cpp
    src/node.cpp
    src/rooted_tree.cpp
    src/sop_stats.cpp
    src/sp_score.cpp
    src/tree_stats.cpp
    src/unrooted_tree.cpp
    src/utils.cpp
    src/w_sop_stats.cpp
)

add_library(sp_alternative_cpp STATIC ${SOURCES})

target_include_directories(sp_alternative_cpp PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(APPLE)
    execute_process(COMMAND xcrun --sdk macosx --show-sdk-path
                    OUTPUT_VARIABLE MACOS_SDK_PATH OUTPUT_STRIP_TRAILING_WHITESPACE)
    target_include_directories(sp_alternative_cpp PUBLIC
        "${MACOS_SDK_PATH}/usr/include/c++/v1")
endif()
