cmake_minimum_required(VERSION 3.14)
project(features_for_msas CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Build the static library first
add_subdirectory(../sp_alternative_cpp ${CMAKE_BINARY_DIR}/sp_alternative_cpp)

add_executable(features_for_msas features_for_msas.cpp)

# Link against the static library
target_link_libraries(features_for_msas PRIVATE sp_alternative_cpp)

# nlohmann/json: use the local json.hpp via the nlohmann/ symlink in this directory
target_include_directories(features_for_msas PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
