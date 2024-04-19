cmake_minimum_required(VERSION 3.27)
project(msa_set_score)

set(CMAKE_CXX_STANDARD 11)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O0")

# Source files
set(SOURCE_FILES
        msa_functions.cpp
        count_functions.cpp
        set_functions.cpp
        msa_set_score.cpp
        count_functions.cpp
)

# Executable
add_executable(msa_set_score ${SOURCE_FILES})
