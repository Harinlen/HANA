cmake_minimum_required(VERSION 3.0)

project(hana_ordering)

# Options
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Enable the lib includes.
include_directories(src)
include_directories(../shared/)

# Construct Binaries.
add_executable(hana_ordering
    ../shared/hmr_algorithm.cpp
    ../shared/hmr_args.cpp
    ../shared/hmr_bin_file.cpp
    ../shared/hmr_contig_graph.cpp
    ../shared/hmr_path.cpp
    ../shared/hmr_ui.cpp
    src/args_ordering.cpp
    src/main.cpp
    src/ordering_descent.cpp
    src/ordering_loader.cpp
)

