cmake_minimum_required(VERSION 2.8)

project(hana_partition)

# Options
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Enable the lib includes.
include_directories(src)
include_directories(../shared/)

# Construct Binaries.
add_executable(hana_partition
    ../shared/hmr_args.cpp
    ../shared/hmr_bin_file.cpp
    ../shared/hmr_contig_graph.cpp
    ../shared/hmr_path.cpp
    ../shared/hmr_ui.cpp
    src/args_partition.cpp
    src/main.cpp
    src/partition.cpp
)
target_link_libraries(hana_partition pthread)

