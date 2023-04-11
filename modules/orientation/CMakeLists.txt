cmake_minimum_required(VERSION 3.0)

project(hana_orientation)

# Options
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Enable the lib includes.
include_directories(src)
include_directories(../shared/)

# Construct Binaries.
add_executable(hana_orientation
    ../shared/hmr_args.cpp
    ../shared/hmr_bin_file.cpp
    ../shared/hmr_contig_graph.cpp
    ../shared/hmr_path.cpp
    ../shared/hmr_ui.cpp
    src/args_orientation.cpp
    src/main.cpp
    src/orientation.cpp
)

