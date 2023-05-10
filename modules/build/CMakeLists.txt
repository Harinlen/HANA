project(hana_build)

# Options
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Enable the lib includes.
include_directories(src)
include_directories(../shared/)

# zlib
find_package(ZLIB)

# Construct Binaries.
add_executable(hana_build
    ../shared/hmr_args.cpp
    ../shared/hmr_bin_file.cpp
    ../shared/hmr_bin_queue.cpp
    ../shared/hmr_contig_graph.cpp
    ../shared/hmr_fasta.cpp
    ../shared/hmr_gz.cpp
    ../shared/hmr_path.cpp
    ../shared/hmr_seq.cpp
    ../shared/hmr_text_file.cpp
    ../shared/hmr_ui.cpp
    src/args_build.cpp
    src/main.cpp
)
target_link_libraries(hana_build pthread z)
