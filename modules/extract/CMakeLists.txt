project(hana_extract)

# Options
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Enable the lib includes.
include_directories(src)
include_directories(../shared/)

# zlib
# find_package(ZLIB)

# Construct Binaries.
add_executable(hana_extract
    ../shared/hmr_algorithm.cpp
    ../shared/hmr_args.cpp
    ../shared/hmr_bam.cpp
    ../shared/hmr_bgzf.cpp
    ../shared/hmr_bin_file.cpp
    ../shared/hmr_bin_queue.cpp
    ../shared/hmr_contig_graph.cpp
    ../shared/hmr_enzyme.cpp
    ../shared/hmr_fasta.cpp
    ../shared/hmr_gz.cpp
    ../shared/hmr_pairs.cpp
    ../shared/hmr_path.cpp
    ../shared/hmr_seq.cpp
    ../shared/hmr_text_file.cpp
    ../shared/hmr_ui.cpp
    src/args_extract.cpp
    src/extract_allele.cpp
    src/extract_fasta.cpp
    src/extract_index_map.cpp
    src/extract_mapping.cpp
    src/main.cpp
)
target_link_libraries(hana_extract pthread z)
