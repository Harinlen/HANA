project(hana_draft)

# Options
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)

# Enable the lib includes.
include_directories(src)
include_directories(../shared/)

# Construct Binaries.
add_executable(hana_draft
    ../shared/hmr_algorithm.cpp
    ../shared/hmr_args.cpp
    ../shared/hmr_bin_file.cpp
    ../shared/hmr_contig_graph.cpp
    ../shared/hmr_path.cpp
    ../shared/hmr_ui.cpp
    src/args_draft.cpp
    src/draft_mappings.cpp
    src/main.cpp
)
target_link_libraries(hana_draft pthread)
