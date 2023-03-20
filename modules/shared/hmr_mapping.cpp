#include "hmr_path.hpp"
#include "hmr_ui.hpp"
#include "hmr_bam.hpp"

#include "hmr_mapping.hpp"

void hmr_mapping_read(const char* filepath, MAPPING_PROC proc, void* user, int threads)
{
    // Based on the suffix of the file path, decide how to load the file.
    std::string mapping_suffix = path_suffix(filepath);
    void* mapping_handle = NULL;
    if (mapping_suffix == ".bam")
    {
        //Read the file as bam.
        hmr_bam_read(filepath, proc, user, threads);
    }
    else
    {
        time_error(-1, "Unknown mapping file suffix: %s", mapping_suffix.data());
    }
}
