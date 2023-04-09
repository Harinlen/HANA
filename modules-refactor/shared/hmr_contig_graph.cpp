#include <cassert>
#include <cstring>

#include "hmr_bin_file.hpp"
#include "hmr_global.hpp"
#include "hmr_ui.hpp"
#include "hmr_path.hpp"

#include "hmr_contig_graph.hpp"

std::string hmr_graph_path_contigs(const char* prefix)
{
    return std::string(prefix) + ".hmr_contigs";
}

std::string hmr_graph_path_reads(const char* prefix)
{
    return std::string(prefix) + ".hmr_reads";
}

std::string hmr_graph_path_edge(const char* prefix)
{
    return std::string(prefix) + ".hmr_edges";
}

std::string hmr_graph_path_contigs_invalid(const char* prefix)
{
    return std::string(prefix) + ".hmr_contigs_invalid";
}

std::string hmr_graph_path_allele_table(const char* prefix)
{
    return std::string(prefix) + ".hmr_allele_table";
}

bool hmr_graph_save_contigs(const char* filepath, const HMR_CONTIGS& contigs)
{
    //Load the file path.
    FILE* contig_file;
    if (!bin_open(filepath, &contig_file, "wb"))
    {
        time_error(-1, "Failed to save contig file %s", filepath);
    }
    //Write the file of the size.
    int32_t contig_size = static_cast<int32_t>(contigs.size());
    fwrite(&contig_size, sizeof(int32_t), 1, contig_file);
    //Now write contig information.
    for (int32_t i = 0; i < contig_size; ++i)
    {
        const HMR_CONTIG& contig = contigs[i];
        //Write the contig size.
        fwrite(&contig.length, sizeof(int32_t), 1, contig_file);
        fwrite(&contig.enzyme_count, sizeof(int32_t), 1, contig_file);
        //Write the contig name.
        fwrite(&contig.name_size, sizeof(int32_t), 1, contig_file);
        fwrite(contig.name, sizeof(char), contig.name_size, contig_file);
    }
    //Close the file.
    fclose(contig_file);
    return true;
}

bool hmr_graph_save_contig_ids(const char* filepath, const HMR_CONTIG_ID_VEC& contig_ids)
{
    FILE* contig_ids_file;
    if (!bin_open(filepath, &contig_ids_file, "wb"))
    {
        time_error(-1, "Failed to save contig ids to file %s", filepath);
        return false;
    }
    //Write the number of invalid ids.
    int32_t id_sizes = static_cast<int32_t>(contig_ids.size());
    fwrite(&id_sizes, sizeof(int32_t), 1, contig_ids_file);
    for (int32_t contig_id : contig_ids)
    {
        fwrite(&contig_id, sizeof(int32_t), 1, contig_ids_file);
    }
    fclose(contig_ids_file);
    return true;
}

bool hmr_graph_save_allele_table(const char* filepath, const HMR_ALLELE_TABLE& allele_table)
{
    FILE* allele_table_file;
    if (!bin_open(filepath, &allele_table_file, "wb"))
    {
        time_error(-1, "Failed to save allele table to file %s", filepath);
        return false;
    }
    //Write the record counts.
    int32_t record_sizes = static_cast<int32_t>(allele_table.size());
    fwrite(&record_sizes, sizeof(int32_t), 1, allele_table_file);
    for (const auto &record : allele_table)
    {
        fwrite(&record.first, sizeof(int32_t), 1, allele_table_file);
        int32_t record_size = static_cast<int32_t>(record.second.size());
        fwrite(&record_size, sizeof(int32_t), 1, allele_table_file);
        for (const int32_t conflict_id : record.second)
        {
            fwrite(&conflict_id, sizeof(int32_t), 1, allele_table_file);
        }
    }
    return true;
}
