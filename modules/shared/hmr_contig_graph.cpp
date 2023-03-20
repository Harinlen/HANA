#include <cassert>

#include "hmr_bin_file.hpp"
#include "hmr_ui.hpp"
#include "hmr_path.hpp"

#include "hmr_contig_graph.hpp"

std::string hmr_graph_path_contig(const char* prefix)
{
    return std::string(prefix) + ".hmr_contig";
}

void hmr_graph_load_contigs_bin(const char* filepath, CONTIG_SIZE_PROC size_parser, CONTIG_PROC parser, void* user)
{
    //Load the file path.
    FILE* contig_file;
    if (!bin_open(filepath, &contig_file, "rb"))
    {
        time_error(-1, "Failed to open contig file %s", filepath);
    }
    //Read the file of the size.
    int32_t contig_size;
    fread(&contig_size, sizeof(int32_t), 1, contig_file);
    size_parser(contig_size, user);
    //Now parsing the parser.
    int32_t length, name_size, name_buff_size = 0;
    char* name_buff = NULL;
    for (int32_t i = 0; i < contig_size; ++i)
    {
        //Read the size.
        fread(&length, sizeof(int32_t), 1, contig_file);
        fread(&name_size, sizeof(int32_t), 1, contig_file);
        //Allocate the name buffer.
        if (name_buff_size < name_size)
        {
            if (name_buff)
            {
                char* update_buff = static_cast<char*>(realloc(name_buff, name_size + 1));
                assert(update_buff);
                name_buff = update_buff;
            }
            else
            {
                name_buff = static_cast<char*>(malloc(name_size + 1));
            }
            name_buff_size = name_size + 1;
        }
        //Read the name.
        assert(name_buff);
        fread(name_buff, sizeof(char), name_size, contig_file);
        //Call the parser function.
        parser(length, name_size, name_buff, user);
    }
    //Close the file.
    fclose(contig_file);
    //Clear the buff.
    if (name_buff)
    {
        free(name_buff);
    }
}

void hmr_graph_load_contigs_text(const char* filepath, CONTIG_SIZE_PROC size_parser, CONTIG_PROC parser, void* user)
{
    //Load the file path.
}

void hmr_graph_load_contigs(const char* filepath, CONTIG_SIZE_PROC size_parser, CONTIG_PROC parser, void* user, bool binary_format)
{
    if (binary_format)
    {
        hmr_graph_load_contigs_bin(filepath, size_parser, parser, user);
    }
    else
    {
        //hmr_graph_load_contigs_text(filepath, size_parser, parser, user);
        time_error(-1, "Not supporting loading text format.");
    }
}

bool hmr_graph_save_contigs_bin(const char* filepath, const HMR_CONTIGS& contigs)
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
        //Write the contig name.
        fwrite(&contig.name_size, sizeof(int32_t), 1, contig_file);
        fwrite(contig.name, sizeof(char), contig.name_size, contig_file);
    }
    //Close the file.
    fclose(contig_file);
    return true;
}

bool hmr_graph_save_contigs_text(const char* filepath, const HMR_CONTIGS& contigs)
{
    time_error(-1, "Not supporting saving text format.");
    return false;
}

bool hmr_graph_save_contigs(const char* filepath, const HMR_CONTIGS& contigs, bool binary_format)
{
    //Save the contig information.
    if (binary_format)
    {
        //Save as binary file, for better performance.
        return hmr_graph_save_contigs_bin(filepath, contigs);
    }
    //Save as text file for debugging and editing.
    return hmr_graph_save_contigs_text(filepath, contigs);
}

std::string hmr_graph_path_reads(const char* prefix)
{
    return std::string(prefix) + ".hmr_reads";
}

std::string hmr_graph_path_edge(const char* prefix)
{
    return std::string(prefix) + ".hmr_edge";
}

void hmr_graph_load_edge_bin(const char* filepath, EDGE_SIZE_PROC size_parser, EDGE_PROC parser, void* user)
{
    //Load the file path.
    FILE* edge_file;
    if (!bin_open(filepath, &edge_file, "rb"))
    {
        time_error(-1, "Failed to open edge information file %s", filepath);
    }
    //Read the edge size information.
    uint64_t edge_sizes = 0;
    fread(&edge_sizes, sizeof(uint64_t), 1, edge_file);
    size_parser(edge_sizes, user);
    HMR_EDGE_COUNTER counter;
    for (uint64_t i = 0; i < edge_sizes; ++i)
    {
        //Read the edge counter structure.
        fread(&counter, sizeof(HMR_EDGE_COUNTER), 1, edge_file);
        //Call the parser.
        parser(counter, user);
    }
    //Close the file.
    fclose(edge_file);
}

void hmr_graph_load_edge(const char* filepath, EDGE_SIZE_PROC size_parser, EDGE_PROC parser, void* user, bool binary_format)
{
    if (binary_format)
    {
        hmr_graph_load_edge_bin(filepath, size_parser, parser, user);
    }
    else
    {
        //hmr_graph_load_contigs_text(filepath, size_parser, parser, user);
        time_error(-1, "Not supporting loading text format.");
    }
}

bool hmr_graph_save_edge(const char* filepath, const HMR_EDGE_COUNTERS& edges)
{
    //Open the contig output file to write the data.
    FILE* edge_file;
    if (!bin_open(filepath, &edge_file, "wb"))
    {
        time_error(-1, "Failed to save edge information from %s", filepath);
        return false;
    }
    //Write the edge information.
    uint64_t edge_sizes = static_cast<uint64_t>(edges.size());
    fwrite(&edge_sizes, sizeof(uint64_t), 1, edge_file);
    for (const HMR_EDGE_COUNTER& edge : edges)
    {
        fwrite(&edge, sizeof(HMR_EDGE_COUNTER), 1, edge_file);
    }
    fclose(edge_file);
    return true;
}

std::string hmr_graph_path_invalid(const char* prefix)
{
    return std::string(prefix) + ".hmr_invalid";
}

bool hmr_graph_save_invalid_bin(const char* filepath, const HMR_CONTIG_INVALID_IDS& ids)
{
    FILE* ids_file;
    if (!bin_open(filepath, &ids_file, "wb"))
    {
        time_error(-1, "Failed to save invalid contig ids to file %s", filepath);
        return false;
    }
    //Write the number of invalid ids.
    size_t id_sizes = ids.size();
    fwrite(&id_sizes, sizeof(size_t), 1, ids_file);
    for (const int32_t& id : ids)
    {
        fwrite(&id, sizeof(int32_t), 1, ids_file);
    }
    fclose(ids_file);
    return true;
}

bool hmr_graph_save_invalid(const char* filepath, const HMR_CONTIG_INVALID_IDS& ids, bool binary_format)
{
    //Save the contig information.
    if (binary_format)
    {
        //Save as binary file, for better performance.
        return hmr_graph_save_invalid_bin(filepath, ids);
    }
    //Save as text file for debugging and editing.
    return false;
}

bool hmr_graph_save_partition(const char* filepath, const CHROMOSOME_GROUP& group_ids)
{
    FILE* group_id_file;
    if (!bin_open(filepath, &group_id_file, "wb"))
    {
        time_error(-1, "Failed to save chromosome partition ids to file %s", filepath);
        return false;
    }
    //Write the total groups.
    size_t group_size = group_ids.size();
    fwrite(&group_size, sizeof(size_t), 1, group_id_file);
    //Write the allele group size.
    for (const auto& allele_groups : group_ids)
    {
        //Write the allele group size.
        group_size = allele_groups.size();
        fwrite(&group_size, sizeof(size_t), 1, group_id_file);
        for (const auto& contig_ids : allele_groups)
        {
            //Write the number of contigs.
            group_size = contig_ids.size();
            fwrite(&group_size, sizeof(size_t), 1, group_id_file);
            for (const int32_t& contig_id : contig_ids)
            {
                fwrite(&contig_id, sizeof(int32_t), 1, group_id_file);
            }
        }
    }
    return true;
}
