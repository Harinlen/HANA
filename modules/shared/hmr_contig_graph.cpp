#include <cassert>
#include <cstring>

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
    int32_t length, enzyme_count, name_size, name_buff_size = 0;
    char* name_buff = NULL;
    for (int32_t i = 0; i < contig_size; ++i)
    {
        //Read the size.
        fread(&length, sizeof(int32_t), 1, contig_file);
        fread(&enzyme_count, sizeof(int32_t), 1, contig_file);
        //Read the name.
        fread(&name_size, sizeof(int32_t), 1, contig_file);
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
        parser(length, enzyme_count, name_size, name_buff, user);
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

void hmr_graph_restore_nodes_size_proc(int32_t node_count, void* user)
{
    HMR_CONTIGS* contigs = reinterpret_cast<HMR_CONTIGS*>(user);
    //Reserve the contig count items.
    contigs->reserve(node_count);
}

//length, enzyme_count, name_size, name_buff, user
void hmr_graph_restore_nodes_node_proc(int32_t length, int32_t enzyme_count, int32_t name_size, char* name_buff, void* user)
{
    HMR_CONTIGS* contigs = reinterpret_cast<HMR_CONTIGS*>(user);
    //Create contig node info.
    char *contig_name = static_cast<char *>(malloc(name_size + 1));
    strncpy(contig_name, name_buff, name_size);
    contig_name[name_size] = '\0';
    contigs->push_back(HMR_CONTIG{ length, enzyme_count, name_size, contig_name });
}

void hmr_graph_restore_contig_data(const char* filepath, HMR_CONTIGS* contigs)
{
    hmr_graph_load_contigs(filepath, hmr_graph_restore_nodes_size_proc, hmr_graph_restore_nodes_node_proc, contigs);
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
        fwrite(&contig.enzyme_count, sizeof(int32_t), 1, contig_file);
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
    HMR_EDGE_INFO counter;
    for (uint64_t i = 0; i < edge_sizes; ++i)
    {
        //Read the edge counter structure.
        fread(&counter, sizeof(HMR_EDGE_INFO), 1, edge_file);
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

void hmr_graph_restore_edges_size_proc(uint64_t edge_sizes, void* user)
{
    HMR_EDGE_COUNTERS *edges = reinterpret_cast<HMR_EDGE_COUNTERS *>(user);
    //Reserve the items.
    edges->reserve(edge_sizes);
}

void hmr_graph_restore_edges_data_proc(const HMR_EDGE_INFO& edge_info, void* user)
{
    HMR_EDGE_COUNTERS *edges = reinterpret_cast<HMR_EDGE_COUNTERS *>(user);
    //Append the data.
    edges->push_back(edge_info);
}

void hmr_graph_restore_edges(const char *filepath, HMR_EDGE_COUNTERS *edges)
{
    hmr_graph_load_edge(filepath, hmr_graph_restore_edges_size_proc, hmr_graph_restore_edges_data_proc, edges);
}

void hmr_graph_restore_edge_map_size_proc(uint64_t edge_sizes, void* user)
{
    HMR_EDGE_MAP* edge_map = reinterpret_cast<HMR_EDGE_MAP*>(user);
    //Reserve the edge info.
    edge_map->reserve(edge_sizes);
}

void hmr_graph_restore_edge_map_data_proc(const HMR_EDGE_INFO& edge_info, void* user)
{
    HMR_EDGE_MAP* edge_map = reinterpret_cast<HMR_EDGE_MAP*>(user);
    //Append the edge information.
    edge_map->insert(std::make_pair(edge_info.edge.data, edge_info.weights));
}

void hmr_graph_restore_edge_map(const char* filepath, HMR_EDGE_MAP *edge_map)
{
    hmr_graph_load_edge(filepath, hmr_graph_restore_edge_map_size_proc, hmr_graph_restore_edge_map_data_proc, edge_map);
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
    for (const HMR_EDGE_INFO& edge : edges)
    {
        fwrite(&edge, sizeof(HMR_EDGE_INFO), 1, edge_file);
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

bool hmr_graph_save_partition(const char* filepath, const CONTIG_ID_VECTOR& contig_ids)
{
    FILE* group_id_file;
    if (!bin_open(filepath, &group_id_file, "wb"))
    {
        time_error(-1, "Failed to save chromosome partition ids to file %s", filepath);
        return false;
    }
    //Write the total groups.
    size_t group_size = contig_ids.size();
    fwrite(&group_size, sizeof(size_t), 1, group_id_file);
    //Write the allele group size.
    for (int32_t contig_id: contig_ids)
    {
        //Write the allele group size.
        fwrite(&contig_id, sizeof(int32_t), 1, group_id_file);
    }
    return true;
}
