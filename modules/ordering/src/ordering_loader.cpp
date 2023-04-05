#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>

#include "hmr_global.hpp"

#include "ordering_type.hpp"

#include "ordering_loader.hpp"

inline bool is_id_in_group(const CONTIG_ID_VECTOR &arr, int32_t contig_id)
{
    return std::find(arr.begin(), arr.end(), contig_id) != arr.end();
}

void ordering_contig_size_proc(int32_t node_count, void *user)
{
    ORDERING_INFO *ordering_info = reinterpret_cast<ORDERING_INFO *>(user);
    //Create the buffer for the node.
    ordering_info->contigs.reserve(node_count);
    ordering_info->contig_size = 0;
}

void ordering_contig_node_proc(int32_t length, int32_t enzyme_count, int32_t name_size, char *name_buff, void *user)
{
    ORDERING_INFO *ordering_info = reinterpret_cast<ORDERING_INFO *>(user);
    //Backup the contig index.
    int32_t contig_index = ordering_info->contig_size;
    ++ordering_info->contig_size;
    if(!is_id_in_group(ordering_info->contig_group, contig_index))
    {
        return;
    }
    //Create contig node info.
    char* contig_name = static_cast<char*>(malloc(static_cast<size_t>(name_size + 2)));
    assert(contig_name);
#ifdef _MSC_VER
    strncpy_s(contig_name, name_size + 1, name_buff, name_size);
#else
    strncpy(contig_name, name_buff, name_size);
#endif
    contig_name[name_size] = '\0';
    ordering_info->contigs.insert(std::make_pair(contig_index, HMR_CONTIG{ length, enzyme_count, name_size, contig_name }));
}

void ordering_edge_map_size_proc(uint64_t edge_sizes, void* user)
{
    HMR_UNUSED(edge_sizes)
    ORDERING_INFO *ordering_info = reinterpret_cast<ORDERING_INFO *>(user);
    //Reserve the edge info.
    ordering_info->edges.reserve(ordering_info->contig_size);
}

inline int32_t contig_enzyme_count(int32_t contig_id, const CONTIG_MAP &contigs)
{
    auto iter = contigs.find(contig_id);
    return iter == contigs.end() ? -1 : iter->second.enzyme_count;
}

inline int32_t contig_length(int32_t contig_id, const CONTIG_MAP &contigs)
{
    auto iter = contigs.find(contig_id);
    return iter == contigs.end() ? -1 : iter->second.length;
}

inline void insert_edge(int32_t a, int32_t b, double weight, ORDERING_COUNTS &map)
{
    //Swap the b and a when necessary.
    if(b < a)
    {
        int32_t temp = b;
        b = a;
        a = temp;
    }
    //Do the insertion.
    auto iter = map.find(a);
    if(iter == map.end())
    {
        CONTIG_COUNTS b_weights;
        b_weights.insert(std::make_pair(b, weight));
        map.insert(std::make_pair(a, b_weights));
    }
    else
    {
        if(iter->second.find(b) != iter->second.end())
        {
            return;
        }
        iter->second.insert(std::make_pair(b, weight));
    }
}

void ordering_edge_map_data_proc(const HMR_EDGE_INFO& edge_info, void* user)
{
    ORDERING_INFO *ordering_info = reinterpret_cast<ORDERING_INFO *>(user);
    //Append the edge information.
    int32_t a_id = edge_info.edge.pos.start, b_id = edge_info.edge.pos.end;
    if(!is_id_in_group(ordering_info->contig_group, a_id) ||
            !is_id_in_group(ordering_info->contig_group, b_id))
    {
        return;
    }
    //Insert the data to edge weights.
    insert_edge(a_id, b_id, edge_info.pairs, ordering_info->edges);
}
