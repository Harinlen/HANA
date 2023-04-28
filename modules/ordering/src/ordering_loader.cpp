#include "hmr_algorithm.hpp"

#include "ordering_loader.hpp"

void ordering_edge_map_size_proc(uint64_t, void*)
{
}

inline void insert_edge(int32_t a, int32_t b, int32_t count, ORDERING_COUNTS& map)
{
    //Swap the b and a when necessary.
    if (a > b)
    {
        hmr_swap(a, b);
    }
    //Do the insertion.
    auto iter = map.find(a);
    if (iter == map.end())
    {
        CONTIG_COUNTS b_weights;
        b_weights.insert(std::make_pair(b, count));
        map.insert(std::make_pair(a, b_weights));
    }
    else
    {
        iter->second.insert(std::make_pair(b, count));
    }
}

void ordering_edge_map_data_proc(HMR_EDGE_INFO* edges, int32_t edge_size, void* user)
{
    ORDERING_EDGE_LOADER* loader = static_cast<ORDERING_EDGE_LOADER*>(user);
    //Build the edge map.
    for (int32_t i = 0; i < edge_size; ++i)
    {
        HMR_EDGE_INFO& edge_info = edges[i];
        //Only save the edges that inside the group.
        int32_t a_id = edge_info.start, b_id = edge_info.end;
        if (hmr_in_ordered_vector(a_id, loader->contig_group) && hmr_in_ordered_vector(b_id, loader->contig_group))
        {
            insert_edge(a_id, b_id, static_cast<int32_t>(edge_info.pairs), loader->edges);
        }
    }
}
