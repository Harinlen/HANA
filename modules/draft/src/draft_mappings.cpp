#include "hmr_algorithm.hpp"
#include "hmr_contig_graph.hpp"

#include "draft_mappings.hpp"

inline void draft_mapping_count_edge(EDGE_COUNTER& counter, uint64_t edge)
{
    const auto iter = counter.find(edge);
    if (iter == counter.end())
    {
        counter.insert(std::make_pair(edge, 1));
    }
    else
    {
        ++iter->second;
    }
}

void draft_mappings_build_edges(HMR_MAPPING* mapping, int32_t buf_size, void* user)
{
    EDGE_COUNTER* edges = static_cast<EDGE_COUNTER*>(user);
    //Loop for all the mapping info.
    for (int32_t i = 0; i < buf_size; ++i)
    {
        const HMR_MAPPING& mapping_info = mapping[i];
        //Increase the counter on the edge.
        HMR_EDGE edge = hmr_graph_edge(mapping_info.refID, mapping_info.next_refID);
        draft_mapping_count_edge(*edges, edge.data);
    }
}
