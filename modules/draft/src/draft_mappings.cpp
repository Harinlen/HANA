#include "hmr_algorithm.hpp"
#include "hmr_contig_graph.hpp"

#include "draft_mappings.hpp"

void draft_mappings_build_edges(HMR_MAPPING* mapping, int32_t buf_size, void* user)
{
    EDGE_COUNTER* edges = static_cast<EDGE_COUNTER*>(user);
    //Loop for all the mapping info.
    for (int32_t i = 0; i < buf_size; ++i)
    {
        const HMR_MAPPING& mapping_info = mapping[i];
        //Only care about inter-connected edges.
        if(mapping_info.refID == mapping_info.next_refID)
        {
            continue;
        }
        //Increase the counter on the edge.
        uint64_t edge = hmr_graph_edge_data(mapping_info.refID, mapping_info.next_refID);
        const auto iter = edges->find(edge);
        if (iter == edges->end())
        {
            edges->insert(std::make_pair(edge, 1));
        }
        else
        {
            ++iter->second;
        }
    }
}

inline void count_edge(EDGE_COUNT_MAP* edge_map, int32_t node_start, int32_t node_end)
{
    const auto iter = edge_map->find(node_start);
    if(iter == edge_map->end())
    {
        //Create a node map for the node start.
        NODE_COUNT_MAP node_map;
        node_map.insert(std::make_pair(node_end, 1));
        edge_map->insert(std::make_pair(node_start, node_map));
    }
    else
    {
        //Find the end node in the target map.
        auto &node_map = iter->second;
        const auto node_iter = node_map.find(node_end);
        if(node_iter == node_map.end())
        {
            node_map.insert(std::make_pair(node_end, 1));
        }
        else
        {
            ++node_iter->second;
        }
    }
}

void draft_mappings_build_edge_counter(HMR_MAPPING *mapping, int32_t buf_size, void *user)
{
    EDGE_COUNT_MAP* edge_map = static_cast<EDGE_COUNT_MAP*>(user);
    //Loop for the all the mapping info.
    for(int32_t i = 0; i < buf_size; ++i)
    {
        const HMR_MAPPING& mapping_info = mapping[i];
        //Find the node start.
        count_edge(edge_map, mapping_info.refID, mapping_info.next_refID);
        count_edge(edge_map, mapping_info.next_refID, mapping_info.refID);
    }
}

inline bool remove_edge(EDGE_COUNT_MAP& edge_map, int32_t node_start, int32_t node_end)
{
    const auto iter = edge_map.find(node_start);
    if(iter == edge_map.end())
    {
        return false;
    }
    auto &node_map = iter->second;
    const auto node_iter = node_map.find(node_end);
    if(node_iter == node_map.end())
    {
        return false;
    }
    node_map.erase(node_iter);
    return true;
}

void draft_mappings_remove_edge(EDGE_COUNT_MAP &edge_map, int32_t node_a, int32_t node_b)
{
    if(remove_edge(edge_map, node_a, node_b))
    {
        remove_edge(edge_map, node_b, node_a);
    }
}
