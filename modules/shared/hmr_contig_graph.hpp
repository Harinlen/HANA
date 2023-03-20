#ifndef HMR_CONTIG_GRAPH_H
#define HMR_CONTIG_GRAPH_H

#include <string>

#include "hmr_contig_graph_type.hpp"

/* Construct edge between node ID */
inline HMR_EDGE hmr_graph_edge(int32_t edge_node_a, int32_t edge_node_b)
{
    HMR_EDGE edge{};
    if (edge_node_a < edge_node_b)
    {
        edge.pos.start = edge_node_a;
        edge.pos.end = edge_node_b;
    }
    else
    {
        edge.pos.start = edge_node_b;
        edge.pos.end = edge_node_a;
    }
    return edge;
}

/* Graph node loader */
std::string hmr_graph_path_contig(const char* prefix);
typedef void (*CONTIG_SIZE_PROC)(int32_t, void*);
typedef void (*CONTIG_PROC)(int32_t, int32_t, char* , void*);
void hmr_graph_load_contigs(const char* filepath, CONTIG_SIZE_PROC size_parser, CONTIG_PROC parser, void *user, bool binary_format = true);
bool hmr_graph_save_contigs(const char* filepath, const HMR_CONTIGS& contigs, bool binary_format = true);

std::string hmr_graph_path_reads(const char* prefix);

std::string hmr_graph_path_edge(const char* prefix);
typedef void (*EDGE_SIZE_PROC)(uint64_t, void*);
typedef void (*EDGE_PROC)(const HMR_EDGE_COUNTER &, void*);
void hmr_graph_load_edge(const char* filepath, EDGE_SIZE_PROC size_parser, EDGE_PROC parser, void *user, bool binary_format = true);
bool hmr_graph_save_edge(const char* filepath, const HMR_EDGE_COUNTERS& edges);

std::string hmr_graph_path_invalid(const char* prefix);
bool hmr_graph_save_invalid(const char* filepath, const HMR_CONTIG_INVALID_IDS& ids, bool binary_format = true);

bool hmr_graph_save_partition(const char* filepath, const CHROMOSOME_GROUP& group_ids);

#endif // HMR_CONTIG_GRAPH_H