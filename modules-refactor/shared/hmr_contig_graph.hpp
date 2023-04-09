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

/* Graph path generated functions */
std::string hmr_graph_path_contigs(const char* prefix);
std::string hmr_graph_path_edge(const char* prefix);
std::string hmr_graph_path_reads(const char* prefix);
std::string hmr_graph_path_contigs_invalid(const char* prefix);
std::string hmr_graph_path_allele_table(const char* prefix);

/* Node operations */
bool hmr_graph_save_contigs(const char* filepath, const HMR_CONTIGS& contigs);

/* Node vector operations */
bool hmr_graph_save_contig_ids(const char* filepath, const HMR_CONTIG_ID_VEC& contig_ids);

/* Allele table operations */
bool hmr_graph_save_allele_table(const char* filepath, const HMR_ALLELE_TABLE& allele_table);

#endif // HMR_CONTIG_GRAPH_H
