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
std::string hmr_graph_path_nodes_invalid(const char* contig_path);
std::string hmr_graph_path_contigs_invalid(const char* prefix);
std::string hmr_graph_path_allele_table(const char* prefix);
std::string hmr_graph_path_cluster_name(const char* prefix, const int32_t index, const int32_t total);

/* Node operations */
void hmr_graph_load_contigs(const char* filepath, HMR_NODES& nodes, HMR_NODE_NAMES *names = NULL);
bool hmr_graph_save_contigs(const char* filepath, const HMR_CONTIGS& nodes);

/* Node vector operations */
void hmr_graph_load_contig_ids(const char* filepath, HMR_CONTIG_ID_VEC& contig_ids);
bool hmr_graph_save_contig_ids(const char* filepath, const HMR_CONTIG_ID_VEC& contig_ids);

/* Edge vector operations */
typedef void (*GRAPH_EDGE_SIZE_PROC)(uint64_t edge_size, void* user);
typedef void (*GRAPH_EDGE_PROC)(HMR_EDGE_INFO* edges, int32_t edge_size, void* user);
void hmr_graph_load_edges(const char* filepath, int32_t buf_size, GRAPH_EDGE_SIZE_PROC size_proc, GRAPH_EDGE_PROC proc, void *user);
bool hmr_graph_save_edges(const char* filepath, const HMR_EDGE_COUNTERS& edges);

/* Allele table operations */
bool hmr_graph_allele_conflict(const HMR_ALLELE_TABLE& allele_table, int32_t contig_id_a, int32_t contig_id_b);
void hmr_graph_load_allele_table(const char* filepath, HMR_ALLELE_TABLE& allele_table);
bool hmr_graph_save_allele_table(const char* filepath, const HMR_ALLELE_TABLE& allele_table);

/* Paired-reads operations */
typedef void (*HMR_READS_PROC)(HMR_MAPPING* mapping, int32_t buf_size, void *user);
void hmr_graph_load_reads(const char* filepath, int32_t buf_size, HMR_READS_PROC proc, void *user);

#endif // HMR_CONTIG_GRAPH_H
