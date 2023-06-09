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

inline uint64_t hmr_graph_edge_data(int32_t edge_node_a, int32_t edge_node_b)
{
    uint64_t edge_a = static_cast<uint64_t>(edge_node_a),
            edge_b = static_cast<uint64_t>(edge_node_b);
    if (edge_node_a < edge_node_b)
    {
        return (edge_a << 32) | edge_b;
    }
    return (edge_b << 32) | edge_a;
}

inline void hmr_graph_edge_pos(uint64_t data, int32_t& start, int32_t& end)
{
    start = static_cast<int32_t>(data >> 32);
    end = static_cast<int32_t>(data & 0xFFFF);
}

/* Graph path generated functions */
std::string hmr_graph_path_contigs(const char* prefix);
std::string hmr_graph_path_edge(const char* prefix);
std::string hmr_graph_path_reads(const char* prefix);
std::string hmr_graph_path_nodes_invalid(const char* contig_path);
std::string hmr_graph_path_contigs_invalid(const char* prefix);
std::string hmr_graph_path_allele_table(const char* prefix);
std::string hmr_graph_path_cluster_name(const char* prefix, const int32_t index, const int32_t total);
std::string hmr_graph_path_chromo_name(const char* seq_path);

/* Node operations */
void hmr_graph_load_contigs(const char* filepath, HMR_NODES& nodes, HMR_NODE_NAMES *names = NULL);
bool hmr_graph_save_contigs(const char* filepath, const HMR_CONTIGS& nodes);

/* Node vector operations */
void hmr_graph_load_contig_ids(const char* filepath, HMR_CONTIG_ID_VEC& contig_ids);
bool hmr_graph_save_contig_ids(const char* filepath, const HMR_CONTIG_ID_VEC& contig_ids);

/* Node table operations */
void hmr_graph_load_contig_table(const char *filepath, HMR_CONTIG_ID_TABLE &contig_table);
bool hmr_graph_save_contig_table(const char *filepath, const HMR_CONTIG_ID_TABLE &contig_table);

/* Convert the contig table to allele map */
void hmr_graph_allele_map_init(HMR_ALLELE_MAP &allele_map, const HMR_CONTIG_ID_TABLE &allele_table);

/* Edge vector operations */
typedef void (*GRAPH_EDGE_SIZE_PROC)(uint64_t edge_size, void* user);
typedef void (*GRAPH_EDGE_PROC)(HMR_EDGE_INFO* edges, int32_t edge_size, void* user);
void hmr_graph_load_edges(const char* filepath, int32_t buf_size, GRAPH_EDGE_SIZE_PROC size_proc, GRAPH_EDGE_PROC proc, void *user);
bool hmr_graph_save_edges(const char* filepath, const HMR_EDGE_COUNTERS& edges);

/* Paired-reads operations */
typedef void (*HMR_READS_PROC)(HMR_MAPPING* mapping, int32_t buf_size, void *user);
void hmr_graph_load_reads(const char* filepath, int32_t buf_size, HMR_READS_PROC proc, void *user);

/* Chromosome sequence operations */
void hmr_graph_load_chromosome(const char* filepath, CHROMOSOME_CONTIGS& seq);
bool hmr_graph_save_chromosome(const char* filepath, const CHROMOSOME_CONTIGS& seq);

#endif // HMR_CONTIG_GRAPH_H
