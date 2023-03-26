#ifndef PARTITION_TYPE_H
#define PARTITION_TYPE_H

#include "hmr_contig_graph_type.hpp"

typedef struct MERGE_OP
{
    CONTIG_ID_VECTOR *a;
    CONTIG_ID_VECTOR *b;
    double weight;
    bool is_valid;
} MERGE_OP;

typedef std::unordered_map<int32_t, double> CONTIG_LINK_DENSITY;
typedef std::unordered_map<int32_t, CONTIG_LINK_DENSITY> MAP_LINK_DENSITY;

typedef struct CLUSTER_INFO
{
    CONTIG_ID_VECTOR** belongs;
    CONTIG_ID_VECTOR** clusters;
    size_t cluster_size;
    MAP_LINK_DENSITY link_densities;
    MERGE_OP *merge;
    size_t merge_size;
} CLUSTER_INFO;

#endif // PARTITION_TYPE_H
