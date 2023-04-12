#ifndef PARTITION_TYPE_H
#define PARTITION_TYPE_H

#include <cstdlib>

#include "hmr_contig_graph_type.hpp"

typedef std::unordered_map<int32_t, double> CONTIG_LINK_DENSITY;
typedef std::vector<CONTIG_LINK_DENSITY> GRAPH_LINK_DENSITY;

typedef struct CLUSTER_MERGE_OP
{
    HMR_CONTIG_ID_VEC *a;
    HMR_CONTIG_ID_VEC* b;
    double score;
    bool is_valid;
} CLUSTER_MERGE_OP;

typedef struct CLUSTER_INFO
{
    HMR_CONTIG_ID_VEC** belongs = NULL;
    HMR_CONTIG_ID_VEC** clusters = NULL;
    size_t cluster_size = 0;
    CLUSTER_MERGE_OP* merge = NULL;
    size_t merge_size = 0;
    HMR_ALLELE_MAP* allele_map = NULL;
    GRAPH_LINK_DENSITY link_densities;
} CLUSTER_INFO;

#endif // PARTITION_TYPE_H
