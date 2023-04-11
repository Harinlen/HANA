#ifndef PARTITION_TYPE_H
#define PARTITION_TYPE_H

#include <cstdlib>

#include "hmr_contig_graph_type.hpp"
#include "hmr_allele_type.hpp"

typedef struct MERGE_OP
{
    CONTIG_ID_VECTOR *a;
    CONTIG_ID_VECTOR *b;
    double weight;
    bool is_valid;
} MERGE_OP;

typedef std::unordered_map<int32_t, double> CONTIG_LINK_DENSITY;
typedef std::vector<CONTIG_LINK_DENSITY> MAP_LINK_DENSITY;

typedef struct LINK_FACTORS
{
    double* factors;
    double avergage_links;
    CONTIG_ID_SET skipped;
    const bool allele_mode;
    const HMR_ALLELE_TABLE &allele_table;
} LINK_FACTORS;

typedef struct CLUSTER_INFO
{
    CONTIG_ID_VECTOR** belongs = NULL;
    CONTIG_ID_VECTOR** clusters = NULL;
    size_t cluster_size = 0;
    MAP_LINK_DENSITY link_densities;
    MERGE_OP *merge = NULL;
    size_t merge_size = 0;
    bool allele_mode = false;
    HMR_ALLELE_TABLE allele_table;
} CLUSTER_INFO;

typedef struct DENSITY_INFO
{
    LINK_FACTORS& link;
    CLUSTER_INFO& cluster;
} DENSITY_INFO;

#endif // PARTITION_TYPE_H
