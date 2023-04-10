#ifndef DRAFT_MAPPINGS_TYPE_H
#define DRAFT_MAPPINGS_TYPE_H

#include "hmr_contig_graph_type.hpp"

typedef std::unordered_map<uint64_t, int32_t> EDGE_COUNTER;

typedef struct EDGE_BUILDER
{
    HMR_ALLELE_TABLE allele_table;
    EDGE_COUNTER counter;
} EDGE_BUILDER;

#endif // DRAFT_MAPPINGS_TYPE_H