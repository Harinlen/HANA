#ifndef ORDERING_TYPE_H
#define ORDERING_TYPE_H

#include "hmr_contig_graph_type.hpp"

typedef std::unordered_map<int32_t, int32_t> CONTIG_COUNTS;
typedef std::unordered_map<int32_t, CONTIG_COUNTS> ORDERING_COUNTS;

typedef std::unordered_map<int32_t, HMR_CONTIG> CONTIG_MAP;

typedef struct ORDERING_TIG ORDERING_TIG;

typedef struct ORDERING_INFO
{
    CONTIG_ID_VECTOR contig_group;
    CONTIG_MAP contigs;
    ORDERING_COUNTS edges;
    int32_t contig_size;

    ORDERING_TIG *init_genome;
    double best;
    uint64_t updated;
} ORDERING_INFO;

#endif // ORDERING_TYPE_H
