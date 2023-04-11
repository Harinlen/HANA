#ifndef ORDERING_TYPE_H
#define ORDERING_TYPE_H

#include "hmr_contig_graph_type.hpp"

typedef struct ORDERING_TIG
{
    int32_t index;
    int32_t length;
} ORDERING_TIG;

typedef std::unordered_map<int32_t, int32_t> CONTIG_COUNTS;
typedef std::unordered_map<int32_t, CONTIG_COUNTS> ORDERING_COUNTS;

typedef struct ORDERING_INFO
{
    ORDERING_TIG* init_genome;
    ORDERING_COUNTS edges;
    int32_t contig_size;

    double best;
    uint64_t updated;
} ORDERING_INFO;

#endif // ORDERING_TYPE_H
