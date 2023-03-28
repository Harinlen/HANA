#ifndef ORDERING_TYPE_H
#define ORDERING_TYPE_H

#include "hmr_contig_graph_type.hpp"

typedef std::unordered_map<int32_t, double> CONTIG_WEIGHTS;
typedef std::unordered_map<int32_t, CONTIG_WEIGHTS> ORDERING_WEIGHTS;

typedef std::unordered_map<int32_t, HMR_CONTIG> CONTIG_MAP;

typedef struct ORDERING_INFO
{
    CONTIG_ID_VECTOR contig_group;
    CONTIG_MAP contigs;
    ORDERING_WEIGHTS edges;
    int32_t contig_index;
} ORDERING_INFO;

#endif // ORDERING_TYPE_H
