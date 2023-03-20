#ifndef HMR_ALLELE_TYPE_H
#define HMR_ALLELE_TYPE_H

#include "hmr_contig_graph_type.hpp"

typedef struct
{
    HMR_CONTIG* source;
    std::vector<HMR_CONTIG*> conflicts;
} ALLELE_RECORD;

#endif // HMR_ALLELE_TYPE_H