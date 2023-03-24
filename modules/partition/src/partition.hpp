#ifndef PARTITION_H
#define PARTITION_H

#include "hmr_contig_graph_type.hpp"

CONTIG_ID_VECTORS partition_contigs(const HMR_CONTIGS &contigs, const HMR_EDGE_COUNTERS &edge_weights, int num_of_groups);

#endif // PARTITION_H
