#ifndef PARTITION_H
#define PARTITION_H

#include "partition_type.hpp"
#include "hmr_contig_graph_type.hpp"

void partition_create_init_clusters(const HMR_CONTIGS &contigs, CLUSTER_INFO &cluster_info);

void partition_init_link_densities(const HMR_CONTIGS &contigs, const HMR_EDGE_MAP &edge_map,
                                   CLUSTER_INFO &cluster_info,
                                   const int32_t min_re, const int32_t max_link_density,
                                   const char *output_prefix);

void partition_init_merge_requests(CLUSTER_INFO &cluster_info);
void partition_cluster(CLUSTER_INFO &cluster_info, int32_t num_of_groups);

#endif // PARTITION_H
