#ifndef PARTITION_H
#define PARTITION_H

#include "partition_type.hpp"

void partition_create_init_clusters(const HMR_CONTIGS &contigs, CLUSTER_INFO &cluster_info);

CONTIG_ID_SET partition_skip_few_res(const HMR_CONTIGS& contigs, const int32_t min_re);
void partition_skip_empty_links(const MAP_LINK_DENSITY& link_densities, CONTIG_ID_SET &skipped);

void partition_init_factors(int32_t num_of_contigs, LINK_FACTORS& factors_info);
void partition_free_factors(LINK_FACTORS& factors_info);
void partition_factors_size_proc(uint64_t, void*);
void partition_factors_edge_proc(const HMR_EDGE_INFO *edge, void* user);

void partition_link_densities_size_proc(uint64_t edge_size, void* user);
void partition_link_densities_edge_proc(const HMR_EDGE_INFO* edge, void* user);

void partition_remove_skipped_contigs(LINK_FACTORS& contig_factors, CLUSTER_INFO& cluster_info, const char* output_prefix);

void partition_cluster(CLUSTER_INFO &cluster_info, int32_t num_of_groups);

#endif // PARTITION_H
