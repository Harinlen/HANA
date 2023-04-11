#ifndef PARTITION_H
#define PARTITION_H

#include "partition_type.hpp"

void partition_init_clusters(const HMR_NODES &nodes, const HMR_CONTIG_ID_VEC &invalid_nodes, CLUSTER_INFO &info);
void partition_free_clusters(CLUSTER_INFO& info);

void partition_edge_size_proc(uint64_t edge_size, void* user);
void partition_edge_proc(HMR_EDGE_INFO* edges, int32_t edge_size, void* user);

void partition_cluster(CLUSTER_INFO& cluster_info, int32_t num_of_groups);

double partition_contig_cluster_linkage(int32_t contig_id, HMR_CONTIG_ID_VEC* cluster, const GRAPH_LINK_DENSITY& link_density, bool* has_linkage);
void partition_recover(const std::vector<HMR_CONTIG_ID_VEC*>& clusters, const HMR_CONTIG_ID_VEC& invalid_ids, const int32_t non_info_ratio, CLUSTER_INFO& info);

#endif // PARTITION_H