#ifndef ORDERING_LOADER_H
#define ORDERING_LOADER_H

#include <cstdint>

#include "hmr_contig_graph_type.hpp"

void ordering_contig_size_proc(int32_t node_count, void* user);
void ordering_contig_node_proc(int32_t length, int32_t enzyme_count, int32_t name_size, char* name_buff, void* user);

void ordering_edge_map_size_proc(uint64_t edge_sizes, void* user);
void ordering_edge_map_data_proc(const HMR_EDGE_INFO& edge_info, void* user);

#endif // ORDERING_LOADER_H
