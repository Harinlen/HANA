#ifndef ORDERING_LOADER_H
#define ORDERING_LOADER_H

#include <cstdint>

#include "ordering_type.hpp"

typedef struct ORDERING_EDGE_LOADER
{
    const HMR_CONTIG_ID_VEC& contig_group;
    ORDERING_COUNTS& edges;
} ORDERING_EDGE_LOADER;

void ordering_edge_map_size_proc(uint64_t edge_size, void* user);
void ordering_edge_map_data_proc(HMR_EDGE_INFO* edges, int32_t edge_size, void* user);

#endif // ORDERING_LOADER_H
