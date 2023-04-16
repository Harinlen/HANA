#ifndef DRAFT_MAPPINGS_H
#define DRAFT_MAPPINGS_H

#include "draft_mappings_type.hpp"

typedef struct ALLELE_NEIGHBOUR
{
    int32_t connected_id;
    int32_t count;
} ALLELE_NEIGHBOUR;

void draft_mappings_build_edges(HMR_MAPPING* mapping, int32_t buf_size, void* user);

void draft_mappings_remove_edge(EDGE_COUNT_MAP &edge_map, int32_t node_a, int32_t node_b);
void draft_mappings_build_edge_counter(HMR_MAPPING* mapping, int32_t buf_size, void* user);

#endif // DRAFT_MAPPINGS_H
