#ifndef DRAFT_MAPPINGS_TYPE_H
#define DRAFT_MAPPINGS_TYPE_H

#include "hmr_contig_graph_type.hpp"

typedef std::unordered_map<uint64_t, int32_t> EDGE_COUNTER;

typedef std::unordered_map<int32_t, int32_t> NODE_COUNT_MAP;
typedef std::vector<NODE_COUNT_MAP> EDGE_COUNT_MAP;

#endif // DRAFT_MAPPINGS_TYPE_H
