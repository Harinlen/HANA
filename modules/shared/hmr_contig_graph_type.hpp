#ifndef HMR_CONTIG_GRAPH_TYPE_H
#define HMR_CONTIG_GRAPH_TYPE_H

#include <cstdint>
#include <vector>
#include <list>
#include <unordered_set>

typedef struct HMR_CONTIG
{
    int32_t length;
    int32_t name_size;
    char *name;
} HMR_CONTIG;

typedef std::vector<HMR_CONTIG> HMR_CONTIGS;

typedef union HMR_EDGE
{
    struct {
        int32_t start;
        int32_t end;
    } pos;
    uint64_t data;
} HMR_EDGE;

typedef struct HMR_EDGE_COUNTER
{
    HMR_EDGE edge;
    int32_t pairs;
    int32_t qualified_pairs;
} HMR_EDGE_COUNTER;

typedef std::vector<HMR_EDGE_COUNTER> HMR_EDGE_COUNTERS;

typedef struct HMR_MAPPING
{
    int32_t refID;
    int32_t pos;
    int32_t next_refID;
    int32_t next_pos;
} HMR_MAPPING;

typedef std::list<int32_t> HMR_CONTIG_INVALID_IDS;
typedef std::unordered_set<int32_t> HMR_CONTIG_INVALID_SET;

typedef std::unordered_set<int32_t> CONTIG_ID_SET;
typedef std::vector<int32_t> CONTIG_ID_VECTOR;

typedef std::vector<CONTIG_ID_VECTOR> CONTIG_ID_VECTORS;
typedef std::vector<CONTIG_ID_VECTORS> CHROMOSOME_GROUP;

#endif // HMR_CONTIG_GRAPH_TYPE_H