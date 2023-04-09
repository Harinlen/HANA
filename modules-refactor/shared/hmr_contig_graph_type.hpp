#ifndef HMR_CONTIG_GRAPH_TYPE_H
#define HMR_CONTIG_GRAPH_TYPE_H

#include <cstdint>
#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>

typedef struct HMR_CONTIG
{
    int32_t length;
    int32_t enzyme_count;
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

typedef struct HMR_EDGE_INFO
{
    int32_t start;
    int32_t end;
    uint64_t pairs;
    double weights;
} HMR_EDGE_INFO;

typedef std::vector<HMR_EDGE_INFO> HMR_EDGE_COUNTERS;

typedef struct HMR_MAPPING
{
    int32_t refID;
    int32_t pos;
    int32_t next_refID;
    int32_t next_pos;
} HMR_MAPPING;

typedef std::vector<int32_t> HMR_CONTIG_ID_VEC;
typedef std::unordered_set<int32_t> HMR_CONTIG_ID_SET;
typedef std::list<int32_t> HMR_CONTIG_ID_CHAIN;
typedef std::unordered_map<int32_t, HMR_CONTIG_ID_SET> HMR_ALLELE_TABLE;

#endif // HMR_CONTIG_GRAPH_TYPE_H
