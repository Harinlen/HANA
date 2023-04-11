#ifndef HMR_CONTIG_GRAPH_TYPE_H
#define HMR_CONTIG_GRAPH_TYPE_H

#include <cstdint>
#include <vector>
#include <list>
#include <unordered_set>
#include <unordered_map>

typedef struct HMR_NODE
{
    int32_t length;
    int32_t enzyme_count;
} HMR_NODE;

typedef std::vector<HMR_NODE> HMR_NODES;

typedef struct HMR_NODE_NAME
{
    int32_t name_size;
    char* name;
} HMR_NODE_NAME;

typedef std::vector<HMR_NODE_NAME> HMR_NODE_NAMES;

typedef struct HMR_CONTIGS
{
    HMR_NODES contigs;
    HMR_NODE_NAMES names;
} HMR_CONTIGS;

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
typedef std::unordered_map<int32_t, HMR_CONTIG_ID_VEC> HMR_ALLELE_TABLE;
typedef std::unordered_map<int32_t, HMR_CONTIG_ID_SET> HMR_ALLELE_MAP;

#endif // HMR_CONTIG_GRAPH_TYPE_H
