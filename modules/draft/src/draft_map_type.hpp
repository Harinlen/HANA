#ifndef DRAFT_MAP_TYPE
#define DRAFT_MAP_TYPE

#include <cstdlib>
#include <cstdio>
#include <unordered_map>
#include <string>
#include <mutex>
#include <thread>

#include "hmr_contig_graph_type.hpp"

#include "draft_fasta_type.hpp"

typedef std::unordered_map<std::string, int32_t> CONTIG_ID_MAP;
typedef std::unordered_map<uint64_t, int32_t> RAW_EDGE_MAP;

typedef struct MAPPING_CONTIG_MAP
{
    int32_t* contig_id_map;
    int32_t contig_idx;
} MAPPING_CONTIG_MAP;

//Mapping thread buffer.
typedef struct MAPPING_FILTER_USER
{
    ENZYME_RANGES* contig_ranges = NULL;
    uint8_t mapq = 0;
    int32_t threads = 0;

    MAPPING_INFO *buf_0 = NULL, *buf_1 = NULL, *proc_buf = NULL, *build_buf = NULL;
    size_t buf_0_size = 0, buf_1_size = 0, buf_max_size = 0, *proc_size = NULL, *build_size = NULL;
    bool activated_0 = true;
    FILE* reads_file = NULL;
    std::mutex reads_mutex, finished_mutex, counter_mutex;
    std::condition_variable finished_cv;
    int32_t finished_counter = 0;
    bool is_working = false;
    bool exit = false;
} MAPPING_FILTER_USER;

typedef struct MAPPING_FILTER_WORKER
{
    RAW_EDGE_MAP edges;
    bool start;
    std::mutex start_mutex;
    std::condition_variable start_cv;
    size_t buf_start, buf_end;
    char* mapping_buf;
    size_t mapping_buf_size, mapping_buf_offset;
} MAPPING_FILTER_WORKER;

//Read position record.
typedef struct MAPPING_DRAFT_USER
{
    const CONTIG_ID_MAP& contig_ids;
    const HMR_CONTIG_INVALID_SET& invalid_ids;

    MAPPING_FILTER_USER& filter;
    MAPPING_FILTER_WORKER* workers;
    MAPPING_CONTIG_MAP& map;
} MAPPING_DRAFT_USER;

#endif // DRAFT_MAP_TYPE
