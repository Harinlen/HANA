#ifndef EXTRACT_MAPPING_TYPE_H
#define EXTRACT_MAPPING_TYPE_H

#include <mutex>
#include <condition_variable>

#include "hmr_contig_graph_type.hpp"

#include "extract_index_map.hpp"

typedef struct MAPPING_BUFFER
{
    HMR_MAPPING* buffer;
    size_t buffer_size, buffer_offset;
} MAPPING_BUFFER;

typedef struct MAPPING_WORKER
{
    MAPPING_BUFFER valid_buffer;
    int32_t start_pos, end_pos;
    int32_t counter;
} MAPPING_WORKER;

typedef struct MAPPING_WORKER_SYNC
{
    std::mutex *start_mutex, complete_mutex, file_mutex;
    std::condition_variable *start_cv, complete_cv;
    bool *start_signal = NULL;
    int32_t completed_worker = 0, total_worker = 0;
    bool exit = false;
    FILE* reads_file;
} MAPPING_WORKER_SYNC;

#endif // EXTRACT_MAPPING_TYPE_H
