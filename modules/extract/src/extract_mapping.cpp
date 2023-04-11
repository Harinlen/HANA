#include <cassert>
#include <cstdlib>
#include <thread>
#include <mutex>

#include "hmr_global.hpp"
#include "hmr_bam.hpp"
#include "hmr_path.hpp"
#include "hmr_ui.hpp"

#include "extract_fasta_type.hpp"

#include "extract_mapping.hpp"

bool position_in_range(int32_t pos, const CONTIG_ENZYME_RANGE& ranges)
{
    for (size_t i = 0; i < ranges.size(); ++i)
    {
        const ENZYME_RANGE& r = ranges[i];
        //Check whether the position is in the range.
        if (r.start <= pos && pos <= r.end)
        {
            return true;
        }
    }
    //No position matched.
    return false;
}

void mapping_buffer_init(MAPPING_BUFFER& buf, size_t buffer_size)
{
    buf.buffer = static_cast<HMR_MAPPING*>(malloc(sizeof(HMR_MAPPING) * buffer_size));
    if (!buf.buffer)
    {
        time_error(-1, "Failed to allocate memory for mapping info buffer.");
    }
    assert(buf.buffer);
    buf.buffer_offset = 0;
    buf.buffer_size = buffer_size;
}

inline bool mapping_buffer_is_full(MAPPING_BUFFER& buf)
{
    return buf.buffer_offset == buf.buffer_size;
}

inline void mapping_buffer_push(MAPPING_BUFFER& buf, const HMR_MAPPING &item)
{
    buf.buffer[buf.buffer_offset++] = item;
}

inline void mapping_buffer_dump(MAPPING_BUFFER& buf, FILE* reads_file)
{
    if (buf.buffer_offset == 0)
    {
        return;
    }
    fwrite(buf.buffer, sizeof(HMR_MAPPING), buf.buffer_offset, reads_file);
    buf.buffer_offset = 0;
}

inline void mapping_buffer_free(MAPPING_BUFFER& buf)
{
    free(buf.buffer);
}

inline void worker_init(MAPPING_WORKER& worker, int32_t worker_index, int32_t buffer_size)
{
    mapping_buffer_init(worker.valid_buffer, buffer_size);
    worker.start_pos = worker_index * buffer_size;
    worker.end_pos = worker.start_pos + buffer_size;
}

inline void worker_free(MAPPING_WORKER& worker)
{
    mapping_buffer_free(worker.valid_buffer);
}

inline void mapping_worker_sync_init(MAPPING_WORKER_SYNC& worker_sync, FILE *reads_file, int32_t num_of_worker)
{
    worker_sync.reads_file = reads_file;
    worker_sync.completed_worker = num_of_worker;
    worker_sync.exit = false;
    worker_sync.total_worker = num_of_worker;
    worker_sync.start_signal = new bool[num_of_worker];
    for (int32_t i = 0; i < num_of_worker; ++i)
    {
        worker_sync.start_signal[i] = false;
    }
    worker_sync.start_mutex = new std::mutex[num_of_worker];
    worker_sync.start_cv = new std::condition_variable[num_of_worker];
}

inline void mapping_worker_sync_start(MAPPING_WORKER_SYNC& worker_sync)
{
    worker_sync.completed_worker = 0;
    for (int32_t i = 0; i < worker_sync.total_worker; ++i)
    {
        worker_sync.start_signal[i] = true;
        worker_sync.start_cv[i].notify_one();
    }
}

inline void mapping_worker_sync_dump(MAPPING_WORKER_SYNC& worker_sync, MAPPING_BUFFER& buf)
{
    worker_sync.file_mutex.lock();
    mapping_buffer_dump(buf, worker_sync.reads_file);
    worker_sync.file_mutex.unlock();
}

inline void mapping_worker_sync_wait_ready(MAPPING_WORKER_SYNC& worker_sync)
{
    if (worker_sync.completed_worker != worker_sync.total_worker)
    {
        std::unique_lock<std::mutex> lock(worker_sync.complete_mutex);
        worker_sync.complete_cv.wait(lock, [&] {return worker_sync.completed_worker == worker_sync.total_worker; });
    }
}

inline void mapping_worker_sync_complete_one(MAPPING_WORKER_SYNC& worker_sync)
{
    worker_sync.complete_mutex.lock();
    ++worker_sync.completed_worker;
    if (worker_sync.completed_worker == worker_sync.total_worker)
    {
        worker_sync.complete_cv.notify_one();
    }
    worker_sync.complete_mutex.unlock();
}

inline void mapping_worker_sync_exit(MAPPING_WORKER_SYNC& worker_sync, std::thread* workers)
{
    //Wait for all the workers are complete their work.
    mapping_worker_sync_wait_ready(worker_sync);
    //Set the exit signal.
    worker_sync.exit = true;
    for (int32_t i = 0; i < worker_sync.total_worker; ++i)
    {
        //Set the start signal.
        worker_sync.start_signal[i] = true;
        //Let the worker start.
        worker_sync.start_cv[i].notify_one();
        //Wait for the worker done.
        workers[i].join();
    }
}

inline void mapping_worker_sync_free(MAPPING_WORKER_SYNC& worker_sync)
{
    delete[] worker_sync.start_signal;
}

// ------ BAM Workers ------

typedef struct BAM_MAPPING_INFO
{
    int32_t refID;
    int32_t pos;
    int32_t next_refID;
    int32_t next_pos;
    uint16_t flag;
    uint8_t mapq;
} BAM_MAPPING_INFO;

typedef struct BAM_MAPPING_BUFFER
{
    BAM_MAPPING_INFO* buffer;
    size_t buffer_size, buffer_offset;
} BAM_MAPPING_BUFFER;

void bam_mapping_buffer_init(BAM_MAPPING_BUFFER& buf, size_t buffer_size)
{
    buf.buffer = static_cast<BAM_MAPPING_INFO*>(malloc(sizeof(BAM_MAPPING_INFO) * buffer_size));
    if (!buf.buffer)
    {
        time_error(-1, "Failed to allocate memory for mapping info buffer.");
    }
    assert(buf.buffer);
    buf.buffer_offset = 0;
    buf.buffer_size = buffer_size;
}

void bam_mapping_buffer_free(BAM_MAPPING_BUFFER& buf)
{
    free(buf.buffer);
}

inline bool bam_mapping_buffer_is_full(BAM_MAPPING_BUFFER* buffer)
{
    return buffer->buffer_offset == buffer->buffer_size;
}

typedef struct BAM_EXTRACTOR
{
    CONTIG_ENZYME_RANGES *contig_enzyme_ranges;
    BAM_MAPPING_BUFFER buf_0, buf_1;
    BAM_MAPPING_BUFFER* buf_filling = NULL, * buf_filtering = NULL;
    MAPPING_WORKER_SYNC sync;
    CONTIG_INDEX_MAP *contig_index_map = NULL;
    int32_t* bam_id_map = NULL;
    int32_t bam_contig_id = 0;
    uint8_t mapq = 0;
    uint16_t check_flag = 0;
} BAM_EXTRACTOR;

void bam_extractor_init(BAM_EXTRACTOR& extractor, CONTIG_INDEX_MAP *contig_index_map, CONTIG_ENZYME_RANGES* contig_enzyme_ranges, uint8_t mapq, FILE *reads_file, int32_t thread_buffer_size, int32_t num_of_worker, uint16_t check_flag)
{
    size_t buffer_size = static_cast<size_t>(thread_buffer_size * num_of_worker);
    mapping_worker_sync_init(extractor.sync, reads_file, num_of_worker);
    bam_mapping_buffer_init(extractor.buf_0, buffer_size);
    bam_mapping_buffer_init(extractor.buf_1, buffer_size);
    extractor.contig_enzyme_ranges = contig_enzyme_ranges;
    extractor.buf_filling = &extractor.buf_0;
    extractor.buf_filtering = &extractor.buf_1;
    extractor.contig_index_map = contig_index_map;
    extractor.bam_id_map = NULL;
    extractor.bam_contig_id = 0;
    extractor.mapq = mapq;
    extractor.check_flag = check_flag;
}

inline int32_t bam_extractor_get_contig_id(const BAM_EXTRACTOR& extractor, int32_t ref_id)
{
    return ref_id < extractor.bam_contig_id ? extractor.bam_id_map[ref_id] : -1;
}

void bam_extractor_free(BAM_EXTRACTOR& extractor)
{
    bam_mapping_buffer_free(extractor.buf_0);
    bam_mapping_buffer_free(extractor.buf_1);
    mapping_worker_sync_free(extractor.sync);
}

void extract_mapping_bam_worker(int32_t id, MAPPING_WORKER &worker, BAM_EXTRACTOR &extractor)
{
    MAPPING_WORKER_SYNC& sync = extractor.sync;
    std::unique_lock<std::mutex> lock(sync.start_mutex[id]);
    while (!sync.exit)
    {
        //Wait for the start signal.
        sync.start_cv[id].wait(lock, [&] {return sync.start_signal[id]; });
        if (sync.exit)
        {
            return;
        }
        //Loop in the worker area.
        BAM_MAPPING_BUFFER* filtering_buf = extractor.buf_filtering;
        for (int32_t i = worker.start_pos; i < worker.end_pos; ++i)
        {
            BAM_MAPPING_INFO& mapping_info = filtering_buf->buffer[i];
            //Check whether the mapping info is valid, then check whether the position is in range.
            int32_t ref_index = bam_extractor_get_contig_id(extractor, mapping_info.refID),
                next_ref_index = bam_extractor_get_contig_id(extractor, mapping_info.next_refID);
            if (ref_index == -1 || next_ref_index == -1 //Reference index cannot be find in the mapping index detection.
                || mapping_info.mapq == 0 || mapping_info.mapq == 255 //Map quality is invalid.
                || mapping_info.mapq < extractor.mapq //Check whether the mapping reaches the minimum quality
                || ((extractor.check_flag & CHECK_FLAG_FLAG) && (mapping_info.flag & 3852)) // Filtered flag from AllHiC.
                || (ref_index == next_ref_index) // We don't care about the pairs on the same contigs.
                || ((extractor.check_flag & CHECK_FLAG_RANGE) && (!position_in_range(mapping_info.pos, (*extractor.contig_enzyme_ranges)[ref_index])))) // Or the position is not in the position.
            {
                continue;
            }
            //Save the mapping info to the buffer.
            if (mapping_buffer_is_full(worker.valid_buffer))
            {
                mapping_worker_sync_dump(extractor.sync, worker.valid_buffer);
            }
            mapping_buffer_push(worker.valid_buffer, HMR_MAPPING{ mapping_info.refID, mapping_info.pos, mapping_info.next_refID, mapping_info.next_pos });
        }
        //Reset the start signal.
        sync.start_signal[id] = false;
        //Increase the sync counter.
        mapping_worker_sync_complete_one(extractor.sync);
    }
}

void extract_bam_num_of_contigs(uint32_t num_of_contigs, void* user)
{
    BAM_EXTRACTOR* bam_extractor = static_cast<BAM_EXTRACTOR*>(user);
    //Initial the BAM contig map.
    bam_extractor->bam_id_map = new int32_t[num_of_contigs];
    if (!bam_extractor->bam_id_map)
    {
        time_error(-1, "No enough memory for BAM contig id map.");
    }
    //Initial all the index mapping to -1.
    for (uint32_t i = 0; i < num_of_contigs; ++i)
    {
        bam_extractor->bam_id_map[i] = -1;
    }
}

void extract_bam_contig(uint32_t name_length, char* name, uint32_t, void* user)
{
    BAM_EXTRACTOR* bam_extractor = static_cast<BAM_EXTRACTOR*>(user);
    //Set the contig.
    bam_extractor->bam_id_map[bam_extractor->bam_contig_id] = extract_contig_index_get(*bam_extractor->contig_index_map, std::string(name, name_length));
    ++bam_extractor->bam_contig_id;
}

void extract_bam_read_align(size_t block_id, const BAM_BLOCK_HEADER* bam_block, void* user)
{
    BAM_EXTRACTOR* bam_extractor = static_cast<BAM_EXTRACTOR*>(user);
    //Fill the data to build extractor.
    if (bam_mapping_buffer_is_full(bam_extractor->buf_filling))
    {
        //Wait for all the workers are ready.
        mapping_worker_sync_wait_ready(bam_extractor->sync);
        //Swap the filling and processing buffer.
        BAM_MAPPING_BUFFER* temp = bam_extractor->buf_filling;
        bam_extractor->buf_filling = bam_extractor->buf_filtering;
        bam_extractor->buf_filtering = temp;
        //Start to process the filling buffer.
        mapping_worker_sync_start(bam_extractor->sync);
        //Reset the filling offset.
        bam_extractor->buf_filling->buffer_offset = 0;
    }
    //Save the data to buffer filling.
    BAM_MAPPING_BUFFER* filling = bam_extractor->buf_filling;
    filling->buffer[filling->buffer_offset] = BAM_MAPPING_INFO
    {
        bam_block->refID,
        bam_block->pos,
        bam_block->next_refID,
        bam_block->next_pos,
        bam_block->flag,
        bam_block->mapq
    };
    ++filling->buffer_offset;
}

// ------ BAM Workers End ------

void extract_mapping_file(const char* filepath, CONTIG_INDEX_MAP* index_map, FILE* reads_file, CONTIG_ENZYME_RANGES* contig_enzyme_ranges, uint16_t check_flag, uint8_t mapq, int32_t thread_buffer_size, int32_t threads)
{
    if (path_ends_with(filepath, ".bam"))
    {
        int32_t num_of_worker = (threads + 1) >> 1;
        //Initialize the extractor.
        BAM_EXTRACTOR bam_extractor;
        bam_extractor_init(bam_extractor, index_map, contig_enzyme_ranges, mapq, reads_file, thread_buffer_size, num_of_worker, check_flag);
        //Prepare the worker buffer.
        MAPPING_WORKER* worker_buffer = new MAPPING_WORKER[num_of_worker];
        //Start BAM filter workers.
        std::thread *workers = new std::thread[num_of_worker];
        for (int32_t i = 0; i < num_of_worker; ++i)
        {
            //Initialize the worker.
            worker_init(worker_buffer[i], i, thread_buffer_size);
            //Start the worker.
            workers[i] = std::thread(extract_mapping_bam_worker, i, std::ref(worker_buffer[i]), std::ref(bam_extractor));
        }
        //Start parsing the bam file.
        hmr_bam_read(filepath, BAM_MAPPING_PROC{ extract_bam_num_of_contigs ,extract_bam_contig, extract_bam_read_align }, &bam_extractor, num_of_worker);
        //Wait for all the workers complete.
        mapping_worker_sync_wait_ready(bam_extractor.sync);
        //Check is there any other data left in the last bam worker.
        BAM_MAPPING_BUFFER* filling = bam_extractor.buf_filling;
        if (filling->buffer_offset)
        {
            //We need to process the last part of work.
            int32_t work_per_workers = (static_cast<int32_t>(filling->buffer_offset) + num_of_worker - 1) / num_of_worker, 
                work_start = 0;
            for (int32_t i = 0; i < num_of_worker; ++i)
            {
                worker_buffer[i].start_pos = work_start;
                work_start += work_per_workers;
                worker_buffer[i].end_pos = hMin(work_start, static_cast<int32_t>(filling->buffer_offset));
            }
            //Set the filtering worker as current worker.
            bam_extractor.buf_filtering = bam_extractor.buf_filling;
            //Start all the workers.
            mapping_worker_sync_start(bam_extractor.sync);
        }
        //Exit all the workers.
        mapping_worker_sync_exit(bam_extractor.sync, workers);
        //Recover the memory and dump the data left in their buffer.
        for (int32_t i = 0; i < num_of_worker; ++i)
        {
            mapping_buffer_dump(worker_buffer[i].valid_buffer, reads_file);
            worker_free(worker_buffer[i]);
        }
        bam_extractor_free(bam_extractor);
        delete[] workers;
        delete[] worker_buffer;
    }
}
