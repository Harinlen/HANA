#include <cassert>
#include <zlib.h>

#include "hmr_bin_file.hpp"
#include "hmr_bin_queue.hpp"
#include "hmr_ui.hpp"
#include "hmr_thread_pool.hpp"
#include "hmr_global.hpp"

#include "hmr_bgzf.hpp"

#define WORK_PER_THREAD (512)

typedef struct BGZF_HEADER
{
    uint8_t ID1;
    uint8_t ID2;
    uint8_t CM;
    uint8_t FLG;
    uint32_t MTIME;
    uint8_t XFL;
    uint8_t OS;
    uint16_t XLEN;
} BGZF_HEADER;

typedef struct BGZF_SUB_HEADER
{
    uint8_t SI1;
    uint8_t SI2;
    uint16_t SLEN;
} BGZF_SUB_HEADER;

typedef struct BGZF_FOOTER
{
    uint32_t CRC32;
    uint32_t ISIZE;
} BGZF_FOOTER;

typedef struct HMR_BGZF_DECOMPRESS
{
    char* cdata;
    uint16_t cdata_size;
    size_t offset;
    size_t raw_size;
} HMR_BGZF_DECOMPRESS;

typedef struct BGZF_UNPACK_PARAM
{
    int id;
    bool finished;
    int32_t max_work;
    HMR_BGZF_DECOMPRESS* pool;
    char* bgzf_raw;
    int32_t *worker_completed;
} BGZF_UNPACK_PARAM;

void hmr_bgzf_decompress(int32_t thread_count, std::mutex &mutex, std::mutex& complete_mutex, std::condition_variable &cv, std::condition_variable &join_cv, BGZF_UNPACK_PARAM *param)
{
    //The thread id and decompress pool would never changed.
    const int id = param->id;
    HMR_BGZF_DECOMPRESS* pool = param->pool;
    while (true)
    {
        //Keep wait.
        std::unique_lock<std::mutex> lock(mutex);
        cv.wait(lock);
        if (param->finished)
        {
            break;
        }
        //Or else, we have to do the work.
        char* bgzf_raw = param->bgzf_raw;
        //Loop and decompress the data.
        for (int i = id * WORK_PER_THREAD, target = hMin(param->max_work, (id + 1) * WORK_PER_THREAD); i < target; ++i)
        {
            z_stream strm;
            strm.zalloc = Z_NULL;
            strm.zfree = Z_NULL;
            strm.opaque = Z_NULL;
            strm.next_in = reinterpret_cast<Bytef*>(pool[i].cdata);
            strm.avail_in = pool[i].cdata_size;
            strm.next_out = reinterpret_cast<Bytef*>(bgzf_raw + pool[i].offset);
            strm.avail_out = static_cast<uInt>(pool[i].raw_size);
            //Add 32 to enable zlib and gzip decoding with header detection.
            if (Z_OK != inflateInit2(&strm, -15))
            {
                time_error(-1, "Failed to initialize decompressor stream.");
            }
            //Decompress the data.
            int error = inflate(&strm, Z_FULL_FLUSH);
            //Recover the compress data memory.
            free(pool[i].cdata);
            //Close the zlib stream.
            inflateEnd(&strm);
        }
        //Okay, mission complete, increase the counter.
        {
            std::unique_lock<std::mutex> counter_lock(complete_mutex);
            ++(*param->worker_completed);
            //Notify the join cv.
            if (*param->worker_completed == thread_count)
            {
                join_cv.notify_all();
            }
        }
    }
}

void hmr_bgzf_parse(FILE* bgzf_file, HMR_BIN_QUEUE* queue, int threads)
{
    //Get the total file size.
    fseek(bgzf_file, 0L, SEEK_END);
#ifdef _MSC_VER
    size_t total_size = _ftelli64(bgzf_file);
#else
    size_t total_size = ftello64(bgzf_file);
#endif
    fseek(bgzf_file, 0L, SEEK_SET);
    //For UI output.
    size_t report_size = (total_size + 9) / 10, report_pos = report_size;
    //Prepare the decompression buffer.
    size_t block_offset = 0;
    BGZF_HEADER header_buf;
    BGZF_FOOTER footer_buf;
    int32_t buf_size = WORK_PER_THREAD * threads, buf_used = 0;
    HMR_BGZF_DECOMPRESS* bgzf_buf = static_cast<HMR_BGZF_DECOMPRESS*>(malloc(sizeof(HMR_BGZF_DECOMPRESS) * buf_size));
    if (!bgzf_buf)
    {
        time_error(-1, "Failed to create BGZF buffer, not enough memory");
    }
    assert(bgzf_buf);
    //Create the working pool.
    int32_t worker_completed = 0;
    std::mutex worker_complete_mutex;
    std::mutex* worker_mutex = new std::mutex[threads], join_mutex;
    if (!worker_mutex)
    {
        time_error(-1, "Not enough memory for creating BGZF mutex.");
    }
    std::condition_variable *worker_cv = new std::condition_variable[threads], join_cv;
    if (!worker_cv)
    {
        time_error(-1, "Not enough memory for creating BGZF condition variable.");
    }
    std::thread* workers = new std::thread[threads];
    if (!worker_cv)
    {
        time_error(-1, "Not enough memory for creating BGZF worker threads.");
    }
    BGZF_UNPACK_PARAM* worker_params = new BGZF_UNPACK_PARAM[threads];
    if (!worker_params)
    {
        time_error(-1, "Not enough memory for creating BGZF worker parameters.");
    }
    for (int32_t i = 0; i < threads; ++i)
    {
        worker_params[i] = BGZF_UNPACK_PARAM { i , false, buf_size, bgzf_buf, NULL, &worker_completed };
        workers[i] = std::thread(hmr_bgzf_decompress, threads, std::ref(worker_mutex[i]), std::ref(worker_complete_mutex), std::ref(worker_cv[i]), std::ref(join_cv), & worker_params[i]);
    }
    //Read while to the end of the file.
    while (!queue->finish && fread(&header_buf, sizeof(BGZF_HEADER), 1, bgzf_file) > 0)
    {
        //Read the Xlen data.
        char* subfield_data = static_cast<char*>(malloc(header_buf.XLEN));
        if (!subfield_data)
        {
            time_error(-1, "Not enough memory for sub field data buffer.");
        }
        //Read the data.
        fread(subfield_data, header_buf.XLEN, 1, bgzf_file);
        //Go through the header.
        uint16_t subfield_left = header_buf.XLEN, bsize = 0;
        char* subfield_pos = subfield_data;
        while (subfield_left > 0)
        {
            BGZF_SUB_HEADER* subfield = reinterpret_cast<BGZF_SUB_HEADER*>(subfield_pos);
            //Check the ID matches the bsize.
            if (subfield->SI1 == 66 && subfield->SI2 == 67 && subfield->SLEN == 2)
            {
                bsize = *(reinterpret_cast<uint16_t*>(subfield_pos + sizeof(BGZF_SUB_HEADER)));
            }
            subfield_pos += subfield->SLEN + sizeof(BGZF_SUB_HEADER);
            subfield_left -= subfield->SLEN + sizeof(BGZF_SUB_HEADER);
        }
        free(subfield_data);
        uint16_t cdata_size = bsize - header_buf.XLEN - 19;
        char* cdata = static_cast<char*>(malloc(cdata_size));
        //Reading the compressed data.
        assert(cdata);
        fread(cdata, cdata_size, 1, bgzf_file);
        //Fetch the footer data.
        fread(&footer_buf, sizeof(BGZF_FOOTER), 1, bgzf_file);
        //Buffer until reach the decompress limit.
        bgzf_buf[buf_used] = HMR_BGZF_DECOMPRESS{ cdata, cdata_size, block_offset, footer_buf.ISIZE };
        ++buf_used;
        //Increase the block offset, and keep going.
        block_offset += footer_buf.ISIZE;
        //Check whether we are reaching the decompression limitation.
        if (buf_used == buf_size)
        {
            //Create the pool.
            char* bgzf_raw = static_cast<char*>(malloc(block_offset));
            //Decompress the data.
            worker_completed = 0;
            for (int i = 0; i < threads; ++i)
            {
                //Update the bgzf raw for the parameter.
                worker_params[i].bgzf_raw = bgzf_raw;
                worker_cv[i].notify_one();
            }
            //Wait for all the threads complete.
            {
                std::unique_lock<std::mutex> join_lock(join_mutex);
                join_cv.wait(join_lock);
            }
            //Push the data to the parsing queue.
            hmr_bin_queue_push(queue, bgzf_raw, block_offset);
            //Reset the buffer used.
            buf_used = 0;
            block_offset = 0;
        }
        //Check should we report the position.
#ifdef _MSC_VER
        size_t bgzf_pos = _ftelli64(bgzf_file);
#else
        size_t bgzf_pos = ftello64(bgzf_file);
#endif
        if (bgzf_pos >= report_pos)
        {
            float percent = static_cast<float>(bgzf_pos) / static_cast<float>(total_size) * 100.0f;
            time_print("BGZF parsed %.1f%%", percent);
            report_pos += report_size;
        }
    }
    //Check whether we still have data left.
    if (buf_used > 0)
    {
        //Create the pool.
        char* bgzf_raw = static_cast<char*>(malloc(block_offset));
        //Decompress the data.
        worker_completed = 0;
        for (int i = 0; i < threads; ++i)
        {
            //Update the bgzf raw for the parameter.
            worker_params[i].max_work = buf_used;
            worker_params[i].bgzf_raw = bgzf_raw;
            worker_cv[i].notify_one();
        }
        //Wait for all the threads complete.
        {
            std::unique_lock<std::mutex> join_lock(join_mutex);
            join_cv.wait(join_lock, [&] { return worker_completed == threads; });
        }
        //Push the data to the parsing queue.
        hmr_bin_queue_push(queue, bgzf_raw, block_offset);
    }
    //Let the workers complete their jobs.
    for (int i = 0; i < threads; ++i)
    {
        worker_params[i].finished = true;
        worker_cv[i].notify_one();
        workers[i].join();
    }
    delete[] workers;
    //Mark BGZF parsing complete.
    hmr_bin_queue_finish(queue);
}

HMR_BGZF_HANDLER* hmr_bgzf_open(const char* filepath, int threads)
{
    //Read the BGZF file.
    HMR_BGZF_HANDLER* bgzf_handler = new HMR_BGZF_HANDLER();
    FILE* bgzf_file = NULL;
    if (!bin_open(filepath, &bgzf_file, "rb"))
    {
        time_error(-1, "Failed to read BGZF file %s\n", filepath);
    }
    bgzf_handler->bgzf_file = bgzf_file;
    //Allocate the processing queue, 3 for triple buffer.
    hmr_bin_queue_create(&(bgzf_handler->queue), 3);
    //Prepare the buffer.
    hmr_bin_buf_create(&bgzf_handler->buffer);
    //Start the BGZF parsing thread.
    bgzf_handler->parse_thread = std::thread(hmr_bgzf_parse, bgzf_file, bgzf_handler->queue, threads);
    //Provide the GZIP handler.
    return bgzf_handler;
}

void hmr_bgzf_close(HMR_BGZF_HANDLER* bgzf_handler)
{
    //Check whether the queue is marked as finished.
    if (!bgzf_handler->queue->finish)
    {
        hmr_bin_queue_finish(bgzf_handler->queue);
    }
    //Wait for parse thread to complete.
    bgzf_handler->parse_thread.join();
    //Free the queue and buffer.
    hmr_bin_buf_free(bgzf_handler->buffer);
    hmr_bin_queue_free(bgzf_handler->queue);
    //Close the file.
    fclose(bgzf_handler->bgzf_file);
}
