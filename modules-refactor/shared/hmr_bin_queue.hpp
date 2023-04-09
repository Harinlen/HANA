#ifndef HMR_BIN_QUEUE_H
#define HMR_BIN_QUEUE_H

#include <condition_variable>

typedef struct HMR_BIN_SLICE
{
    char *data;
    size_t data_size;
} HMR_BIN_SLICE;

typedef struct HMR_BIN_QUEUE
{
    bool finish, force_end;
    HMR_BIN_SLICE *slices;
    size_t head, tail, size;
    std::mutex mutex;
    std::condition_variable pop_cv, push_cv;
} HMR_BIN_QUEUE;

typedef struct HMR_BIN_DATA_BUF
{
    char *data;
    size_t size, reserve, offset;
} HMR_BIN_DATA_BUF;

void hmr_bin_queue_create(HMR_BIN_QUEUE **queue, size_t size);
void hmr_bin_queue_free(HMR_BIN_QUEUE *queue);

void hmr_bin_queue_push(HMR_BIN_QUEUE *queue, char *raw_data, size_t raw_data_size);
HMR_BIN_SLICE hmr_bin_queue_pop(HMR_BIN_QUEUE *queue);
void hmr_bin_queue_finish(HMR_BIN_QUEUE* queue);

void hmr_bin_buf_create(HMR_BIN_DATA_BUF **buf);
char *hmr_bin_buf_fetch(HMR_BIN_DATA_BUF *buf, HMR_BIN_QUEUE *queue, size_t size);
inline uint32_t hmr_bin_buf_fetch_uint32(HMR_BIN_DATA_BUF* buf, HMR_BIN_QUEUE* queue)
{
    return *(reinterpret_cast<uint32_t*>(hmr_bin_buf_fetch(buf, queue, 4)));
}
char hmr_bin_buf_getc(HMR_BIN_DATA_BUF* buf, HMR_BIN_QUEUE* queue);
void hmr_bin_buf_free(HMR_BIN_DATA_BUF* buf);

#endif // HMR_BIN_QUEUE_H
