#include <cassert>
#include <cstring>

#include "hmr_ui.hpp"

#include "hmr_bin_queue.hpp"

void hmr_bin_queue_create(HMR_BIN_QUEUE **queue, size_t size)
{
    //Allocate the queue.
    (*queue) = new HMR_BIN_QUEUE();
    if(!(*queue))
    {
        time_error(-1, "Failed to create binary buffer queue, no enough memory.");
    }
    //Create the slices buffer.
    (*queue)->slices = static_cast<HMR_BIN_SLICE *>(malloc(sizeof(HMR_BIN_SLICE) * size));
    if (!((*queue)->slices))
    {
        time_error(-1, "Failed to create binary buffer queue slices, no enough memory.");
    }
    (*queue)->size = size;
    (*queue)->head = 0;
    (*queue)->tail = 0;
    (*queue)->finish = false;
}

void hmr_bin_queue_free(HMR_BIN_QUEUE *queue)
{
    //Clear the queue slices.
    free(queue->slices);
    //Clear the queue.
    delete queue;
}

void hmr_bin_queue_push(HMR_BIN_QUEUE *queue, char *raw_data, size_t raw_data_size)
{
    //Check whether the queue is full.
    std::unique_lock<std::mutex> push_lock(queue->mutex);
    queue->push_cv.wait(push_lock, [queue]
    {
        return !((queue->tail+1==queue->head) || (queue->head==0&&queue->tail==queue->size - 1));
    });
    //Push the data to the queue.
    queue->slices[queue->tail] = HMR_BIN_SLICE {raw_data, raw_data_size};
    queue->tail = (queue->tail+1 == queue->size) ? 0 : (queue->tail+1);
    //Notify the conditional variable.
    queue->pop_cv.notify_one();
}

HMR_BIN_SLICE hmr_bin_queue_pop(HMR_BIN_QUEUE *queue)
{
    //Wait until the queue is not empty.
    std::unique_lock<std::mutex> pop_lock(queue->mutex);
    queue->pop_cv.wait(pop_lock, [queue]
    {
        return (queue->head != queue->tail) || queue->finish;
    });
    //Extract the data.
    HMR_BIN_SLICE slice = queue->slices[queue->head];
    if(queue->finish && queue->head == queue->tail)
    {
        slice.data = NULL;
        slice.data_size = 0;
    }
    else
    {
        queue->head = (queue->head+1 == queue->size) ? 0 : (queue->head+1);
    }
    //Notify the push variable.
    queue->push_cv.notify_one();
    return slice;
}

void hmr_bin_queue_finish(HMR_BIN_QUEUE* queue)
{
    std::unique_lock<std::mutex> finish_lock(queue->mutex);
    //Mark queue is finished using.
    queue->finish = true;
    queue->pop_cv.notify_one();
}

void hmr_bin_buf_create(HMR_BIN_DATA_BUF** buf)
{
    //Prepare the buffer.
    HMR_BIN_DATA_BUF* buffer = static_cast<HMR_BIN_DATA_BUF*>(malloc(sizeof(HMR_BIN_DATA_BUF)));
    if (!buffer)
    {
        time_error(-1, "Failed to create binary buffer structure, no enough memory.");
    }
    //Initial the structure.
    buffer->data = NULL;
    buffer->offset = 0;
    buffer->reserve = 0;
    buffer->size = 0;
    *buf = buffer;
}

char *hmr_bin_buf_fetch(HMR_BIN_DATA_BUF *buf, HMR_BIN_QUEUE *queue, size_t size)
{
    //Check the left data is enough.
    size_t residual = buf->size - buf->offset;
    if(residual >= size)
    {
        //Enough to hold the data.
        char *data = buf->data + buf->offset;
        buf->offset += size;
        return data;
    }
    //Well, we have to fill the buffer with the queue data size.
    while(residual < size)
    {
        //Fetch the data from the queue.
        HMR_BIN_SLICE bin_slice = hmr_bin_queue_pop(queue);
        if(bin_slice.data == NULL)
        {
            //Reach the end of the data, still cannot fulfill, failed to fetch.
            return NULL;
        }
        //Check whether we can fast swap to the slice memory.
        if(residual == 0)
        {
            //Directly use the slice memory.
            if(buf->data)
            {
                free(buf->data);
            }
            //Update the buffer.
            buf->data = bin_slice.data;
            buf->size = bin_slice.data_size;
            buf->reserve = bin_slice.data_size;
        }
        else
        {
            //Copy the residual data.
            assert(buf->data + buf->offset);
            memcpy(buf->data, buf->data + buf->offset, residual);
            //Calculate the new buffer size.
            size_t data_size = residual + bin_slice.data_size;
            //Increase the buffer size when necessary.
            if(data_size > buf->reserve)
            {
                char *expected = static_cast<char *>(realloc(buf->data, data_size));
                if (!expected)
                {
                    time_error(-1, "Failed to expand the binary buffer, no enough memory.");
                }
                assert(expected);
                buf->data = expected;
                buf->reserve = data_size;
            }
            //Copy the slice data to buffer.
            memcpy(buf->data + residual, bin_slice.data, bin_slice.data_size);
            //Update the data size.
            buf->size = data_size;
            //Release the slice data.
            free(bin_slice.data);
        }
        //Update the residual data.
        buf->offset = 0;
        residual = buf->size;
    }
    //Now the data should be enough.
    char *data = buf->data;
    buf->offset = size;
    return data;
}

char hmr_bin_buf_getc(HMR_BIN_DATA_BUF* buf, HMR_BIN_QUEUE* queue)
{
    //Fetch 1 bytes from the queue.
    char* result = hmr_bin_buf_fetch(buf, queue, 1);
    if (NULL == result)
    {
        return -1;
    }
    //Or else, parse the current char.
    return result[0];
}

void hmr_bin_buf_free(HMR_BIN_DATA_BUF* buf)
{
    if (buf->data)
    {
        free(buf->data);
    }
    free(buf);
}
