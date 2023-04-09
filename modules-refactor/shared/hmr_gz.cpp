#include <cstdlib>
#include <cassert>
#include <cstring>
#include <zlib.h>

#include "hmr_ui.hpp"
#include "hmr_global.hpp"
#include "hmr_bin_queue.hpp"

#include "hmr_gz.hpp"

// 4MB Data chunk
#define DATA_CHUNK (4194304)

typedef struct GZIP_HEADER
{
    char magic[2];
    uint8_t compress_method;
    uint8_t flag;
    char mod_time[4];
    uint8_t extra_flags;
    uint8_t os;
} GZIP_HEADER;

void hmr_gzip_parse(FILE *gz_file, HMR_BIN_QUEUE *queue)
{
    //Prepare the zstream.
    z_stream strm;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    //Add 32 to enable zlib and gzip decoding with header detection.
    if (Z_OK != inflateInit2(&strm, 32 + MAX_WBITS))
    {
        time_error(-1, "Failed to initialize decompressor stream.");
    }
    //Get the file size.
    fseek(gz_file, 0L, SEEK_END);
#ifdef _MSC_VER
    size_t total_size = _ftelli64(gz_file);
#else
    size_t total_size = ftello64(gz_file);
#endif
    fseek(gz_file, 0L, SEEK_SET);
    //Allocate the memory for binary reading.
    char *gz_buffer = static_cast<char *>(malloc(DATA_CHUNK));
    if (NULL == gz_buffer)
    {
        time_error(-1, "Failed to allocate GZIP reading buffer.");
    }
    //Read the file to buffer.
    assert(NULL != gz_buffer);
    size_t gz_file_offset = 0;
    assert(gz_buffer != NULL);
    size_t gz_buffer_size = gz_file_offset = fread(gz_buffer, 1, DATA_CHUNK, gz_file);
    //Configure the stream.
    strm.next_in = reinterpret_cast<Bytef *>(gz_buffer);
    strm.avail_in = static_cast<uInt>(gz_buffer_size);
    int error = Z_OK;
    while (Z_STREAM_END != error && !queue->finish)
    {
        //Assume all the compression ratio to be 2x.
        size_t slice_reserved = DATA_CHUNK << 2;
        char* slice_data = static_cast<char*>(malloc(slice_reserved));
        strm.next_out = reinterpret_cast<Bytef*>(slice_data);
        strm.avail_out = static_cast<uInt>(slice_reserved);
        //Decompress the stream.
        error = inflate(&strm, Z_NO_FLUSH);
        //While the data is valid, push the output data to data slice.
        if (Z_OK == error || Z_STREAM_END == error)
        {
            //Send the data.
            hmr_bin_queue_push(queue, slice_data, slice_reserved - strm.avail_out);
        }
        else
        {
            time_error(-1, "Error happens when reading GZIP file.");
        }
        //Check whether we reach the end of the buffer.
        if (strm.avail_in > 0)
        {
            //Move the data to the beginning of the header.
            size_t gz_buffer_used = DATA_CHUNK - strm.avail_in;
            memcpy(gz_buffer, gz_buffer + gz_buffer_used, strm.avail_in);
        }
        //Reset the next input at the start of the buffer.
        strm.next_in = reinterpret_cast<Bytef*>(gz_buffer);
        //Check whether we can still fill the buffer.
        if (gz_file_offset < total_size)
        {
            size_t bytes_expected = hMin(total_size - gz_file_offset, static_cast<size_t>(DATA_CHUNK - strm.avail_in));
            //Fill the buffer.
            fread(gz_buffer + strm.avail_in, 1, bytes_expected, gz_file);
            strm.avail_in += static_cast<uInt>(bytes_expected);
        }
    }
    //Mark GZIP parsing complete.
    hmr_bin_queue_finish(queue);
    // The parsing is completed, free the buffer.
    free(gz_buffer);
    //Close the zlib stream.
    inflateEnd(&strm);
}

HMR_GZ_HANDLER *hmr_gz_open_read(const char *filepath)
{
    //Read the GZIP file.
    HMR_GZ_HANDLER *gz_handler = new HMR_GZ_HANDLER();
#ifdef _MSC_VER
    FILE *gz_file = NULL;
    fopen_s(&gz_file, filepath, "rb");
#else
    FILE *gz_file = fopen(filepath, "rb");
#endif
    if(!gz_file)
    {
        time_error(1, "Failed to open GZIP file %s", filepath);
    }
    gz_handler->gz_file = gz_file;
    //Allocate the processing queue, 3 for triple buffer.
    hmr_bin_queue_create(&(gz_handler->queue), 3);
    //Prepare the buffer.
    hmr_bin_buf_create(&gz_handler->buffer);
    //Start the GZIP parsing thread.
    gz_handler->parse_thread = std::thread(hmr_gzip_parse, gz_file, gz_handler->queue);
    //Provide the GZIP handler.
    return gz_handler;
}

void hmr_gz_close_read(HMR_GZ_HANDLER* gz_handler)
{
    //Check whether the queue is marked as finished.
    if (!gz_handler->queue->finish)
    {
        hmr_bin_queue_finish(gz_handler->queue);
    }
    //Wait for parse thread to complete.
    gz_handler->parse_thread.join();
    //Free the queue and buffer.
    hmr_bin_buf_free(gz_handler->buffer);
    hmr_bin_queue_free(gz_handler->queue);
    //Close the file.
    fclose(gz_handler->gz_file);
}
