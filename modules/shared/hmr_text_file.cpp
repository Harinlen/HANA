#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cstring>
#include <cassert>

#include "hmr_path.hpp"
#include "hmr_gz.hpp"
#include "hmr_bin_queue.hpp"
#include "hmr_ui.hpp"

#include "hmr_text_file.hpp"

#define TEXT_LINE_BUFFER_SIZE    (33554432) // 32MB

ssize_t text_getline_file(char** line, size_t* line_size, TEXT_LINE_BUF* line_buf, void* file_handle)
{
    //Check whether the buffer has and data rest.
    if (line_buf->buf_size > line_buf->offset)
    {
        //Try to find '\n' inside the buffer.
        size_t residual = line_buf->buf_size - line_buf->offset;
        char* buf = line_buf->buf + line_buf->offset,
            * pos = static_cast<char*>(memchr(buf, '\n', residual));
        if (pos != NULL)
        {
            //We found the '\n'
            *line = buf;
            *line_size = pos - buf + 1;
            //Update the buffer offset.
            line_buf->offset += *line_size;
            //Get the line size.
            return *line_size;
        }
        //Move the residual data to the front of the buffer.
        memcpy(line_buf->buf, buf, residual);
        //Reset the offset and buffer size.
        line_buf->buf_size = residual;
        line_buf->offset = 0;
    }
    //Need to fetch the data.
    FILE* fp = static_cast<FILE*>(file_handle);
    char* target = line_buf->buf + line_buf->buf_size;
    for(;;)
    {
        //Fill the buffer.
        size_t inc_size = fread(target, 1, line_buf->reserved - line_buf->buf_size, fp);
        if (inc_size == 0)
        {
            break;
        }
        //Add the readed bytes into buffer size.
        line_buf->buf_size += inc_size;
        //Search the '\n' char in the new buffer area.
        char *pos = static_cast<char*>(memchr(target, '\n', inc_size));
        if (pos != NULL)
        {
            //We found the line end!
            *line = line_buf->buf;
            *line_size = pos - line_buf->buf;
            line_buf->offset = *line_size;
            return *line_size;
        }
        //Double the buffer.
        size_t needed_buffer = line_buf->reserved << 1;
        char* needed = static_cast<char*>(realloc(line_buf->buf, needed_buffer));
        assert(needed);
        line_buf->buf = needed;
        line_buf->reserved = needed_buffer;
        //Update the target.
        target = line_buf->buf + line_buf->buf_size;
    }
    //Check whether we have any data in the buffer.
    if (line_buf->buf_size > line_buf->offset)
    {
        //Pass the rest of the data.
        *line = target;
        *line_size = line_buf->buf_size;
        line_buf->offset = *line_size;
        return *line_size;
    }
    return -1;
}

ssize_t text_getline_gz(char** line, size_t* line_size, TEXT_LINE_BUF* line_buf, void* file_handle)
{
    //Check whether the buffer has and data rest.
    if (line_buf->buf_size > line_buf->offset)
    {
        //Try to find '\n' inside the buffer.
        size_t residual = line_buf->buf_size - line_buf->offset;
        char* buf = line_buf->buf + line_buf->offset,
            * pos = static_cast<char*>(memchr(buf, '\n', residual));
        if (pos != NULL)
        {
            //We found the '\n'
            *line = buf;
            *line_size = pos - buf + 1;
            //Update the buffer offset.
            line_buf->offset += *line_size;
            //Get the line size.
            return *line_size;
        }
        //Move the residual data to the front of the buffer.
        memcpy(line_buf->buf, buf, residual);
        //Reset the offset and buffer size.
        line_buf->buf_size = residual;
        line_buf->offset = 0;
    }
    else
    {
        //Free the line buffer.
        free(line_buf->buf);
        line_buf->buf = NULL;
        line_buf->buf_size = 0;
        line_buf->offset = 0;
    }
    //Get the GZ handler.
    HMR_GZ_HANDLER* gz_handle = reinterpret_cast<HMR_GZ_HANDLER*>(file_handle);
    auto queue = gz_handle->queue;
    //Fetch the data from the queue.
    char* target;
    for (;;)
    {
        HMR_BIN_SLICE bin_slice = hmr_bin_queue_pop(queue);
        if (bin_slice.data == NULL)
        {
            //Reach the end of the data, still cannot fulfill, failed to fetch.
            break;
        }
        //Append the data to the buffer.
        if (line_buf->buf == NULL)
        {
            //Directly assign the slice data as buffer.
            line_buf->buf = bin_slice.data;
            line_buf->buf_size = bin_slice.data_size;
            target = line_buf->buf;
        }
        else
        {
            //Extend the size of the current buf.
            size_t expected_size = line_buf->buf_size + bin_slice.data_size;
            if (line_buf->reserved < expected_size)
            {
                char* needed = static_cast<char*>(realloc(line_buf->buf, expected_size));
                assert(needed);
                line_buf->buf = needed;
                line_buf->reserved = expected_size;
            }
            //Update the target.
            target = line_buf->buf + line_buf->buf_size;
            //Copy the data.
            memcpy(target, bin_slice.data, bin_slice.data_size);
            //Update the buffer size.
            line_buf->buf_size = expected_size;
            //Free the data slice.
            free(bin_slice.data);
        }
        //Search for '\n'.
        char* pos = static_cast<char*>(memchr(target, '\n', bin_slice.data_size));
        if (pos != NULL)
        {
            //We found the line end!
            *line = line_buf->buf;
            *line_size = pos - line_buf->buf;
            line_buf->offset = *line_size;
            return *line_size;
        }
    }
    //Check whether we have any data in the buffer.
    if (line_buf->buf && line_buf->buf_size > line_buf->offset)
    {
        //Pass the rest of the data.
        *line = line_buf->buf;
        *line_size = line_buf->buf_size;
        line_buf->offset = *line_size;
        return *line_size;
    }
    return -1;
}

std::string text_open_read(const char *filepath, void **handle)
{
    //Check the arguments.
    std::string suffix = path_suffix(filepath);
    if (suffix == ".gz")
    {
        //Use gzip module to open the file.
        *handle = hmr_gz_open_read(filepath);
        return "gz";
    }
    //Open as a normal text file.
#ifdef _MSC_VER
    FILE* text_file = NULL;
    if (fopen_s(&text_file, filepath, "r") != 0)
    {
        return "";
    }
#else
    FILE* text_file = fopen(filepath, "r");
    if (text_file == NULL)
    {
        return "";
    }
#endif
    *handle = text_file;
    return "txt";
}

bool text_open_read_line(const char* filepath, TEXT_LINE_HANDLE* handle)
{
    //Initial the buffer status.
    TEXT_LINE_BUF& buf = handle->buf;
    buf.buf = NULL;
    buf.offset = 0;
    buf.buf_size = 0;
    //Check the arguments.
    std::string mode = text_open_read(filepath, &(handle->file_handle));
    if (mode == "gz")
    {
        //Use gzip readline to read the file.
        handle->parser = text_getline_gz;
        buf.buf = NULL;
        buf.buf_size = 0;
        buf.reserved = 0;
        return true;
    }
    if (mode == "txt")
    {
        //Use text mode to read the file.
        handle->parser = text_getline_file;
        //Allocate the line reading buffer.
        buf.buf = static_cast<char*>(malloc(TEXT_LINE_BUFFER_SIZE));
        buf.buf_size = 0;
        buf.reserved = TEXT_LINE_BUFFER_SIZE;
        return true;
    }
    return false;
}

void text_close_read_line(TEXT_LINE_HANDLE *handle)
{
    //Recover the buffer.
    TEXT_LINE_BUF& buf = handle->buf;
    if (buf.buf)
    {
        free(buf.buf);
    }
    //Close the file handle.
    if(text_getline_file == handle->parser)
    {
        fclose(static_cast<FILE *>(handle->file_handle));
        return;
    }
    if (text_getline_gz == handle->parser)
    {
        hmr_gz_close_read(static_cast<HMR_GZ_HANDLER*>(handle->file_handle));
        return;
    }
    time_error(-1, "Unknown text getline handler to close.");
}

bool text_open_write(const char* filepath, FILE** handle)
{
    //Open as a normal text file.
#ifdef _MSC_VER
    FILE* text_file = NULL;
    if (fopen_s(&text_file, filepath, "w") != 0)
    {
        return false;
    }
#else
    FILE* text_file = fopen(filepath, "w");
    if (text_file == NULL)
    {
        return false;
    }
#endif
    * handle = text_file;
    return true;
}
