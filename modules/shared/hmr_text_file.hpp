#ifndef HMR_TEXT_FILE_H
#define HMR_TEXT_FILE_H

#include "hmr_char.hpp"

#ifdef _MSC_VER
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

typedef struct TEXT_LINE_BUF
{
    char* buf;
    size_t buf_size, reserved, offset;
} TEXT_LINE_BUF;

/* Text line parsing function type. */
typedef ssize_t(*TEXT_GETLINE)(char** line, size_t* line_size, TEXT_LINE_BUF *buf, void* file_handle);

typedef struct TEXT_LINE_HANDLE
{
    TEXT_GETLINE parser;
    void* file_handle;
    TEXT_LINE_BUF buf;
} TEXT_LINE_HANDLE;

std::string text_open_read(const char *filepath, void** handle);
bool text_open_read_line(const char* filepath, TEXT_LINE_HANDLE *handle);
void text_close_read_line(TEXT_LINE_HANDLE* handle);

bool text_open_write(const char* filepath, FILE** handle);

#endif // HMR_TEXT_FILE_H
