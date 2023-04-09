#ifndef HMR_GZ_H
#define HMR_GZ_H

#include <cstdio>
#include <thread>

typedef struct HMR_BIN_QUEUE HMR_BIN_QUEUE;
typedef struct HMR_BIN_DATA_BUF HMR_BIN_DATA_BUF;

typedef struct HMR_GZ_HANDLER
{
    HMR_BIN_QUEUE *queue;
    HMR_BIN_DATA_BUF *buffer;
    FILE* gz_file;
    std::thread parse_thread;
} HMR_GZ_HANDLER;

HMR_GZ_HANDLER *hmr_gz_open_read(const char *filepath);
void hmr_gz_close_read(HMR_GZ_HANDLER* gz_handler);

#endif // HMR_GZ_H
