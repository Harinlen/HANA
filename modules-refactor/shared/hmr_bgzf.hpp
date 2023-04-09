#ifndef HMR_BGZF_H
#define HMR_BGZF_H

#include <cstdio>
#include <thread>

typedef struct HMR_BIN_QUEUE HMR_BIN_QUEUE;
typedef struct HMR_BIN_DATA_BUF HMR_BIN_DATA_BUF;

/* BGFZ file handler */
typedef struct HMR_BGZF_HANDLER
{
    HMR_BIN_QUEUE* queue;
    HMR_BIN_DATA_BUF* buffer;
    FILE* bgzf_file;
    std::thread parse_thread;
} HMR_BGZF_HANDLER;

/* BGFZ file process functions */
HMR_BGZF_HANDLER* hmr_bgzf_open(const char* filepath, int threads = 1);
void hmr_bgzf_close(HMR_BGZF_HANDLER* bgzf_handler);

#endif // HMR_BGZF_H