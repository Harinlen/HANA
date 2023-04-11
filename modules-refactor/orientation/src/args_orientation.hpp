#ifndef ARGS_ORIENTATION_H
#define ARGS_ORIENTATION_H

#include <cstdlib>
#include <cstdint>
#include <vector>

typedef struct HMR_ARGS
{
    const char* nodes = NULL;
    const char* reads = NULL;
    std::vector<char*> seq;
    int read_buffer_size = 512;
} HMR_ARGS;

#endif // ARGS_ORIENTATION_H
