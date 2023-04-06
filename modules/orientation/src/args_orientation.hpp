#ifndef ARGS_ORIENTATION_H
#define ARGS_ORIENTATION_H

#include <cstdlib>
#include <cstdint>

typedef struct HMR_ARGS
{
    const char* nodes = NULL;
    const char* reads = NULL;
    const char* seq = NULL;
    const char* output = NULL;
    int buffer_size = 65535;
} HMR_ARGS;

#endif // ARGS_ORIENTATION_H
