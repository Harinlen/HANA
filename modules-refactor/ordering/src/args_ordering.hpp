#ifndef ARGS_ORDERING_H
#define ARGS_ORDERING_H

#include <cstdlib>
#include <cstdint>

typedef struct HMR_ARGS
{
    const char* nodes = NULL;
    const char* edge = NULL;
    const char* group = NULL;
    const char* output = NULL;
    double mutapb = 0.2;
    int ngen = 5000, npop = 100, threads = 1, read_buffer_size = 512, max_gen = 1000000;
    uint64_t seed = 0;
} HMR_ARGS;

#endif // ARGS_ORDERING_H
