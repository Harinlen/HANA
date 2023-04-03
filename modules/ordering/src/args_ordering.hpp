#ifndef ARGS_PARTITION_H
#define ARGS_PARTITION_H

#include <cstdlib>
#include <cstdint>

typedef struct HMR_ARGS
{
    const char* nodes = NULL;
    const char* edge = NULL;
    const char* group = NULL;
    const char* output = NULL;
    double mutapb = 0.2;
    int ngen = 5000, npop = 100;
    uint64_t seed = 42;
} HMR_ARGS;

#endif // ARGS_PARTITION_H
