#ifndef ARGS_PARTITION_H
#define ARGS_PARTITION_H

#include <cstdlib>
#include <vector>

typedef struct HMR_ARGS
{
    const char* nodes = NULL;
    const char* edges = NULL;
    const char* allele = NULL;
    const char* output = NULL;
    int groups = -1, read_buffer_size = 512, non_informative_ratio = 3;
} HMR_ARGS;

#endif // ARGS_PARTITION_H
