#ifndef ARGS_PARTITION_H
#define ARGS_PARTITION_H

#include <cstdlib>
#include <vector>

typedef struct HMR_ARGS
{
    const char* nodes = NULL;
    const char* edge = NULL;
    const char* allele_table = NULL;
    const char* output = NULL;
    int groups = -1, threads = 1, min_re = 25, max_link_density = 2;
} HMR_ARGS;

#endif // ARGS_PARTITION_H
