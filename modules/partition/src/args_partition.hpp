#ifndef ARGS_PARTITION_H
#define ARGS_PARTITION_H

#include <vector>

typedef struct HMR_ARGS
{
    const char* nodes = NULL;
    const char* edge = NULL;
    const char* allele_table = NULL;
    const char* output = NULL;
    int groups = -1, allele_groups = -1, threads = 1;
} HMR_ARGS;

#endif // ARGS_PARTITION_H
