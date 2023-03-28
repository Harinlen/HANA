#ifndef ARGS_DUMP_H
#define ARGS_DUMP_H

#include <vector>

typedef struct HMR_ARGS
{
    char* bam_file = nullptr;
    char* output = nullptr;
    char* edge_file = nullptr;
    char* group_file = nullptr;
    int threads = 1;
} HMR_ARGS;

#endif // ARGS_DUMP_H
