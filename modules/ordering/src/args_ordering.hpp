#ifndef ARGS_PARTITION_H
#define ARGS_PARTITION_H

#include <cstdlib>

typedef struct HMR_ARGS
{
    const char* nodes = NULL;
    const char* edge = NULL;
    const char* group = NULL;
    const char* output = NULL;
} HMR_ARGS;

#endif // ARGS_PARTITION_H
