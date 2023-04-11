#ifndef ARGS_DRAFT_H
#define ARGS_DRAFT_H

#include <cstdlib>
#include <vector>

typedef struct HMR_ARGS
{
    const char* nodes = NULL;
    const char* reads = NULL;
    const char* allele_table = NULL;
    const char* output = NULL;
    double max_density = 2.0;
    int min_links = 3, min_re = 10, read_buffer_size = 512;
} HMR_ARGS;

#endif // ARGS_DRAFT_H
