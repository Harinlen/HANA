#ifndef ARGS_DUMP_H
#define ARGS_DUMP_H

#include <cstdlib>
#include <vector>

typedef struct HMR_ARGS
{
    const char* fasta = NULL;
    const char* output = NULL;
    const char* nodes = NULL;
    const char* allele = NULL;
} HMR_ARGS;

#endif // ARGS_DUMP_H
