#ifndef ARGS_BUILD_H
#define ARGS_BUILD_H

#include <cstdlib>
#include <vector>

typedef struct HMR_ARGS
{
    const char* fasta = NULL;
    std::vector<char*> chromosomes;
    const char* output = NULL;
} HMR_ARGS;

#endif // ARGS_BUILD_H
