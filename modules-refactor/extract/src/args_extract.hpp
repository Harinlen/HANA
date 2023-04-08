#ifndef ARGS_EXTRACT_H
#define ARGS_EXTRACT_H

#include <cstdlib>
#include <vector>

typedef struct HMR_ARGS
{
    const char* fasta = NULL;
    const char* output = NULL;
    std::vector<char*> mappings;
    char* enzyme = nullptr;
    int mapq = 40, threads = 1;
} HMR_ARGS;

#endif // ARGS_EXTRACT_H
