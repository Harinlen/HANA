#ifndef ARGS_EXTRACT_H
#define ARGS_EXTRACT_H

#include <cstdlib>
#include <vector>

typedef struct HMR_ARGS
{
    const char* fasta = NULL;
    const char* output = NULL;
    const char* allele = NULL;
    std::vector<char*> mappings;
    char* enzyme = nullptr;
    int mapq = 40, threads = 1, range = 500, min_enzymes = 0, fasta_pool = 32, mapping_pool = 512;
} HMR_ARGS;

#endif // ARGS_EXTRACT_H
