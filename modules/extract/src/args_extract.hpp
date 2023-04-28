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
    int mapq = 40, threads = 1, range = 500, fasta_pool = 32, mapping_pool = 512, pairs_read_len = 150;
    bool skip_flag = false, skip_range = false;
} HMR_ARGS;

#endif // ARGS_EXTRACT_H
