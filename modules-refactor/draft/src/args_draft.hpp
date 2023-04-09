#ifndef ARGS_DRAFT_H
#define ARGS_DRAFT_H

#include <cstdlib>
#include <vector>

typedef struct HMR_ARGS
{
    const char* contigs = NULL;
    const char* reads = NULL;
    const char* allele_table = NULL;
    const char* output = nullptr;
    int enzyme_nuc_length = 0, mapq = 40, threads = 1, range = 500, min_enzymes = 0, fasta_pool = 32, mapping_pool = 1;
} HMR_ARGS;

#endif // ARGS_DRAFT_H
