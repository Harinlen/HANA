#ifndef EXTRACT_FASTA_H
#define EXTRACT_FASTA_H

#include <list>

#include "hmr_contig_graph_type.hpp"
#include "hmr_thread_pool.hpp"

#include "extract_fasta_type.hpp"

/* enzyme searching */
typedef struct ENZYME_SEARCH_PARAM
{
    const char* enzyme;
    int32_t enzyme_length;
    int32_t offset;
} ENZYME_SEARCH_PARAM;

void extract_enzyme_search_start(const char* enzyme, const int32_t enzyme_length, ENZYME_SEARCH_PARAM& search);
void extract_enzyme_search_end(ENZYME_SEARCH_PARAM& search);

typedef struct CONTIG_RANGE_RESULT
{
    CONTIG_ENZYME_RANGE ranges;
    int32_t counter;
    int32_t contig_index;
} CONTIG_RANGE_RESULT;

typedef std::list<CONTIG_RANGE_RESULT> CONTIG_RANGE_RESULTS;

typedef struct ENZYME_RANGE_SEARCH
{
    ENZYME_SEARCH_PARAM init_search_param;
    char* seq;
    size_t seq_size;
    int32_t half_range;
    int32_t contig_index;
    CONTIG_RANGE_RESULTS* results;
} ENZYME_RANGE_SEARCH;

void contig_range_search(const ENZYME_RANGE_SEARCH& param);
typedef hmr::thread_pool<ENZYME_RANGE_SEARCH> RANGE_SEARCH_POOL;

typedef std::list<HMR_CONTIG> CONTIG_CHAIN;

typedef struct EXTRACT_FASTA_USER
{
    const ENZYME_SEARCH_PARAM& init_search_param;
    CONTIG_CHAIN& nodes;
    RANGE_SEARCH_POOL& pool;
    const int32_t half_range;
    CONTIG_RANGE_RESULTS& results;
} EXTRACT_FASTA_USER;

void extract_fasta_search_proc(int32_t, char* name, size_t name_size, char* seq, size_t seq_size, void* user);

#endif // EXTRACT_FASTA_H