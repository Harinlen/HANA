#ifndef EXTRACT_FASTA_H
#define EXTRACT_FASTA_H

#include <deque>

#include "hmr_contig_graph_type.hpp"
#include "hmr_thread_pool.hpp"

#include "extract_fasta_type.hpp"

/* enzyme searching */
typedef struct ENZYME_SEARCH_PARAM
{
    const char* enzyme;
    int32_t enzyme_length;
    int32_t offset;
    int32_t last_result;
} ENZYME_SEARCH_PARAM;

typedef std::vector<ENZYME_SEARCH_PARAM> CANDIDATE_ENZYMES;

void extract_enzyme_search_start(const std::vector<char*> &enzymes, CANDIDATE_ENZYMES& search);
void extract_enzyme_search_end(CANDIDATE_ENZYMES& );

typedef struct CONTIG_RANGE_RESULT
{
    CONTIG_ENZYME_RANGE ranges;
    int32_t counter;
    int32_t contig_index;
} CONTIG_RANGE_RESULT;

typedef std::deque<CONTIG_RANGE_RESULT> CONTIG_RANGE_RESULTS;

typedef struct ENZYME_RANGE_SEARCH
{
    CANDIDATE_ENZYMES init_search_param;
    char* seq;
    size_t seq_size;
    int32_t half_range;
    int32_t contig_index;
    CONTIG_RANGE_RESULTS* results;
} ENZYME_RANGE_SEARCH;

void contig_range_search(const ENZYME_RANGE_SEARCH& param);
typedef hmr::thread_pool<ENZYME_RANGE_SEARCH> RANGE_SEARCH_POOL;

typedef std::deque<HMR_NODE> CONTIG_CHAIN;
typedef std::deque<HMR_NODE_NAME> CONTIG_NAME_CHAIN;

typedef struct EXTRACT_FASTA_USER
{
    const CANDIDATE_ENZYMES& init_search_param;
    CONTIG_CHAIN& nodes;
    CONTIG_NAME_CHAIN& node_names;
    RANGE_SEARCH_POOL& pool;
    const int32_t half_range;
    CONTIG_RANGE_RESULTS& results;
} EXTRACT_FASTA_USER;

void extract_fasta_search_proc(int32_t, char* name, size_t name_size, char* seq, size_t seq_size, void* user);

#endif // EXTRACT_FASTA_H
