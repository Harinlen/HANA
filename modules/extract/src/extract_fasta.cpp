#include <cstring>

#include "hmr_global.hpp"

#include "extract_fasta.hpp"

static std::mutex enzyme_search_result_mutex;

void extract_enzyme_search_start(const ENZYME_VEC &enzymes, CANDIDATE_ENZYMES& search)
{
    //Assign the value to structure.
    search.reserve(enzymes.size());
    for(const std::string &enzyme: enzymes)
    {
        ENZYME_SEARCH_PARAM enzyme_search;
        enzyme_search.enzyme = enzyme.data();
        enzyme_search.enzyme_length = static_cast<int32_t>(enzyme.size());
        enzyme_search.offset = 0;
        search.push_back(enzyme_search);
    }
}

void extract_enzyme_search_end(CANDIDATE_ENZYMES &)
{
}

int32_t contig_draft_search(const char* seq, int32_t seq_size, ENZYME_SEARCH_PARAM& search)
{
    // Check range.
    if (search.offset + search.enzyme_length >= seq_size)
    {
        return -1;
    }
    // Perform searching.
    const char* result = strstr(seq + search.offset, search.enzyme);
    if (result == NULL)
    {
        return -1;
    }
    //Calculate the offset.
    return static_cast<int32_t>(result - seq);
}

int32_t contig_draft_search_next(const char* seq, int32_t seq_size, CANDIDATE_ENZYMES &search_param)
{
    //Get the next minimum value.
    int32_t next_valid_id = -1, next_valid_min = -1;
    for(size_t i=0; i<search_param.size(); ++i)
    {
        if(search_param[i].last_result != -1 &&
                next_valid_min < search_param[i].last_result)
        {
            next_valid_min = search_param[i].last_result;
            next_valid_id = static_cast<int32_t>(i);
        }
    }
    //If we still have a new id.
    if(next_valid_min != -1)
    {
        //Perform the search on the specific id.
        search_param[next_valid_id].offset = search_param[next_valid_id].last_result + search_param[next_valid_id].enzyme_length;
        search_param[next_valid_id].last_result = contig_draft_search(seq, seq_size, search_param[next_valid_id]);
    }
    return next_valid_min;
}

void contig_range_search(const ENZYME_RANGE_SEARCH& param)
{
    std::deque<ENZYME_RANGE> ranges;
    CANDIDATE_ENZYMES range_param = param.init_search_range;
    const char* seq = param.seq;
    const int32_t seq_size = static_cast<int32_t>(param.seq_size);
    const int32_t start_border = param.half_range, end_border = seq_size - param.half_range;
    //Do the init search.
    for (ENZYME_SEARCH_PARAM& search: range_param)
    {
        search.last_result = contig_draft_search(seq, seq_size, search);
    }
    //Fetch the enzyme position.
    int32_t enzyme_pos = contig_draft_search_next(seq, seq_size, range_param), counter = 0;
    while (enzyme_pos != -1)
    {
        //Increase the counter.
        ++counter;
        //Record the enzyme position.
        int32_t range_start = enzyme_pos, range_end = enzyme_pos;
        range_start = (range_start < start_border) ? 0 : range_start - start_border;
        range_end = (range_end > end_border) ? seq_size : range_end + start_border;
        //Check shall we merge the result to the last record.
        if (!ranges.empty() && range_start < ranges.back().end)
        {
            //Extend the last enzyme range.
            ranges.back().end = range_end;
        }
        else
        {
            //Record as a new range block.
            ranges.push_back(ENZYME_RANGE{ range_start, range_end });
        }
        //Perform the next search.
        enzyme_pos = contig_draft_search_next(seq, seq_size, range_param);
    }
    //Check range search param.
    if(!param.init_search_calc.empty())
    {
        //Perform weight counter calculator search.
        CANDIDATE_ENZYMES calc_param = param.init_search_calc;
        for (ENZYME_SEARCH_PARAM& search: calc_param)
        {
            search.last_result = contig_draft_search(seq, seq_size, search);
        }
        //Fetch the enzyme position.
        int32_t calc_enzyme_pos = contig_draft_search_next(seq, seq_size, calc_param);
        //Reset the counter.
        counter = 0;
        while (calc_enzyme_pos != -1)
        {
            //Increase the counter.
            ++counter;
            //Perform the next search.
            calc_enzyme_pos = contig_draft_search_next(seq, seq_size, calc_param);
        }
    }
    //Recover the sequence memory.
    free(param.seq);
    //Construct the enzyme range result.
    CONTIG_RANGE_RESULT contig_result;
    contig_result.contig_index = param.contig_index;
    contig_result.counter = counter;
    hDequeListToVector(ranges, contig_result.ranges);
    {
        enzyme_search_result_mutex.lock();
        param.results->push_back(contig_result);
        enzyme_search_result_mutex.unlock();
    }
}

void extract_fasta_search_proc(int32_t, char* name, size_t name_size, char* seq, size_t seq_size, void* user)
{
    EXTRACT_FASTA_USER *node_user = reinterpret_cast<EXTRACT_FASTA_USER*>(user);
    //Create the node information.
    int32_t contig_id = static_cast<int32_t>(node_user->nodes.size());
    node_user->nodes.push_back(HMR_NODE{ static_cast<int32_t>(seq_size), -1 });
    node_user->node_names.push_back(HMR_NODE_NAME{ static_cast<int32_t>(name_size), name });
    //Start to search the enzyme.
    node_user->pool.push_task(ENZYME_RANGE_SEARCH
        {
            node_user->init_search_range,
            node_user->init_search_calc,
            seq, 
            seq_size, 
            node_user->half_range, 
            contig_id, 
            &(node_user->results)
        });
}
