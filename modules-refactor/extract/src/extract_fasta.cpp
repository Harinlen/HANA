#include "hmr_global.hpp"

#include "extract_fasta.hpp"

static std::mutex enzyme_search_result_mutex;

void extract_enzyme_search_start(const char* enzyme, const int32_t enzyme_length, ENZYME_SEARCH_PARAM& search)
{
    //Assign the value to structure.
    search.enzyme = enzyme;
    search.enzyme_length = enzyme_length;
    search.offset = 0;
}

void extract_enzyme_search_end(ENZYME_SEARCH_PARAM & search)
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

void contig_range_search(const ENZYME_RANGE_SEARCH& param)
{
    std::list<ENZYME_RANGE> ranges;
    ENZYME_SEARCH_PARAM search_param = param.init_search_param;
    const char* seq = param.seq;
    const int32_t seq_size = static_cast<int32_t>(param.seq_size);
    const int32_t start_border = param.half_range, end_border = seq_size - param.half_range;
    //Search starts here.
    int32_t enzyme_pos = contig_draft_search(seq, seq_size, search_param), counter = 0;
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
        //Update the offset.
        search_param.offset = enzyme_pos + search_param.enzyme_length;
        //Perform the next search.
        enzyme_pos = contig_draft_search(seq, seq_size, search_param);
    }
    //Construct the enzyme range result.
    CONTIG_RANGE_RESULT contig_result;
    contig_result.contig_index = param.contig_index;
    contig_result.counter = counter;
    hMoveListToVector(ranges, contig_result.ranges);
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
    node_user->nodes.push_back(HMR_CONTIG{ static_cast<int32_t>(seq_size), -1, static_cast<int32_t>(name_size), name });
    //Start to search the enzyme.
    node_user->pool.push_task(ENZYME_RANGE_SEARCH
        {
            node_user->init_search_param, 
            seq, 
            seq_size, 
            node_user->half_range, 
            contig_id, 
            &(node_user->results)
        });
}
