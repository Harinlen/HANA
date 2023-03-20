#include <cstring>

#include "draft_fasta_index.hpp"

void contig_draft_search_start(const char* enzyme, int32_t enzyme_length, ENZYME_SEARCH& search)
{
    //Assign the value to structure.
    search.enzyme = enzyme;
    search.enzyme_length = enzyme_length;
    search.offset = 0;
}

void contig_draft_search_end(ENZYME_SEARCH& search)
{
    search.enzyme_length = 0;
}

int32_t contig_draft_search(const char* seq, size_t seq_size, ENZYME_SEARCH* search)
{
    // Perform searching.
    const char *result = strstr(seq + search->offset, search->enzyme);
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
    //Get the search parameter.
    ENZYME_SEARCH* search = param.search;
    char* seq = param.seq;
    int32_t seq_size = param.seq_size, offset = 0;
    int32_t enzyme_pos = contig_draft_search(seq, seq_size, search);
    const int32_t half_range = param.range, end_range = param.seq_size - half_range;
    size_t counter = 0;
    while (enzyme_pos != -1)
    {
        //Increase the counter.
        ++counter;
        //Record the enzyme position.
        int32_t range_start = offset + enzyme_pos, range_end = range_start;
        //Calculate the range end.
        range_start = (range_start < half_range) ? 0 : range_start - half_range;
        range_end = (range_end > end_range) ? param.seq_size : (range_end + half_range);
        //Check shall we merged to last ranges.
        if (!ranges.empty() && range_start <= ranges.back().end)
        {
            //Update the back result.
            ranges.back().end = range_end;
        }
        else
        {
            //Append the new range.
            ranges.push_back(ENZYME_RANGE{ range_start, range_end });
        }
        //Update the offset.
        offset += enzyme_pos + search->enzyme_length;
        //Search the next position.
        enzyme_pos = contig_draft_search(seq + offset, seq_size - offset, search);
    }
    //Convert the enzyme range to array.
    ENZYME_RANGES& chain_ranges = param.chain_node->data;
    chain_ranges.counter = counter;
    chain_ranges.length = ranges.size();
    chain_ranges.ranges = static_cast<ENZYME_RANGE*>(malloc(sizeof(ENZYME_RANGE) * ranges.size()));
    if (!chain_ranges.ranges)
    {
        time_error(-1, "No enough memory for chain range allocation.");
    }
    int32_t range_index = 0;
    for (auto i = ranges.begin(); i != ranges.end(); ++i)
    {
        chain_ranges.ranges[range_index] = *i;
        ++range_index;
    }
    //Free the sequence.
    free(param.seq);
}

void contig_draft_build(int32_t index, char* seq_name, size_t seq_name_size, char* seq, size_t seq_size, void* user)
{
    DRAFT_NODES_USER* node_user = reinterpret_cast<DRAFT_NODES_USER*>(user);
    //Append the sequence information.
    node_user->nodes->push_back(HMR_CONTIG{ static_cast<int32_t>(seq_size), static_cast<int32_t>(seq_name_size), seq_name });
    //Create the search chain.
    ENZYME_RANGE_CHAIN* chain_node = new ENZYME_RANGE_CHAIN();
    chain_node->next = NULL;
    //Append the chain node to the chain.
    if (node_user->chain_head == NULL)
    {
        //Then set the head to be the node.
        node_user->chain_tail = node_user->chain_head = chain_node;
    }
    else
    {
        //Append the chain node to the chail tail.
        node_user->chain_tail->next = chain_node;
        node_user->chain_tail = chain_node;
    }
    //Push the search request into search pool.
    node_user->pool->push_task(ENZYME_RANGE_SEARCH{ node_user->search, chain_node, seq, static_cast<int32_t>(seq_size), node_user->range });
}
