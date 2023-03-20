#include <cassert>

#include "hmr_contig_graph.hpp"

#include "draft_mapping.hpp"

void mapping_draft_n_contig(uint32_t n_ref, void* user)
{
    MAPPING_DRAFT_USER* mapping_user = reinterpret_cast<MAPPING_DRAFT_USER*>(user);
    //Allocate the contig id maps.
    mapping_user->contig_id_map = new int32_t[n_ref];
    mapping_user->contig_idx = 0;
    //Prepare the output buffer. (<< 20) = 1 MiB
    mapping_user->output_size = sizeof(HMR_MAPPING) << 20;
    mapping_user->output_buffer = static_cast<char*>(malloc(mapping_user->output_size));
    mapping_user->output_offset = 0;
    //Check memory allocation.
    assert(mapping_user->output_buffer);
}

void mapping_draft_contig(uint32_t name_length, char* name, uint32_t length, void* user)
{
    MAPPING_DRAFT_USER* mapping_user = reinterpret_cast<MAPPING_DRAFT_USER*>(user);
    //Build the contig mapping index.
    std::string contig_name(name, name_length);
    auto name_finder = mapping_user->contig_ids.find(contig_name);
    int32_t contig_id = (name_finder == mapping_user->contig_ids.end()) ? -1 : name_finder->second;
    //Check the contig id.
    if (contig_id != -1)
    {
        //If the the contig id is in invalid set, then set the id to be -1.
        if (mapping_user->invalid_ids.find(contig_id) != mapping_user->invalid_ids.end())
        {
            contig_id = -1;
        }
    }
    //Record the contig id at the map.
    mapping_user->contig_id_map[mapping_user->contig_idx] = contig_id;
    ++mapping_user->contig_idx;
}

bool position_in_range(int32_t pos, const ENZYME_RANGES& ranges)
{
    for (size_t i = 0; i < ranges.length; ++i)
    {
        const ENZYME_RANGE& r = ranges.ranges[i];
        //Check whether the position is in the range.
        if (r.start <= pos && pos <= r.end)
        {
            return true;
        }
    }
    //No position matched.
    return false;
}

bool find_reads_pair(READ_POS_MAP& pos_map, const READ_POS& read_pos, const READ_POS& paired_read_pos)
{
    //Find the paired pos.
    auto finder = pos_map.find(paired_read_pos.data);
    if (finder == pos_map.end())
    {
        return false;
    }
    //We found the paired position, check whether the read pos is in the set.
    auto read_finder = finder->second.find(read_pos.data);
    if (read_finder == finder->second.end())
    {
        return false;
    }
    finder->second.erase(read_pos.data);
    return true;
}

inline int32_t contig_id(MAPPING_DRAFT_USER* mapping_user, int32_t refId)
{
    return (refId > -1 && refId < mapping_user->contig_idx) ? mapping_user->contig_id_map[refId] : -1;
}

void mapping_draft_read_align(size_t id, const MAPPING_INFO& mapping_info, void* user)
{
    MAPPING_DRAFT_USER* mapping_user = reinterpret_cast<MAPPING_DRAFT_USER*>(user);
    //Check whether the mapping info is valid, then check whether the position is in range.
    int32_t ref_index = contig_id(mapping_user, mapping_info.refID),
        next_ref_index = contig_id(mapping_user, mapping_info.next_refID);
    if (ref_index == -1 || next_ref_index == -1 //Reference index invalid detection.
        || mapping_info.mapq == 255 //Map quality is invalid
        //Position in range.
        || !position_in_range(mapping_info.pos, mapping_user->contig_ranges[ref_index]))
    {
        return;
    }
    //Construct the edge info.
    READ_POS read_pos{ ref_index, mapping_info.pos }, paired_read_pos{ next_ref_index, mapping_info.next_pos };
    HMR_EDGE edge = hmr_graph_edge(read_pos.read.id, paired_read_pos.read.id);
    //Create the edge information, no matter how.
    RAW_EDGE_MAP& edge_map = mapping_user->edges;
    auto edge_iter = edge_map.find(edge.data);
    if (edge_map.find(edge.data) == edge_map.end())
    {
        //Insert a new edge record.
        edge_map.insert(std::make_pair(edge.data, MAPPING_COUNT {1, 0}));
        //Update the edge iterator.
        edge_iter = edge_map.find(edge.data);
    }
    else
    {
        //Increase the edge counter.
        ++(edge_iter->second.pairs);
    }
    //Check whether the mapping reaches the minimum quality,
    if (mapping_info.mapq < mapping_user->mapq ||
        //the index are coming from the same contig.
        ref_index == next_ref_index)
    {
        return;
    }
    //The current mapping information matches the request, check whether the paired read is in the record.
    if (find_reads_pair(mapping_user->records, read_pos, paired_read_pos))
    {
        //Paired information are found.
        ++(edge_iter->second.qualified_pairs);
        //Save the edge to the buffer.
        if (mapping_user->output_offset == mapping_user->output_size)
        {
            fwrite(mapping_user->output_buffer, mapping_user->output_size, 1, mapping_user->reads_file);
            mapping_user->output_offset = 0;
        }
        //Construct and write the mapping info to the reads file.
        HMR_MAPPING* mapping = reinterpret_cast<HMR_MAPPING*>(mapping_user->output_buffer + mapping_user->output_offset);
        *mapping = HMR_MAPPING{ read_pos.read.id, mapping_info.pos, paired_read_pos.read.id, mapping_info.next_pos };
        mapping_user->output_offset += sizeof(HMR_MAPPING);
    }
    else
    {
        READ_POS_MAP& pos_map = mapping_user->records;
        //Insert the mapping pair info to pos map.
        auto finder = pos_map.find(read_pos.data);
        if (finder == pos_map.end())
        {
            //Create a new record.
            READ_POS_SET pos_set;
            pos_set.insert(paired_read_pos.data);
            pos_map.insert(std::make_pair(read_pos.data, pos_set));
        }
        else
        {
            finder->second.insert(paired_read_pos.data);
        }
    }
}
