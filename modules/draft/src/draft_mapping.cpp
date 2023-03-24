#include <cassert>

#include "hmr_contig_graph.hpp"
#include "hmr_global.hpp"

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
    HMR_UNUSED(length)
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

inline int32_t contig_id(MAPPING_DRAFT_USER* mapping_user, int32_t refId)
{
    return (refId > -1 && refId < mapping_user->contig_idx) ? mapping_user->contig_id_map[refId] : -1;
}

void mapping_draft_read_align(size_t id, const MAPPING_INFO& mapping_info, void* user)
{
    HMR_UNUSED(id)
    MAPPING_DRAFT_USER* mapping_user = reinterpret_cast<MAPPING_DRAFT_USER*>(user);
    //Check whether the mapping info is valid, then check whether the position is in range.
    int32_t ref_index = contig_id(mapping_user, mapping_info.refID),
        next_ref_index = contig_id(mapping_user, mapping_info.next_refID);
    // 3852 stands for:
    // - read unmapped (0x4)
    // - mate unmapped (0x8)
    // - not primary alignment (0x100)
    // - read fails platform/vendor quality checks (0x200)
    // - read is PCR or optical duplicate (0x400)
    // - supplementary alignment (0x800)
    if (ref_index == -1 || next_ref_index == -1 //Reference index invalid detection.
        || mapping_info.mapq == 255 //Map quality is invalid.
        || mapping_info.mapq < mapping_user->mapq //Check whether the mapping reaches the minimum quality
        || (mapping_info.flag & 3852) != 0 // Filtered source code.
        || ref_index == next_ref_index // We don't care about the pairs on the same contigs.
        || (!position_in_range(mapping_info.pos, mapping_user->contig_ranges[ref_index]))) // Or the position is not in the position.
    {
        return;
    }
    //Output the edge to the output buffer.
    if (mapping_user->output_size == mapping_user->output_offset)
    {
        //Flush the entire buffer to the file.
        fwrite(mapping_user->output_buffer, mapping_user->output_size, 1, mapping_user->reads_file);
        //Reset the offset back to beginning.
        mapping_user->output_offset = 0;
    }
    //Save the edges to buffer.
    HMR_MAPPING* output_read = reinterpret_cast<HMR_MAPPING*>(mapping_user->output_buffer + mapping_user->output_offset);
    (*output_read) = HMR_MAPPING{ ref_index, mapping_info.pos, next_ref_index, mapping_info.next_pos };
    mapping_user->output_offset += sizeof(HMR_MAPPING);
    //Construct the edge info.
    HMR_EDGE edge = hmr_graph_edge(ref_index, next_ref_index);
    //A new pair is found.
    auto& edge_map = mapping_user->edges;
    auto edge_iter = edge_map.find(edge.data);
    if (edge_iter == edge_map.end())
    {
        //Insert one count to the data.
        edge_map.insert(std::make_pair(edge.data, 1));
    }
    else
    {

        ++(edge_iter->second);
    }
}
