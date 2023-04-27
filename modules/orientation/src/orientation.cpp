#include <cassert>

#include "hmr_contig_graph.hpp"
#include "hmr_ui.hpp"

#include "orientation.hpp"

void orientation_init(std::vector<char*> seq_paths, const HMR_NODES& nodes, ORIENTATION_INFO& info)
{
    //Expand the info.
    info.belongs.resize(nodes.size());
    info.sequences.resize(seq_paths.size());
    //Clear the belongs.
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        info.belongs[i] = -1;
    }
    //Load the sequences from path.
    for (int32_t i = 0, total_seqs = static_cast<int32_t>(seq_paths.size()); i < total_seqs; ++i)
    {
        time_print("Loading sequence from %s", seq_paths[i]);
        auto &seq = info.sequences[i];
        hmr_graph_load_contig_ids(seq_paths[i], seq.contig_id_seq);
        //Set belongs of the sequence.
        for (int32_t contig_id : seq.contig_id_seq)
        {
            info.belongs[contig_id] = i;
        }
        //Allocate memory for the sequences.
        const size_t buffer_size = seq.contig_id_seq.size() << 1;
        seq.cost_buffer = static_cast<double*>(malloc(sizeof(double) * buffer_size));
        if (!seq.cost_buffer)
        {
            time_error(-1, "Failed to allocate memory for the cost matrix.");
        }
        assert(seq.cost_buffer);
        for (size_t j = 0; j < buffer_size; ++j)
        {
            seq.cost_buffer[j] = 0.0;
        }
        seq.cost_matrix[DIRECTION_POSITIVE] = seq.cost_buffer;
        seq.cost_matrix[DIRECTION_NEGATIVE] = seq.cost_buffer + seq.contig_id_seq.size();
        //Prepare the memory for the contig length.
        seq.contig_size = static_cast<double*>(malloc(sizeof(double) * seq.contig_id_seq.size()));
        if (!seq.contig_size)
        {
            time_error(-1, "Failed to allocate memory for contig sizes.");
        }
        assert(seq.contig_size);
        //Build the contig->index hash map for storage.
        for (size_t j = 0; j < seq.contig_id_seq.size(); ++j)
        {
            seq.contig_id_map.insert(std::make_pair(seq.contig_id_seq[j], static_cast<int32_t>(j)));
            seq.contig_size[j] = static_cast<double>(nodes[seq.contig_id_seq[j]].length);
        }
    }
}

void orientation_calc_gradient(HMR_MAPPING* mapping, int32_t buf_size, void* user)
{
    ORIENTATION_INFO* info = static_cast<ORIENTATION_INFO*>(user);
    for (int32_t i = 0; i < buf_size; ++i)
    {
        const HMR_MAPPING& pair = mapping[i];
        //Find whether all the contig ids are in the sequence.
        if (info->belongs[pair.refID] != info->belongs[pair.next_refID] || info->belongs[pair.refID] == -1)
        {
            continue;
        }
        //Check the position of contig a and b.
        auto& seq = info->sequences[info->belongs[pair.refID]];
        int32_t pos_a_in_seq = seq.contig_id_map.find(pair.refID)->second,
            pos_b_in_seq = seq.contig_id_map.find(pair.next_refID)->second;
        if (pos_a_in_seq > pos_b_in_seq)
        {
            seq.cost_matrix[DIRECTION_POSITIVE][pos_a_in_seq] += static_cast<double>(pair.pos);
            seq.cost_matrix[DIRECTION_NEGATIVE][pos_a_in_seq] += seq.contig_size[pos_a_in_seq] - static_cast<double>(pair.pos);
            seq.cost_matrix[DIRECTION_POSITIVE][pos_b_in_seq] += seq.contig_size[pos_b_in_seq] - static_cast<double>(pair.next_pos);
            seq.cost_matrix[DIRECTION_NEGATIVE][pos_b_in_seq] += static_cast<double>(pair.next_pos);
        }
        else
        {
            seq.cost_matrix[DIRECTION_POSITIVE][pos_a_in_seq] += seq.contig_size[pos_a_in_seq] - static_cast<double>(pair.pos);
            seq.cost_matrix[DIRECTION_NEGATIVE][pos_a_in_seq] += static_cast<double>(pair.pos);
            seq.cost_matrix[DIRECTION_POSITIVE][pos_b_in_seq] += static_cast<double>(pair.next_pos);
            seq.cost_matrix[DIRECTION_NEGATIVE][pos_b_in_seq] += seq.contig_size[pos_b_in_seq] - static_cast<double>(pair.next_pos);
        }
    }
}

CHROMOSOME_CONTIGS orientation_extract(const ORIENTATION_SEQUENCE& sequence)
{
    const auto& ids = sequence.contig_id_seq;
    const auto& cost_matrix = sequence.cost_matrix;
    CHROMOSOME_CONTIGS result;
    result.reserve(ids.size());
    for (size_t i = 0; i < ids.size(); ++i)
    {
        result.push_back(HMR_DIRECTED_CONTIG{ ids[i], cost_matrix[DIRECTION_NEGATIVE][i] < cost_matrix[DIRECTION_POSITIVE][i] });
    }
    return result;
}
