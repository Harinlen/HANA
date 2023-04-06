#include <cassert>
#include <cstdlib>
#include <cstdio>

#include "hmr_args.hpp"
#include "hmr_contig_graph.hpp"
#include "hmr_path.hpp"
#include "hmr_ui.hpp"

#include "args_orientation.hpp"

extern HMR_ARGS opts;

constexpr uint8_t DIRECTION_POSITIVE = 0;
constexpr uint8_t DIRECTION_NEGATIVE = 1;

typedef struct ORIENTATION_INFO
{
    CONTIG_ID_VECTOR contig_id_seq;
    std::unordered_map<int32_t, int32_t> contig_id_map;
    double *contig_size;
    double *cost_matrix[2];
    double *cost_buffer;
    int32_t contig_id;
} ORIENTATION_INFO;

void orientation_contig_load(int32_t length, int32_t, int32_t, char* , void* user)
{
    ORIENTATION_INFO *info = reinterpret_cast<ORIENTATION_INFO *>(user);
    //Check whether the contig id is in the sequence.
    auto iter = info->contig_id_map.find(info->contig_id);
    if(iter != info->contig_id_map.end())
    {
        //Find out the position in the sequence.
        info->contig_size[iter->second] = static_cast<double>(length);
    }
    //Increase the id for the next contig.
    ++info->contig_id;
}

void orientation_reads_pair_load(const HMR_MAPPING *mapping, void *user)
{
    ORIENTATION_INFO *info = reinterpret_cast<ORIENTATION_INFO *>(user);
    //Find whether all the contig ids are in the sequence.
    auto iter_a = info->contig_id_map.find(mapping->refID);
    if(iter_a == info->contig_id_map.end())
    {
        return;
    }
    auto iter_b = info->contig_id_map.find(mapping->next_refID);
    if(iter_b == info->contig_id_map.end())
    {
        return;
    }
    //Check the position a and position b.
    int32_t pos_a = iter_a->second, pos_b = iter_b->second;
    if(pos_a > pos_b)
    {
        info->cost_matrix[DIRECTION_POSITIVE][pos_a] += info->contig_size[pos_a] - mapping->pos;
        info->cost_matrix[DIRECTION_NEGATIVE][pos_a] += mapping->pos;
        info->cost_matrix[DIRECTION_POSITIVE][pos_b] += mapping->next_pos;
        info->cost_matrix[DIRECTION_NEGATIVE][pos_b] += info->contig_size[pos_b] - mapping->next_pos;
    }
    else
    {
        info->cost_matrix[DIRECTION_POSITIVE][pos_a] += mapping->pos;
        info->cost_matrix[DIRECTION_NEGATIVE][pos_a] += info->contig_size[pos_a] - mapping->pos;
        info->cost_matrix[DIRECTION_POSITIVE][pos_b] += info->contig_size[pos_b] - mapping->next_pos;
        info->cost_matrix[DIRECTION_NEGATIVE][pos_b] += mapping->next_pos;
    }
}

int main(int argc, char *argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.nodes) { help_exit(-1, "Missing HMR graph contig information file path."); }
    if (!path_can_read(opts.nodes)) { time_error(-1, "Cannot read HMR graph contig file %s", opts.nodes); }
    if (!opts.reads) { help_exit(-1, "Missing HMR paired-reads file path."); }
    if (!path_can_read(opts.reads)) { time_error(-1, "Cannot read HMR paired-reads file %s", opts.reads); }
    if (!opts.seq) { help_exit(-1, "Missing HMR paired-reads file path."); }
    if (!path_can_read(opts.reads)) { time_error(-1, "Cannot read HMR paired-reads file %s", opts.reads); }
    if (!opts.output) { help_exit(-1, "Missing HMR chromosome sequence output file path."); }
    if (opts.buffer_size < 1) { help_exit(-1, "Buffer size should be greater than 1."); }
    time_print("Execution configuration:");
    time_print("\tReads-pair buffer units: %d", opts.buffer_size);
    //Load the sequence file.
    ORIENTATION_INFO info;
    time_print("Loading group contig index from %s", opts.seq);
    hmr_graph_load_partition(opts.seq, info.contig_id_seq);
    time_print("%zu contig indices loaded.", info.contig_id_seq.size());
    //Prepare the direction cost matrix.
    info.cost_buffer = static_cast<double *>(malloc(sizeof(double) * (info.contig_id_seq.size() << 1)));
    if(!info.cost_buffer)
    {
        time_error(-1, "Failed to allocate memory for the cost matrix.");
    }
    assert(info.cost_buffer);
    for(size_t i=0, i_max = (info.contig_id_seq.size() << 1); i<i_max; ++i)
    {
        info.cost_buffer[i] = 0.0;
    }
    info.cost_matrix[DIRECTION_POSITIVE] = info.cost_buffer;
    info.cost_matrix[DIRECTION_NEGATIVE] = info.cost_buffer + info.contig_id_seq.size();
    //Prepare the memory for the contig length.
    info.contig_size = static_cast<double *>(malloc(sizeof(double) * info.contig_id_seq.size()));
    if(!info.contig_size)
    {
        time_error(-1, "Failed to allocate memory for contig sizes.");
    }
    //Build the contig id-position map.
    for(size_t i=0; i<info.contig_id_seq.size(); ++i)
    {
        info.contig_id_map.insert(std::make_pair(info.contig_id_seq[i], static_cast<int32_t>(i)));
    }
    //Load the contig information.
    time_print("Loading contig information from %s", opts.nodes);
    info.contig_id = 0;
    hmr_graph_load_contigs(opts.nodes, NULL, orientation_contig_load, &info);
    time_print("%zu contig(s) loaded.", info.contig_id);
    //Load the reads information.
    time_print("Loading reads pair information from %s", opts.reads);
    hmr_graph_load_reads(opts.reads, static_cast<size_t>(opts.buffer_size), orientation_reads_pair_load, &info);
    time_print("Reads information loaded, direction calculated.");
    //Find out the best direction of the contig expected.
    CHROMOSOME_CONTIGS chromosome_contigs;
    chromosome_contigs.reserve(info.contig_id_seq.size());
    for(size_t i=0; i<info.contig_id_seq.size(); ++i)
    {
        chromosome_contigs.push_back(HMR_DIRECTED_CONTIG
                                     {
                                         info.contig_id_seq[i],
                                         info.cost_matrix[DIRECTION_NEGATIVE][i] < info.cost_matrix[DIRECTION_POSITIVE][i]
                                     });
    }
    //Dump the data to the output file.
    time_print("Writing direction to %s", opts.output);
    hmr_graph_save_chromosome(opts.output, chromosome_contigs);
    time_print("Orientation complete.");
    free(info.contig_size);
    free(info.cost_buffer);
    return 0;
}
