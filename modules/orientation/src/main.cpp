#include <cassert>
#include <cstdlib>
#include <cstdio>

#include "hmr_args.hpp"
#include "hmr_contig_graph.hpp"
#include "hmr_path.hpp"
#include "hmr_ui.hpp"

#include "args_orientation.hpp"
#include "orientation.hpp"

extern HMR_ARGS opts;

int main(int argc, char *argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.nodes) { help_exit(-1, "Missing HMR graph contig information file path."); }
    if (!path_can_read(opts.nodes)) { time_error(-1, "Cannot read HMR graph contig file %s", opts.nodes); }
    if (!opts.reads) { help_exit(-1, "Missing HMR paired-reads file path."); }
    if (!path_can_read(opts.reads)) { time_error(-1, "Cannot read HMR paired-reads file %s", opts.reads); }
    if (opts.seq.empty()) { help_exit(-1, "Missing HMR ordered sequence file path."); }
    for (const char* seq_path : opts.seq)
    {
        if (!path_can_read(seq_path)) { time_error(-1, "Cannot read HMR ordered sequence path file %s", seq_path); }
    }
    //if (!opts.output) { help_exit(-1, "Missing HMR chromosome sequence output file path."); }
    time_print("Execution configuration:");
    time_print("\tPaired-reads buffer: %dK", opts.read_buffer_size);
    opts.read_buffer_size <<= 10;
    time_print("\tOptimized sequences: %zu", opts.seq.size());
    //Load the sequence file.
    ORIENTATION_INFO info;
    {
        time_print("Loading contig information from %s", opts.nodes);
        HMR_NODES nodes;
        hmr_graph_load_contigs(opts.nodes, nodes);
        time_print("%zu contig(s) loaded.", nodes.size());
        //Load the sequence file.
        time_print("Loading sequence indices...");
        orientation_init(opts.seq, nodes, info);
        time_print("%zu sequence(s) loaded.", opts.seq.size());
    }
    //Loading the reads and parse the sequence, it automatically find the best orientation.
    time_print("Loading reads pair information from %s", opts.reads);
    hmr_graph_load_reads(opts.reads, opts.read_buffer_size, orientation_calc_gradient, &info);
    time_print("Reads information loaded, direction gradient calculated.");
    //Based on the gradient, extract the direction.
    time_print("Extracting direction results...");
    std::vector<CHROMOSOME_CONTIGS> chromosomes;
    chromosomes.reserve(info.sequences.size());
    //Extract the chromosome sequence from the result.
    for (const auto& chromosome_sequence : info.sequences)
    {
        chromosomes.push_back(orientation_extract(chromosome_sequence));
    }
    time_print("%zu sequence of orientation generated.", chromosomes.size());
    //Dump the data to the output file.
    for (size_t i = 0; i < opts.seq.size(); ++i)
    {
        std::string chromo_path = hmr_graph_path_chromo_name(opts.seq[i]);
        const auto& chromo_data = chromosomes[i];
        time_print("Writing %zu records to %s", chromo_data.size(), chromo_path.c_str());
        hmr_graph_save_chromosome(chromo_path.c_str(), chromo_data);
    }
    time_print("Orientation complete.");
    return 0;
}
