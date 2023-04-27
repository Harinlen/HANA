#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cassert>

#include "hmr_args.hpp"
#include "hmr_contig_graph.hpp"
#include "hmr_ui.hpp"
#include "hmr_path.hpp"

#include "ordering_loader.hpp"
#include "ordering_ea.hpp"

#include "args_ordering.hpp"

extern HMR_ARGS opts;

int main(int argc, char* argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.nodes) { help_exit(-1, "Missing HMR graph contig information file path."); }
    if (!path_can_read(opts.nodes)) { time_error(-1, "Cannot read HMR graph contig file %s", opts.nodes); }
    if (!opts.edge) { help_exit(-1, "Missing HMR graph edge file path."); }
    if (!path_can_read(opts.edge)) { time_error(-1, "Cannot read HMR graph edge file %s", opts.edge); }
    if (!opts.group) { help_exit(-1, "Missing HMR contig group file path."); }
    if (!path_can_read(opts.group)) { time_error(-1, "Cannot read HMR contig group file %s", opts.group); }
    if (!opts.output) { help_exit(-1, "Missing HMR sequence output file path."); }
    time_print("Execution configuration:");
    time_print("\tMaximum idle generations: %d", opts.ngen);
    time_print("\tMaximum possible generations: %d", opts.max_gen);
    time_print("\tMutation probability: %lf", opts.mutapb);
    time_print("\tNum of candidate sequences: %d", opts.npop);
    if (opts.seed == 0)
    {
        //Use randome device to generate a randome seed.
        opts.seed = std::random_device()();
    }
    time_print("\tRandom seed: %lu", opts.seed);
    time_print("\tThreads: %d", opts.threads);
    time_print("\tEdge buffer: %dK", opts.read_buffer_size);
    opts.read_buffer_size <<= 10;
    //Read the group file for the index.
    time_print("Loading group contig index from %s", opts.group);
    HMR_CONTIG_ID_VEC contig_group;
    hmr_graph_load_contig_ids(opts.group, contig_group);
    std::sort(contig_group.begin(), contig_group.end());
    time_print("%zu contig indices loaded.", contig_group.size());
    //Load the contig information.
    time_print("Loading contig information from %s", opts.nodes);
    HMR_NODES contigs;
    hmr_graph_load_contigs(opts.nodes, contigs);
    time_print("%zu contig(s) loaded.", contigs.size());
    //Initialize the ordering info.
    time_print("Loading edge information from %s", opts.edge);
    ORDERING_INFO info;
    info.contig_size = static_cast<int32_t>(contig_group.size());
    info.edges.reserve(info.contig_size);
    ORDERING_EDGE_LOADER loader{ contig_group, info.edges };
    hmr_graph_load_edges(opts.edge, opts.read_buffer_size, ordering_edge_map_size_proc, ordering_edge_map_data_proc, &loader);
    time_print("Group edges are loaded.");
    //Shuffle the initial order for good luck.
    time_print("Generating initial contig orders...");
    std::mt19937_64 rng(opts.seed);
    std::shuffle(contig_group.begin(), contig_group.end(), rng);
    info.init_genome = static_cast<ORDERING_TIG*>(malloc(sizeof(ORDERING_TIG) * info.contig_size));
    if (!info.init_genome)
    {
        time_error(-1, "Failed to allocate memory for initial order sequence.");
    }
    assert(info.init_genome);
    //Loop for 2 phase, the first phase is generating a ordered sequence, the second phase is to test whether it could be better.
    for (int32_t phase = 0; phase < 2; ++phase)
    {
        time_print("Starting evaluation algorithm phase %d...", phase + 1);
        ordering_ea_init(contig_group, contigs, info);
        contig_group = ordering_ea_optimize(phase + 1, opts.npop, opts.ngen, opts.max_gen, opts.mutapb, info, rng, opts.threads);
    }
    free(info.init_genome);
    //Dump the data to output file.
    time_print("Writing ordered contig indices to %s", opts.output);
    hmr_graph_save_contig_ids(opts.output, contig_group);
    time_print("Ordering complete.");
    return 0;
}