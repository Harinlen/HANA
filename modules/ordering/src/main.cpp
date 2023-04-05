#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <random>

#include "hmr_args.hpp"
#include "hmr_ui.hpp"
#include "hmr_path.hpp"
#include "hmr_contig_graph.hpp"
#include "ordering_type.hpp"

#include "args_ordering.hpp"
#include "ordering_loader.hpp"
#include "ordering_descent.hpp"

extern HMR_ARGS opts;

int main(int argc, char *argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.nodes) { help_exit(-1, "Missing HMR graph contig information file path."); }
    if (!path_can_read(opts.nodes)) { time_error(-1, "Cannot read HMR graph contig file %s", opts.nodes); }
    if (!opts.edge) { help_exit(-1, "Missing HMR graph edge weight file path."); }
    if (!path_can_read(opts.edge)) { time_error(-1, "Cannot read HMR graph edge weight file %s", opts.edge); }
    if (!opts.group) { help_exit(-1, "Missing HMR contig group file path."); }
    if (!path_can_read(opts.group)) { time_error(-1, "Cannot read HMR contig group file %s", opts.group); }
    if (!opts.output) { help_exit(-1, "Missing HMR sequence output file path."); }
    time_print("Execution configuration:");
    time_print("\tMaximum idle generations: %d", opts.ngen);
    time_print("\tMutation probability: %lf", opts.mutapb);
    time_print("\tNum of candidate sequences: %d", opts.npop);
    if (opts.seed == 0)
    {
        //Use randome device to generate a randome seed.
        opts.seed = std::random_device()();
    }
    time_print("\tRandom seed: %lu", opts.seed);
    time_print("\tThreads: %d", opts.threads);
    //Read the group file for the index.
    ORDERING_INFO order_info;
    time_print("Loading group contig index from %s", opts.group);
    hmr_graph_load_partition(opts.group, order_info.contig_group);
    std::sort(order_info.contig_group.begin(), order_info.contig_group.end());
    time_print("%zu contig indices loaded.", order_info.contig_group.size());
    //Load the contig information.
    time_print("Loading contig information from %s", opts.nodes);
    hmr_graph_load_contigs(opts.nodes, ordering_contig_size_proc, ordering_contig_node_proc, &order_info);
    time_print("%zu contig(s) loaded.", order_info.contigs.size());
    time_print("Loading edge information from %s", opts.edge);
    hmr_graph_load_edge(opts.edge, ordering_edge_map_size_proc, ordering_edge_map_data_proc, &order_info);
    time_print("Group edges are loaded.");
    //Create the initial ordering.
    time_print("Generating initial contig orders...");
    order_info.contig_size = static_cast<int32_t>(order_info.contig_group.size());
    //Shuffle the initial order for good luck.
    std::mt19937_64 rng(opts.seed);
    std::shuffle(order_info.contig_group.begin(), order_info.contig_group.end(), rng);
    //Start to reduce the gradient of the groups.
    time_print("Running evolutionary algorithm...");
    order_info.init_genome = ordering_init_alloc(order_info.contig_size);
    for(int32_t phase=1; phase<3; ++phase)
    {
        //Optimize the current phase.
        time_print("Starting phase %d...", phase);
        ordering_init(order_info, order_info.contig_group);
        ordering_optimize_phase(phase, opts.npop, opts.ngen, opts.mutapb, order_info, rng, opts.threads);
    }
    free(order_info.init_genome);
    //Dump the data to output file.
    time_print("Writing ordered contig indices to %s", opts.output);
    hmr_graph_save_partition(opts.output, order_info.contig_group);
    return 0;
}
