#include <algorithm>
#include <cstdlib>
#include <cstdio>

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
    time_print("Generating initial genomes...");
    order_info.contig_size = static_cast<int32_t>(order_info.contig_group.size());
    static const int my_order[] = {65, 61, 22, 74, 27, 75, 47, 19, 64, 63, 59, 13, 9, 67, 28, 23, 20, 58, 57, 72, 42, 60, 16, 50, 56, 71, 40, 48, 51, 1, 62, 33, 14, 77, 6, 4, 78, 2, 30, 53, 5, 69, 24, 41, 52, 34, 8, 45, 29, 25, 36, 38, 43, 3, 7, 37, 39, 66, 17, 46, 70, 0, 21, 49, 26, 18, 35, 54, 15, 76, 11, 44, 31, 68, 10, 73, 32, 12, 55};
    int my_idx[79];
    for(auto i=0; i<79; ++i)
    {
        my_idx[i] = order_info.contig_group[my_order[i]];
    }
    order_info.init_genome = ordering_init(order_info, true, my_idx);
    //Start to reduce the gradient of the groups.
    for(int32_t phase=1; phase<3; ++phase)
    {
        //Optimize the current phase.
        ordering_optimize_phase(phase, opts.npop, opts.ngen, opts.mutapb, order_info);
        break;
    }
    return 0;
}
