#include <algorithm>
#include <cstdlib>
#include <cstdio>

#include "hmr_args.hpp"
#include "hmr_ui.hpp"
#include "hmr_path.hpp"
#include "hmr_contig_graph.hpp"
#include "partition.hpp"

#include "args_partition.hpp"

extern HMR_ARGS opts;

int main(int argc, char* argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.nodes) { help_exit(-1, "Missing HMR graph contig information file path."); }
    if (!path_can_read(opts.nodes)) { time_error(-1, "Cannot read HMR graph contig file %s", opts.nodes); }
    if (!opts.edge) { help_exit(-1, "Missing HMR graph edge weight file path."); }
    if (!path_can_read(opts.edge)) { time_error(-1, "Cannot read HMR graph edge weight file %s", opts.edge); }
    if (opts.groups < 1) { time_error(-1, "Please specify the group to be separated."); }
    if (!opts.output) { help_exit(-1, "Missing output HMR partition file path."); }
    //Print the execution configuration.
    time_print("Execution configuration:");
    time_print("\tNumber of Partitions: %d", opts.groups);
    time_print("\tThreads: %d", opts.threads);
    time_print("\tAllele mode: %s", opts.allele_table ? "Yes" : "No");
    //Load the contig node information.
    HMR_CONTIGS contigs;
    time_print("Loading contig information from %s", opts.nodes);
    hmr_graph_restore_contig_data(opts.nodes, &contigs);
    time_print("%zu contig(s) information loaded.", contigs.size());
    //Read the edge information.
    time_print("Loading edge information from %s", opts.edge);
    HMR_EDGE_COUNTERS edge_weights;
    hmr_graph_restore_edges(opts.edge, &edge_weights);
    time_print("Contig edges loaded.");
    //Sort the edge values.
    time_print("Sorting edges based on their weights...");
    std::sort(edge_weights.begin(), edge_weights.end(),
              [](const HMR_EDGE_INFO &l, const HMR_EDGE_INFO &r)
    {
        return l.weights > r.weights;
    });
    time_print("%zu edges sorted.", edge_weights.size());
    //Partition mission start.
    time_print("Dividing contigs into %d groups...", opts.groups);
    CONTIG_ID_VECTORS group_result = partition_contigs(contigs, edge_weights, opts.groups);
    time_print("Partition complete.");
    char output_path[4096];
    int32_t group_count = static_cast<int32_t>(group_result.size());
    for(int32_t i=0; i<group_count; ++i)
    {
        time_print("Writing group %zu result (%zu contigs)...", i+1, group_result[i].size());
        sprintf(output_path, "%s_%dg%d.hmr_group", opts.output, group_count, i+1);
        hmr_graph_save_partition(output_path, group_result[i]);
    }
    time_print("Partition complete.");
    return 0;
}
