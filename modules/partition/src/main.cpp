#include <algorithm>
#include <cstdlib>
#include <cstdio>

#include "hmr_args.hpp"
#include "hmr_ui.hpp"
#include "hmr_path.hpp"
#include "hmr_contig_graph.hpp"
#include "partition_type.hpp"
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
    //Building clusters for the contigs.
    CLUSTER_INFO cluster_user;
    partition_create_init_clusters(contigs, cluster_user);
    {
        //Read the edge information.
        time_print("Loading edge information from %s", opts.edge);
        HMR_EDGE_MAP graph_edges;
        hmr_graph_restore_edge_map(opts.edge, &graph_edges);
        time_print("Contig edges loaded.");
        //Generate link densities based on graph edge information.
        partition_init_link_densities(contigs, graph_edges, cluster_user, opts.min_re, opts.max_link_density, opts.output);
    }
    //Start clustering.
    time_print("Start clustering...");
    partition_cluster(cluster_user, opts.groups);
    time_print("%zu group(s) finally generated.", cluster_user.cluster_size);
    //Dumping the result to the output file.
    char group_file_path[4096];
    for(size_t i=0; i<cluster_user.cluster_size; ++i)
    {
        CONTIG_ID_VECTOR &cluster = *cluster_user.clusters[i];
        //Construct the file path.
        sprintf(group_file_path, "%s_%zug%zu.hmr_group", opts.output, cluster_user.cluster_size, i+1);
        time_print("Saving cluster %zu (%zu contigs) to %s...", i+1, cluster.size(), group_file_path);
        //Sort the partition.
        std::sort(cluster.begin(), cluster.end());
        hmr_graph_save_partition(group_file_path, cluster);
    }
    time_print("Partition complete.");
    return 0;
}
