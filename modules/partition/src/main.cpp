#include <algorithm>

#include "hmr_args.hpp"
#include "hmr_path.hpp"
#include "hmr_ui.hpp"

#include "args_partition.hpp"
#include "hmr_contig_graph.hpp"

#include "partition.hpp"

extern HMR_ARGS opts;

int main(int argc, char* argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.nodes) { help_exit(-1, "Missing HMR graph contig information file path."); }
    if (!path_can_read(opts.nodes)) { time_error(-1, "Cannot read HMR graph contig file %s", opts.nodes); }
    if (!opts.edges) { help_exit(-1, "Missing HMR graph edge weight file path."); }
    if (!path_can_read(opts.edges)) { time_error(-1, "Cannot read HMR graph edge weight file %s", opts.edges); }
    if (opts.groups < 1) { time_error(-1, "Please specify the group to be separated."); }
    if (!opts.output) { help_exit(-1, "Missing output HMR partition file path."); }
    //Print the execution configuration.
    time_print("Execution configuration:");
    time_print("\tNumber of Partitions: %d", opts.groups);
    time_print("\tAllele mode: %s", opts.allele ? "Yes" : "No");
    time_print("\tEdge buffer: %dK", opts.read_buffer_size);
    opts.read_buffer_size <<= 10;
    time_print("\tNon informative ratio: %d", opts.non_informative_ratio);
    //Load the contig node information.
    HMR_NODES nodes;
    time_print("Loading contig information from %s", opts.nodes);
    hmr_graph_load_contigs(opts.nodes, nodes);
    time_print("%zu contig(s) information loaded.", nodes.size());
    HMR_CONTIG_ID_VEC invalid_nodes;
    {
        std::string invalid_node_path = hmr_graph_path_nodes_invalid(opts.nodes);
        if (path_can_read(invalid_node_path.c_str()))
        {
            //Load the invalid nodes.
            time_print("Loading invalid contig id from %s", invalid_node_path.c_str());
            hmr_graph_load_contig_ids(invalid_node_path.c_str(), invalid_nodes);
            time_print("%zu invalid id(s) loaded.", invalid_nodes.size());
        }
    }
    //Load the allele table when needed.
    HMR_ALLELE_MAP allele_map;
    if (opts.allele)
    {
        time_print("Loading allele table from %s", opts.allele);
        HMR_CONTIG_ID_TABLE allele_table;
        hmr_graph_load_contig_table(opts.allele, allele_table);
        //Convert the allele table into allele map.
        hmr_graph_allele_map_init(allele_map, allele_table);
        time_print("Allele table loaded.");
    }
    //Construct the partition clusters.
    CLUSTER_INFO partition_info;
    time_print("Initialize the partition information...");
    partition_init_clusters(nodes, invalid_nodes, partition_info);
    if (opts.allele)
    {
        partition_info.allele_map = &allele_map;
    }
    hmr_graph_load_edges(opts.edges, opts.read_buffer_size, partition_edge_size_proc, partition_edge_proc, &partition_info);
    time_print("%zu merge operations built.", partition_info.merge_size);
    //Start clustering.
    time_print("Clustering %zu informative contigs with target of %d groups...", partition_info.cluster_size, opts.groups);
    partition_cluster(partition_info, opts.groups);
    time_print("Merge stage complete, %zu cluster(s) left.", partition_info.cluster_size);
    //Ignore the clusters only have 1 contig.
    time_print("Filtering individual node clusters...");
    std::vector<HMR_CONTIG_ID_VEC*> clusters;
    clusters.reserve(partition_info.cluster_size);
    for (size_t i = 0; i < partition_info.cluster_size; ++i)
    {
        //Ignore the individual node cluster.
        if (partition_info.clusters[i]->size() == 1)
        {
            int32_t contig_id = (*partition_info.clusters[i])[0];
            //Add contig id to invalid contig.
            invalid_nodes.emplace_back(contig_id);
            //Reset the belongs flag.
            partition_info.belongs[contig_id] = NULL;
            free(partition_info.clusters[i]);
            continue;
        }
        //Save the clusters.
        clusters.emplace_back(partition_info.clusters[i]);
    }
    partition_info.cluster_size = 0;
    time_print("%zu clusters remain.", clusters.size());
    //Try to recover previously skipped contigs.
    time_print("Recovering skipped contigs...");
    //Find out the best matched cluster, but do not add them in.
    std::sort(invalid_nodes.begin(), invalid_nodes.end());
    partition_recover(clusters, invalid_nodes, opts.non_informative_ratio, partition_info);
    time_print("Gathering contig clusters...");
    for (int32_t contig_id : invalid_nodes)
    {
        //Add the contig id to the cluster groups.
        if (partition_info.belongs[contig_id])
        {
            partition_info.belongs[contig_id]->emplace_back(contig_id);
        }
    }
    //Dumping the result to the output file.
    int32_t cluster_total = static_cast<int32_t>(clusters.size());
    for (int32_t i = 0; i < cluster_total; ++i)
    {
        std::string cluster_path = hmr_graph_path_cluster_name(opts.output, i+1, cluster_total);
        HMR_CONTIG_ID_VEC* group_cluster = clusters[i];
        int32_t total_re = 0;
        int64_t total_bp = 0;
        for(const int32_t contig_id: *group_cluster)
        {
            total_re += nodes[contig_id].enzyme_count;
            total_bp += nodes[contig_id].length;
        }
        time_print("Saving cluster %zu (%zu contigs, %d RE total, avg 1 per %d bp) to %s", i + 1, group_cluster->size(), total_re, total_bp / total_re, cluster_path.c_str());
        std::sort(clusters[i]->begin(), clusters[i]->end());

        hmr_graph_save_contig_ids(cluster_path.c_str(), *clusters[i]);
    }
    partition_free_clusters(partition_info);
    time_print("Partition complete.");
    return 0;
}
