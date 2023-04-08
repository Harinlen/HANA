#include <algorithm>
#include <cstdlib>
#include <cstdio>

#include "hmr_args.hpp"
#include "hmr_ui.hpp"
#include "hmr_path.hpp"
#include "hmr_contig_graph.hpp"
#include "hmr_allele.hpp"
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
    bool allele_mode = opts.allele_table;
    time_print("Execution configuration:");
    time_print("\tNumber of Partitions: %d", opts.groups);
    time_print("\tAllele mode: %s", allele_mode ? "Yes" : "No");
    //Load the contig node information.
    HMR_CONTIGS contigs;
    time_print("Loading contig information from %s", opts.nodes);
    hmr_graph_restore_contig_data(opts.nodes, &contigs);
    time_print("%zu contig(s) information loaded.", contigs.size());
    //Building clusters for the contigs.
    CLUSTER_INFO cluster_user;
    //Read the allele table when provided.
    cluster_user.allele_mode = allele_mode;
    if(cluster_user.allele_mode)
    {
        time_print("Loading allele table from %s", opts.allele_table);
        hmr_allele_table_conflict_set(opts.allele_table, contigs, cluster_user.allele_table);
        time_print("Allele table loaded.");
    }
    partition_create_init_clusters(contigs, cluster_user);
    {
        //Prepare the link density user.
        LINK_FACTORS contig_factors{ NULL, 0.0, CONTIG_ID_SET(), allele_mode, cluster_user.allele_table };
        partition_init_factors(static_cast<int32_t>(contigs.size()), contig_factors);
        //Read the edge information.
        time_print("Calculating contig factors...");
        hmr_graph_load_edge(opts.edge, partition_factors_size_proc, partition_factors_edge_proc, &contig_factors);
        contig_factors.avergage_links = contig_factors.avergage_links / static_cast<double>(contigs.size()) * 2.0;
        time_print("Contig factors ready.");
        time_print("Skipping contigs likely from repetitive regions with multiplicity %d...", opts.max_link_density);
        int32_t num_of_contigs = static_cast<int32_t>(contigs.size());
        for (int32_t i = 0; i < num_of_contigs; ++i)
        {
            contig_factors.factors[i] /= contig_factors.avergage_links;
            if (contig_factors.factors[i] >= static_cast<double>(opts.max_link_density))
            {
                contig_factors.skipped.insert(i);
            }
        }
        //Generate link densities based on graph edge information.
        time_print("Skipping contigs whose enzyme count is less than %d...", opts.min_re);
        {
            CONTIG_ID_SET skipped = partition_skip_few_res(contigs, opts.min_re);
            contig_factors.skipped.insert(skipped.begin(), skipped.end());
        }
        time_print("Initiating merge requests and link densities...");
        DENSITY_INFO density_info{ contig_factors, cluster_user };
        cluster_user.link_densities.resize(num_of_contigs);
        hmr_graph_load_edge(opts.edge, partition_link_densities_size_proc, partition_link_densities_edge_proc, &density_info);
        time_print("%zu merge requests created.", cluster_user.merge_size);
        time_print("Skipping empty link contigs...");
        partition_skip_empty_links(cluster_user.link_densities, contig_factors.skipped);
        time_print("%zu contig(s) are skipped total.", contig_factors.skipped.size());
        //Delete all the skipped sets.
        if (!contig_factors.skipped.empty())
        {
            //Remove skipped nodes.
            partition_remove_skipped_contigs(contig_factors, cluster_user, opts.output);
        }
        partition_free_factors(contig_factors);
    }
    //Start clustering.
    time_print("Start clustering %zu groups...", cluster_user.cluster_size);
    partition_cluster(cluster_user, opts.groups);
    time_print("%zu group(s) finally generated.", cluster_user.cluster_size);
    //Dumping the result to the output file.
    char group_file_path[4097];
    for(size_t i=0; i<cluster_user.cluster_size; ++i)
    {
        CONTIG_ID_VECTOR &cluster = *cluster_user.clusters[i];
        //Construct the file path.
#ifdef _MSC_VER
        sprintf_s(group_file_path, 4096, "%s_%zug%zu.hmr_group", opts.output, cluster_user.cluster_size, i + 1);
#else
        sprintf(group_file_path, "%s_%zug%zu.hmr_group", opts.output, cluster_user.cluster_size, i+1);
#endif
        time_print("Saving cluster %zu (%zu contigs) to %s...", i+1, cluster.size(), group_file_path);
        //Sort the partition.
        std::sort(cluster.begin(), cluster.end());
        hmr_graph_save_partition(group_file_path, cluster);
    }
    time_print("Partition complete.");
    return 0;
}
