#include <cstdlib>
#include <cstdio>

#include "hmr_args.hpp"
#include "hmr_ui.hpp"
#include "hmr_path.hpp"
#include "hmr_contig_graph.hpp"

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
    if (opts.allele_groups > -1 && opts.allele_groups < 2) { time_error(-1, "Allele groups must be greater than 2."); }
    if (opts.allele_groups > 0 && (!opts.allele_table)) { time_error(-1, "Allele table must be provided for allele group division."); }
    if (!opts.output) { help_exit(-1, "Missing output HMR partition file path."); }
    //Print the execution configuration.
    time_print("Execution configuration:");
    time_print("\tNumber of Partitions: %d", opts.groups);
    time_print("\tThreads: %d", opts.threads);
    time_print("\tAllele mode: %s", opts.allele_groups > 0 ? "Yes" : "No");
    //Load the contig node information.
    HMR_CONTIGS contigs;
    time_print("Loading contig information from %s", opts.nodes);
    hmr_graph_restore_contig_data(opts.nodes, &contigs);
    time_print("%zu contig(s) information loaded.", contigs.size());
    //Read the edge information.
    time_print("Loading edge information from %s", opts.edge);
    HMR_EDGE_COUNTERS edges;
    
    partition_load_edges(opts.edge, edges);
    time_print("Contig edges loaded.");
    //Partition mission start.
    time_print("Dividing contigs into %d groups...", opts.groups);
    auto partition_result = partition_run(contigs, edges, opts.groups, opts.allele_groups, allele_table, opts.threads);
    time_print("%zu group(s) of contigs generated.", partition_result.size());
    //Dump the result to output file.

    FILE* test = fopen("E:\\result.txt", "w");
    for (const auto &i : partition_result)
    {
        fprintf(test, "----------\n");
        for (const auto& j : i[0])
        {
            fprintf(test, "%s\n", contigs[j].name);
        }
    }
    fclose(test);
    return 0;
}