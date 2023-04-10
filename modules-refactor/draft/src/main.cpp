#include <algorithm>

#include "hmr_algorithm.hpp"
#include "hmr_args.hpp"
#include "hmr_contig_graph.hpp"
#include "hmr_global.hpp"
#include "hmr_path.hpp"
#include "hmr_ui.hpp"

#include "args_draft.hpp"
#include "draft_mappings.hpp"

extern HMR_ARGS opts;

int main(int argc, char* argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.nodes) { help_exit(-1, "Missing HMR graph contig information file path."); }
    if (!path_can_read(opts.nodes)) { time_error(-1, "Cannot read HMR graph contig file %s", opts.nodes); }
    if (!opts.reads) { help_exit(-1, "Missing HMR paired-reads file path."); }
    if (!path_can_read(opts.reads)) { time_error(-1, "Cannot read HMR paired-reads file %s", opts.reads); }
    //Print the execution configuration.
    bool allele_mode = opts.allele_table;
    time_print("Execution configuration:");
    time_print("\tAllele mode: %s", allele_mode ? "Yes" : "No");
    time_print("\tMinimum edge links: %d", opts.min_links);
    time_print("\tMinimum RE sites: %d", opts.min_re);
    time_print("\tMaximum link density: %.2lf", opts.max_density);
    time_print("\tPaired-reads buffer: %dK", opts.read_buffer_size);
    opts.read_buffer_size <<= 10;
    //Load the contig node information.
    HMR_NODES nodes;
    HMR_NODE_NAMES node_names;
    time_print("Loading contig information from %s", opts.nodes);
    hmr_graph_load_contigs(opts.nodes, nodes, &node_names);
    time_print("%zu contig(s) information loaded.", nodes.size());
    int32_t num_of_contigs = static_cast<int32_t>(nodes.size());
    //Load the invalid node information when necessary.
    std::vector<int32_t> contig_invalid;
    contig_invalid.resize(num_of_contigs);
    //Filter out the minimum REs.
    time_print("Skip contig(s) with few RE sites and finding the maximum RE site count...");
    int32_t max_enzyme_count = 0;
    for (int32_t i = 0; i < num_of_contigs; ++i)
    {
        //Check the enzyme count reach the limits.
        contig_invalid[i] = (nodes[i].enzyme_count < opts.min_re);
        //Find out the maximum enzyme count.
        max_enzyme_count = hMax(max_enzyme_count, nodes[i].enzyme_count);
    }
    time_print("Maximum enzyme counts in contigs: %d", max_enzyme_count);
    //Load the allele table when necessary.
    HMR_ALLELE_TABLE allele_table;
    if (allele_mode)
    {
        time_print("Loading allele table from %s", opts.allele_table);
        hmr_graph_load_allele_table(opts.allele_table, allele_table);
        time_print("%zu allele record(s) loaded.", allele_table.size());
    }
    //Read through the reads file, calculate the pairs.
    time_print("Reading paired-reads from %s", opts.reads);
    HMR_EDGE_COUNTERS edges;
    std::vector<double> node_factors;
    {
        EDGE_BUILDER edge_builder{ allele_table, EDGE_COUNTER() };
        hmr_graph_load_reads(opts.reads, opts.read_buffer_size, draft_mappings_build_edges, &edge_builder);
        time_print("Paired-reads have been counted, %zu edge(s) generated.", edge_builder.counter.size());
        //Filter out the valid links.
        time_print("Removing edges failed to reach the minimum links...");
        time_print("Calculating the node repetitive factors...");
        edges.reserve(edge_builder.counter.size());
        node_factors.resize(num_of_contigs);
        for (int32_t i = 0; i < num_of_contigs; ++i)
        {
            node_factors[i] = 0.0;
        }
        double links_average = 0.0;
        int64_t max_re_square = hSquare(static_cast<int64_t>(max_enzyme_count));
        for (const auto& iter : edge_builder.counter)
        {
            //Skip the edge.
            int32_t edge_num_of_links = iter.second;
            if (edge_num_of_links < opts.min_links)
            {
                continue;
            }
            //Create the edge, count to the factors.
            HMR_EDGE edge {};
            edge.data = iter.first;
            // Calculat the edge weights..
            double edge_weights = static_cast<double>(max_re_square * edge_num_of_links / nodes[edge.pos.start].enzyme_count / nodes[edge.pos.end].enzyme_count);
            edges.emplace_back(HMR_EDGE_INFO{ edge.pos.start, edge.pos.end, static_cast<uint64_t>(edge_num_of_links), edge_weights });
            edges.emplace_back(HMR_EDGE_INFO{ edge.pos.end, edge.pos.start, static_cast<uint64_t>(edge_num_of_links), edge_weights });
            node_factors[edge.pos.start] += edge_weights;
            node_factors[edge.pos.end] += edge_weights;
            links_average += edge_weights;
        }
        //Calculate the links average.
        links_average = 2.0 * links_average / static_cast<double>(num_of_contigs);
        time_print("%zu edge(s) generated, average links = %.2lf", edges.size(), links_average);
        time_print("Skip contig(s) with maximum multiplicity...");
        //Calculate the node factors.
        int32_t invalid_counter = 0;
        FILE* factor_skip;
        fopen_s(&factor_skip, "E:\\Downloads\\sampleDataForALLHiC\\Sspon-hind3\\hana_contig_factor_skip.txt", "w");
        for (int32_t i = 0; i < num_of_contigs; ++i)
        {
            if (contig_invalid[i])
            {
                ++invalid_counter;
                continue;
            }
            //Calculate the node factors.
            node_factors[i] /= links_average;
            //Check whether the node factors is valid or not.
            contig_invalid[i] = (static_cast<int32_t>(node_factors[i]) >= opts.max_density);
            if (contig_invalid[i])
            {
                ++invalid_counter;
            }
        }
        fclose(factor_skip);
        //Extract invalid node vector.
        HMR_CONTIG_ID_VEC invalid_nodes;
        invalid_nodes.reserve(invalid_counter);
        for (int32_t i = 0; i < num_of_contigs; ++i)
        {
            if (contig_invalid[i])
            {
                invalid_nodes.emplace_back(i);
            }
        }
        //Loop for all the edges, update their weights.
        time_print("%zu node(s) are skipped.", invalid_nodes.size());
    }
    //Based on the node factors, change the graph to directed-graph.
    time_print("Adjust link densities by repetitive factors...");
    size_t edge_removed = 0;
    for (size_t i = 0; i < edges.size(); ++i)
    {
        auto& edge_info = edges[i];
        if (contig_invalid[edge_info.start] || contig_invalid[edge_info.end])
        {
            ++edge_removed;
            continue;
        }
        edge_info.weights /= node_factors[edge_info.start];
        if (edge_removed)
        {
            //Need to move the current info to n previous.
            edges[i - edge_removed] = edge_info;
        }
    }
    edges.resize(edges.size() - edge_removed);
    time_print("%zu edge(s) generated.", edges.size());
    //Write the edge information to the target file.
    std::string edge_path = hmr_graph_path_edge(opts.output);
    time_print("Saving edge information to %s", edge_path.c_str());
    hmr_graph_save_edges(edge_path.c_str(), edges);
    time_print("Draft complete.");

    FILE* edge_dump;
    fopen_s(&edge_dump, "E:\\Downloads\\sampleDataForALLHiC\\Sspon-hind3\\hana_edge_info.txt", "w");
    for (const auto& edge : edges)
    {
        fprintf(edge_dump, "%s\t%s\t%lld\n", node_names[edge.start].name, node_names[edge.end].name, static_cast<int64_t>(edge.weights));
    }
    fclose(edge_dump);

    return 0;
}