#include <algorithm>
#include <cstring>
#include <cmath>

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
    time_print("Loading contig information from %s", opts.nodes);
    hmr_graph_load_contigs(opts.nodes, nodes);
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
    //Read through the reads file, calculate the pairs.
    time_print("Reading paired-reads from %s", opts.reads);
    HMR_EDGE_COUNTERS edges;
    std::vector<double> node_factors;
    {
        EDGE_COUNTER edge_pairs;
        //Load the allele table when necessary.
        if (allele_mode)
        {
            //First build the edge counter map.
            EDGE_COUNT_MAP edge_pair_map;
            hmr_graph_load_reads(opts.reads, opts.read_buffer_size, draft_mappings_build_edge_counter, &edge_pair_map);
            time_print("Paired-reads map has been built.");
            //Read the allele table.
            HMR_CONTIG_ID_TABLE allele_table;
            time_print("Loading allele table from %s", opts.allele_table);
            hmr_graph_load_contig_table(opts.allele_table, allele_table);
            time_print("%zu allele record(s) loaded.", allele_table.size());
            //Loop through the allele table, remove the edge pairs.
            time_print("Removing allele edges...");
            //Remove the inter-connected edges inside a row.
            for(auto &conflict_row: allele_table)
            {
                int32_t row_size = static_cast<int32_t>(conflict_row.size());
                if(row_size < 2)
                {
                    continue;
                }
                //Sort the row for better operation.
                //Remove the record inter-connected edges.
                for(int32_t i=0; i<row_size-1; ++i)
                {
                    for(int32_t j=i+1; j<row_size; ++j)
                    {
                        draft_mappings_remove_edge(edge_pair_map, conflict_row[i], conflict_row[j]);
                    }
                }
            }
            //For edges between rows, leave only the strongest connected edges.
            size_t num_of_allele_rows = allele_table.size();
            for(size_t i=0; i<num_of_allele_rows - 1; ++i)
            {
                for(size_t j=i+1; j<num_of_allele_rows; ++j)
                {
                    //Copy two allele row.
                    HMR_CONTIG_ID_VEC row_i(allele_table[i]), row_j(allele_table[j]);
                    //Only care about the different contigs.
                    hmr_remove_common(row_i, row_j);
                    //We want to find the maximum paired edges between these two rows.
                    //And only one contig is matching the other.
                    std::list<HMR_EDGE_INFO> edge_list;
                    for(const int32_t contig_i: row_i)
                    {
                        auto iter_i = edge_pair_map.find(contig_i);
                        if(iter_i == edge_pair_map.end())
                        {
                            continue;
                        }
                        NODE_COUNT_MAP &contig_i_map = iter_i->second;
                        for(const int32_t contig_j: row_j)
                        {
                            auto iter_j = contig_i_map.find(contig_j);
                            if(iter_j == contig_i_map.end())
                            {
                                continue;
                            }
                            edge_list.push_back(HMR_EDGE_INFO{contig_i, contig_j, static_cast<uint64_t>(iter_j->second), 0.0});
                        }
                    }
                    //If there is no edges for these two rows, we don't care about that.
                    if(edge_list.empty())
                    {
                        continue;
                    }
                    std::vector<HMR_EDGE_INFO> edges;
                    hMoveListToVector(edge_list, edges);
                    std::sort(edges.begin(), edges.end(),
                              [](const HMR_EDGE_INFO &lhs, const HMR_EDGE_INFO &rhs)
                    {
                        return lhs.pairs > rhs.pairs;
                    });
                    //Pick out the first edges, check whether these two contigs are still in the row.
                    //If not, remove the edge.
                    for(const HMR_EDGE_INFO &edge: edges)
                    {
                        if(!hmr_in_ordered_vector(edge.start, row_i) ||
                                !hmr_in_ordered_vector(edge.end, row_j))
                        {
                            //Remove the edge.
                            draft_mappings_remove_edge(edge_pair_map, edge.start, edge.end);
                            continue;
                        }
                        //Or else, we found an edge that is valid, remove the contig from the vector.
                        hmr_remove_one(row_i, edge.start);
                        hmr_remove_one(row_j, edge.end);
                    }
                }
            }
            //Convert the map to edge pairs.
            for(auto node_map_iter: edge_pair_map)
            {
                int32_t node_start = node_map_iter.first;
                for(auto end_iter: node_map_iter.second)
                {
                    HMR_EDGE edge = hmr_graph_edge(node_start, end_iter.first);
                    //Only add the edge that is not existed before.
                    if(edge_pairs.find(edge.data) == edge_pairs.end())
                    {
                        edge_pairs.insert(std::make_pair(edge.data, end_iter.second));
                    }
                }
            }
            time_print("%zu edges remain.", edge_pairs.size());
        }
        else
        {
            //Directly build the edge counter.
            hmr_graph_load_reads(opts.reads, opts.read_buffer_size, draft_mappings_build_edges, &edge_pairs);
            time_print("Paired-reads have been counted, %zu edge(s) generated.", edge_pairs.size());
        }
        //Filter out the valid links.
        time_print("Removing edges failed to reach the minimum links...");
        time_print("Calculating the node repetitive factors...");
        edges.reserve(edge_pairs.size());
        node_factors.resize(num_of_contigs);
        for (int32_t i = 0; i < num_of_contigs; ++i)
        {
            node_factors[i] = 0.0;
        }
        double links_average = 0.0;
        int64_t max_re_square = hSquare(static_cast<int64_t>(max_enzyme_count));
        int32_t node_start, node_end;
        for (const auto& iter : edge_pairs)
        {
            //Create the edge, count to the factors.
            hmr_graph_edge_pos(iter.first, node_start, node_end);
            //Skip the edge.
            int32_t edge_num_of_links = iter.second;
            if (edge_num_of_links < opts.min_links)
            {
                continue;
            }

            // Calculat the edge weights..
            double edge_weights = static_cast<double>(max_re_square * edge_num_of_links / nodes[node_start].enzyme_count / nodes[node_end].enzyme_count);
            edges.emplace_back(HMR_EDGE_INFO{ node_start, node_end, static_cast<uint64_t>(edge_num_of_links), edge_weights });
            edges.emplace_back(HMR_EDGE_INFO{ node_end, node_start, static_cast<uint64_t>(edge_num_of_links), edge_weights });
            node_factors[node_start] += edge_weights;
            node_factors[node_end] += edge_weights;
            links_average += edge_weights;
        }
        //Calculate the links average.
        links_average = 2.0 * links_average / static_cast<double>(num_of_contigs);
        time_print("%zu edge(s) generated, average links = %.2lf", edges.size(), links_average);
        time_print("Skip contig(s) with maximum multiplicity...");
        //Calculate the node factors.
        int32_t invalid_counter = 0;
        for (int32_t i = 0; i < num_of_contigs; ++i)
        {
            //Calculate the node factors.
            node_factors[i] /= links_average;
            //Check whether the node factors is valid or not.
            contig_invalid[i] |= (node_factors[i] >= opts.max_density);
            if (contig_invalid[i])
            {
                ++invalid_counter;
            }
        }
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
        if (!invalid_nodes.empty())
        {
            //Write the skipped nodes.
            std::string invalid_node_path = hmr_graph_path_contigs_invalid(opts.output);
            time_print("Saving invalid contig ids to %s", invalid_node_path.c_str());
            hmr_graph_save_contig_ids(invalid_node_path.c_str(), invalid_nodes);
            time_print("Done");
        }
    }
    //Based on the node factors, change the graph to directed-graph.
    time_print("Adjust link densities by repetitive factors...");
    for (size_t i = 0; i < edges.size(); ++i)
    {
        auto& edge_info = edges[i];
        edge_info.weights = ceil(edge_info.weights / node_factors[edge_info.start]);
    }
    time_print("%zu edge(s) generated.", edges.size());
    //Write the edge information to the target file.
    std::string edge_path = hmr_graph_path_edge(opts.output);
    time_print("Saving edge information to %s", edge_path.c_str());
    hmr_graph_save_edges(edge_path.c_str(), edges);
    time_print("Draft complete.");
    return 0;
}
