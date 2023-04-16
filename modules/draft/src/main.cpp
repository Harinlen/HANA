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
            //Build the node-to-node map.
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
                //Find out the best edge of the connected neighbours.
                std::unordered_map<int32_t, ALLELE_NEIGHBOUR> best_neighbors;
                //Collect the best matching node id.
                for(size_t i=0; i<conflict_row.size(); ++i)
                {
                    int32_t contig_id = conflict_row[i];
                    auto node_map_iter = edge_pair_map.find(contig_id);
                    if(node_map_iter == edge_pair_map.end())
                    {
                        continue;
                    }
                    std::vector<int32_t> contig_neighbours;
                    NODE_COUNT_MAP &node_map = node_map_iter->second;
                    contig_neighbours.reserve(node_map.size());
                    for(const auto iter: node_map)
                    {
                        contig_neighbours.push_back(iter.first);
                    }
                    //Check all the neighbour of the current node.
                    for(const auto neighbour_id: contig_neighbours)
                    {
                        //Try to find the node id in neighour counts.
                        auto current_iter = node_map_iter->second.find(neighbour_id);
                        if(current_iter == node_map_iter->second.end())
                        {
                            continue;
                        }
                        auto existed_count_iter = best_neighbors.find(neighbour_id);
                        if(existed_count_iter == best_neighbors.end())
                        {
                            //Create a new record for the contig.
                            best_neighbors.insert(std::make_pair(neighbour_id, ALLELE_NEIGHBOUR { contig_id, current_iter->second }));
                        }
                        else
                        {
                            //Need to compare with the row ID's record.
                            auto &existed_count_info = existed_count_iter->second;
                            if(existed_count_info.count > current_iter->second)
                            {
                                //Current edge need to be removed.
                                draft_mappings_remove_edge(edge_pair_map, neighbour_id, contig_id);
                            }
                            else
                            {
                                //Remove the neighour with its original connected contig.
                                draft_mappings_remove_edge(edge_pair_map, neighbour_id, existed_count_info.connected_id);
                                //Update the connected id and counts to current contig.
                                existed_count_info.connected_id = contig_id;
                                existed_count_info.count = current_iter->second;
                            }
                        }
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
            //Skip the edge.
            int32_t edge_num_of_links = iter.second;
            if (edge_num_of_links < opts.min_links)
            {
                continue;
            }
            //Create the edge, count to the factors.
            hmr_graph_edge_pos(iter.first, node_start, node_end);
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
