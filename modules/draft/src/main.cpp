#include <cstdlib>
#include <cstdio>
#include <cassert>

#include "hmr_args.hpp"
#include "hmr_path.hpp"
#include "hmr_ui.hpp"
#include "hmr_enzyme.hpp"
#include "hmr_fasta.hpp"
#include "hmr_mapping.hpp"
#include "hmr_bin_file.hpp"
#include "hmr_text_file.hpp"
#include "hmr_global.hpp"

#include "hmr_contig_graph_type.hpp"
#include "hmr_contig_graph.hpp"

#include "args_draft.hpp"
#include "draft_fasta_index.hpp"
#include "draft_mapping.hpp"

extern HMR_ARGS opts;

int main(int argc, char* argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.fasta) { help_exit(-1, "Missing FASTA file path."); }
    if (!path_can_read(opts.fasta)) { time_error(-1, "Cannot read FASTA file %s", opts.fasta); }
    if (opts.mappings.empty()) { help_exit(-1, "Missing Hi-C mapping file path."); }
    if (!opts.enzyme) { help_exit(-1, "Missing restriction enzyme cutting site."); }
    if (!opts.output) { help_exit(-1, "Missing output file prefix."); }
    //Convert the enzyme into sequence.
    hmr_enzyme_formalize(opts.enzyme, &opts.enzyme_nuc, &opts.enzyme_nuc_length);
    time_print("Execution configuration:");
    time_print("\tMinimum map quality: %d", opts.mapq);
    time_print("\tRestriction enzyme: %s", opts.enzyme_nuc);
    time_print("\tMinimum restriction enzyme count: %d", opts.min_enzymes);
    time_print("\tHalf of enzyme range: %d", opts.range);
    time_print("\tMapping file count: %zu", opts.mappings.size());
    time_print("\tThreads: %d", opts.threads);
    //Load the FASTA and find the enzyme ranges in sequences.
    HMR_CONTIGS contigs;
    HMR_CONTIG_INVALID_SET invalid_id_set;
    ENZYME_RANGES* contig_ranges;
    size_t max_enzyme_counts = 0;
    {
        //Prepare the enzyme for searching.
        ENZYME_SEARCH search;
        contig_draft_search_start(opts.enzyme_nuc, opts.enzyme_nuc_length, search);
        //Read the FASTA and build the enzyme ranges.
        DRAFT_NODES_USER node_user{ opts.range, &contigs, &search, NULL, NULL, NULL };
        {
            //Read the FASTA file and build the enzyme index.
            RANGE_SEARCH_POOL search_pool(contig_range_search, opts.threads * 32, opts.threads);
            node_user.pool = &search_pool;
            time_print("Searching enzyme in %s", opts.fasta);
            hmr_fasta_read(opts.fasta, contig_draft_build, &node_user);
        }
        contig_draft_search_end(search);
        //Convert the search node information.
        contig_ranges = static_cast<ENZYME_RANGES*>(malloc(sizeof(ENZYME_RANGES) * contigs.size()));
        ENZYME_RANGE_CHAIN* chain_node = node_user.chain_head;
        int32_t range_index = 0;
        while (chain_node != NULL)
        {
            //Calculate the max enzyme count.
            max_enzyme_counts = hMax(max_enzyme_counts, chain_node->data.counter);
            //Save the data to contig range.
            contig_ranges[range_index].counter = chain_node->data.counter;
            contig_ranges[range_index] = chain_node->data;
            ++range_index;
            //Free the chain node.
            ENZYME_RANGE_CHAIN* next = chain_node->next;
            free(chain_node);
            chain_node = next;
        }
        if (range_index != contigs.size())
        {
            time_error(-1, "Unexpected error: contig size and range index size mismatch.");
        }
        //Search complete.
        time_print("%zu contig(s) indexed.", range_index);
        //Dump the node data to target file.
        {
            std::string path_contig = hmr_graph_path_contig(opts.output);
            time_print("Save contig information to %s", path_contig.data());
            hmr_graph_save_contigs(path_contig.data(), contigs);
            time_print("Done");
        }
        //Checking which contig is valid, if invalid, generate the invalid list.
        time_print("Checking invalid contig(s)...");
        HMR_CONTIG_INVALID_IDS invalid_ids;
        for (int32_t i = 0; i < range_index; ++i)
        {
            if (contig_ranges[i].counter < opts.min_enzymes)
            {
                invalid_ids.push_back(i);
            }
        }
        time_print("%zu invalid contig(s) detected.", invalid_ids.size());
        if (!invalid_ids.empty())
        {
            std::string path_invalid = hmr_graph_path_invalid(opts.output);
            time_print("Save invalid contig indices to %s", path_invalid.data());
            hmr_graph_save_invalid(path_invalid.data(), invalid_ids);
            time_print("Done");
            //Construct the invalid set.
            invalid_id_set = HMR_CONTIG_INVALID_SET(invalid_ids.begin(), invalid_ids.end());
        }
    }
    //Build the contig edge counters.
    HMR_EDGE_COUNTERS edge_counters;
    double max_enzyme_counts_square = static_cast<double>(hSquare(max_enzyme_counts));
    {
        time_print("Constructing contig name -> id map...");
        //Map contig name to an index.
        CONTIG_ID_MAP contig_ids;
        for (size_t i = 0; i < contigs.size(); ++i)
        {
            contig_ids.insert(std::make_pair(std::string(contigs[i].name, contigs[i].name_size), static_cast<int>(i)));
        }
        time_print("Contig name -> id map has been built.");
        //Prepare the read-pair information output.
        std::string path_reads = hmr_graph_path_reads(opts.output);
        FILE* reads_file = NULL;
        if (!bin_open(path_reads.data(), &reads_file, "wb"))
        {
            time_error(-1, "Failed to create read information file %s", path_reads.data());
        }
        time_print("Writing reads summary information to %s", path_reads.data());
        //Prepare the output buffer.
        size_t output_size = sizeof(HMR_MAPPING) * 1024;
        char* output_buffer = NULL;
        if (output_size)
        {
            output_buffer = static_cast<char*>(malloc(output_size));
            assert(output_buffer);
        }
        //Read all the mapping files and count the edges.
        MAPPING_DRAFT_USER mapping_user { contig_ids, invalid_id_set, contig_ranges, static_cast<uint8_t>(opts.mapq), NULL, 0, reads_file, output_buffer, 0, output_size, RAW_EDGE_MAP()};
        for (char* mapping_path : opts.mappings)
        {
            time_print("Loading reads from %s", mapping_path);
            //Build the reads mapping.
            hmr_mapping_read(mapping_path, MAPPING_PROC{ mapping_draft_n_contig, mapping_draft_contig, mapping_draft_read_align }, &mapping_user, opts.threads);
            //Recover the mapping array.
            delete[] mapping_user.contig_id_map;
            //Reset the contig map.
            mapping_user.contig_id_map = NULL;
            mapping_user.contig_idx = 0;
        }
        time_print("Contig edges built from %zu file(s).", opts.mappings.size());
        //Check reads file buffer is complete.
        if (mapping_user.output_offset)
        {
            fwrite(mapping_user.output_buffer, mapping_user.output_offset, 1, reads_file);
        }
        fclose(reads_file);
        free(mapping_user.output_buffer);
        time_print("Reads summary information saved.");
        //Build the edge map.
        time_print("Calculating %zu edge weights...", mapping_user.edges.size());
        edge_counters.reserve(mapping_user.edges.size());
        for (const auto& edge_info : mapping_user.edges)
        {
            //Construct the edge counter.
            HMR_EDGE edge;
            edge.data = edge_info.first;
            int32_t pair_count = edge_info.second;
            double weight = max_enzyme_counts_square / static_cast<double>(contig_ranges[edge.pos.start].counter) / static_cast<double>(contig_ranges[edge.pos.end].counter) * pair_count;
            edge_counters.push_back(HMR_EDGE_INFO{ edge, pair_count, weight });
        }
        time_print("Done");
    }
    //Dump the edge information into files.
    std::string path_edge = hmr_graph_path_edge(opts.output);
    time_print("Save contig edge information to %s", path_edge.data());
    hmr_graph_save_edge(path_edge.data(), edge_counters);
    time_print("Done");
    return 0;
}
