#include <algorithm>
#include <cassert>
#include <cstring>

#include "hmr_args.hpp"
#include "hmr_bin_file.hpp"
#include "hmr_contig_graph.hpp"
#include "hmr_enzyme.hpp"
#include "hmr_fasta.hpp"
#include "hmr_global.hpp"
#include "hmr_path.hpp"
#include "hmr_ui.hpp"

#include "args_extract.hpp"
#include "extract_fasta.hpp"
#include "extract_mapping.hpp"
#include "extract_allele.hpp"

extern HMR_ARGS opts;

int main(int argc, char* argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.fasta) { help_exit(-1, "Missing FASTA file path."); }
    if (!path_can_read(opts.fasta)) { time_error(-1, "Cannot read FASTA file %s", opts.fasta); }
    if (opts.allele && !path_can_read(opts.allele)) { time_error(-1, "Cannot read allele table %s", opts.allele); }
    for (const char* mapping_path : opts.mappings)
    {
        if (!path_can_read(mapping_path))
        {
            time_error(-1, "Cannot read mapping file %s", mapping_path);
        }
    }
    if (!opts.output) { help_exit(-1, "Missing output file prefix."); }
    //Check enzyme validation.
    if (opts.enzyme.empty()) { help_exit(-1, "Missing restriction enzyme cutting site(s)."); }
    ENZYME_VEC enzymes = hmr_enzyme_formalize(opts.enzyme), weighted_enzymes;
    if (!opts.weight_enzyme.empty())
    {
        //Formalize the weighted enzyme.
        weighted_enzymes = hmr_enzyme_formalize(opts.weight_enzyme);
    }
    time_print("Execution configuration:");
    time_print("\tBAM minimum map quality: %d", opts.mapq);
    {
        std::string enzyme_total = enzymes[0];
        for (size_t i=1; i<enzymes.size(); ++i)
        {
            enzyme_total += ", " + std::string(enzymes[i]);
        }
        time_print("\tRestriction enzyme(s): %s", enzyme_total.c_str());
    }
    if(!weighted_enzymes.empty())
    {
        std::string enzyme_total = weighted_enzymes[0];
        for (size_t i=1; i<weighted_enzymes.size(); ++i)
        {
            enzyme_total += ", " + std::string(weighted_enzymes[i]);
        }
        time_print("\tWeight calculated restriction enzyme(s): %s", enzyme_total.c_str());
    }
    time_print("\tValid enzyme distance of Hi-C pairs: %dx2", opts.range);
    time_print("\tRead length used for pairs file: %d", opts.pairs_read_len);
    time_print("\tMapping file(s): %zu", opts.mappings.size());
    time_print("\tThreads: %d", opts.threads);
    time_print("\tFASTA search buffer records size / thread: %d", opts.fasta_pool);
    time_print("\tMapping filter buffer records size / thread: %dK", opts.mapping_pool);
    opts.mapping_pool <<= 10;
    //Construct the check flag.
    uint16_t check_flag = 0xFFFF;
    time_print("\tChecking flag settings...");
    if (opts.skip_flag) { check_flag &= ~CHECK_FLAG_FLAG; time_print("\tSkip FLAG checking."); }
    if (opts.skip_range) { check_flag &= ~CHECK_FLAG_RANGE; time_print("\tSkip range checking."); }
    //Load the FASTA and find the enzyme ranges in sequences.
    HMR_CONTIGS nodes;
    CONTIG_ENZYME_RANGES contig_enzyme_ranges;
    {
        //Prepare the enzyme for searching.
        CANDIDATE_ENZYMES search_range, search_calc;
        extract_enzyme_search_start(enzymes, search_range);
        if(!weighted_enzymes.empty())
        {
            extract_enzyme_search_start(weighted_enzymes, search_calc);
        }
        CONTIG_CHAIN node_chain;
        CONTIG_NAME_CHAIN node_name_chain;
        CONTIG_RANGE_RESULTS node_ranges;
        {
            //Start the enzyme search pool.
            RANGE_SEARCH_POOL pool(contig_range_search, opts.threads * opts.fasta_pool, opts.threads);
            //Construct the contig info build user.
            EXTRACT_FASTA_USER node_build_user {search_range, search_calc, node_chain, node_name_chain, pool, opts.range, node_ranges };
            time_print("Searching enzyme in %s", opts.fasta);
            hmr_fasta_read(opts.fasta, extract_fasta_search_proc, &node_build_user);
        }
        extract_enzyme_search_end(search_range);
        //Convert the node chain into node vector.
        hDequeListToVector(node_name_chain, nodes.names);
        hDequeListToVector(node_chain, nodes.contigs);
        //Construct the enzyme range vector and find invalid ids.
        contig_enzyme_ranges.resize(node_ranges.size());
        while (!node_ranges.empty())
        {
            const auto& node_range = node_ranges.front();
            int32_t node_id = node_range.contig_index;
            //Increase 1 to avoid divided by zero.
            nodes.contigs[node_id].enzyme_count = node_range.counter + 1;
            contig_enzyme_ranges[node_id] = node_range.ranges;
            node_ranges.pop_front();
        }
    }
    time_print("%zu contig(s) are indexed.", nodes.contigs.size());
    //Dump the node data to target file.
    {
        std::string path_contig = hmr_graph_path_contigs(opts.output);
        time_print("Save contig information to %s", path_contig.data());
        hmr_graph_save_contigs(path_contig.data(), nodes);
        time_print("Done");
    }
    time_print("Constructing contig index map...");
    CONTIG_INDEX_MAP contig_index_map;
    for (size_t i = 0; i < nodes.names.size(); ++i)
    {
        contig_index_map.insert(std::make_pair(std::string(nodes.names[i].name, nodes.names[i].name_size), static_cast<int>(i)));
    }
    time_print("Contig index map has been built.");
    if(opts.allele)
    {
        time_print("Building allele table from %s", opts.allele);
        HMR_CONTIG_ID_TABLE allele_table = extract_allele_table(opts.allele, &contig_index_map);
        std::string path_allele_table = hmr_graph_path_allele_table(opts.output);
        time_print("Saving allele table to %s", path_allele_table.data());
        hmr_graph_save_contig_table(path_allele_table.data(), allele_table);
        time_print("Done");
    }
    //Prepare the read-pair information output.
    std::string path_reads = hmr_graph_path_reads(opts.output);
    FILE* reads_file = NULL;
    if (!bin_open(path_reads.data(), &reads_file, "wb"))
    {
        time_error(-1, "Failed to create read information file %s", path_reads.data());
    }
    time_print("Writing reads information to %s", path_reads.data());
    for (char* mapping_path : opts.mappings)
    {
        time_print("Loading reads from %s", mapping_path);
        extract_mapping_file(mapping_path, &contig_index_map, reads_file, &contig_enzyme_ranges, check_flag, opts.pairs_read_len, opts.mapq, opts.mapping_pool, opts.threads);
    }
    fclose(reads_file);
    time_print("Extract complete.");
    return 0;
}
