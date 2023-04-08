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
    //Ensure the threads is an even number.
    opts.threads = ((opts.threads + 1) >> 1) << 1;
    //Convert the enzyme into sequence.
    hmr_enzyme_formalize(opts.enzyme, &opts.enzyme_nuc, &opts.enzyme_nuc_length);
    time_print("Execution configuration:");
    time_print("\tMinimum map quality: %d", opts.mapq);
    time_print("\tRestriction enzyme: %s", opts.enzyme_nuc);
    time_print("\tMinimum restriction enzyme count: %d", opts.min_enzymes);
    time_print("\tHalf of enzyme range: %d", opts.range);
    time_print("\tMapping file count: %zu", opts.mappings.size());
    time_print("\tThreads: %d", opts.threads);
    time_print("\tFASTA search buffer per thread: %d", opts.fasta_pool);
    time_print("\tMapping cache buffer size: %d M", opts.mapping_pool);
    opts.mapping_pool <<= 13;
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
            RANGE_SEARCH_POOL search_pool(contig_range_search, opts.threads * opts.fasta_pool, opts.threads);
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
        if (static_cast<size_t>(range_index) != contigs.size())
        {
            time_error(-1, "Unexpected error: contig size and range index size mismatch.");
        }
        //Search complete.
        time_print("%zu contig(s) indexed.", range_index);
        //Checking which contig is valid, if invalid, generate the invalid list.
        time_print("Checking invalid contig(s)...");
        HMR_CONTIG_INVALID_IDS invalid_ids;
        for (int32_t i = 0; i < range_index; ++i)
        {
            contigs[i].enzyme_count = static_cast<int32_t>(contig_ranges[i].counter);
            if (contig_ranges[i].counter < static_cast<size_t>(opts.min_enzymes))
            {
                invalid_ids.push_back(i);
            }
        }
        time_print("%zu invalid contig(s) detected.", invalid_ids.size());
        //Dump the node data to target file.
        {
            std::string path_contig = hmr_graph_path_contig(opts.output);
            time_print("Save contig information to %s", path_contig.data());
            hmr_graph_save_contigs(path_contig.data(), contigs);
            time_print("Done");
        }
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
    //Initialize filter user and mapping workers.
    int32_t num_of_workers = opts.threads >> 1;
    MAPPING_CONTIG_MAP contig_map{ NULL, 0 };
    MAPPING_FILTER_USER filter_user;
    filter_user.contig_ranges = contig_ranges;
    filter_user.mapq = static_cast<uint8_t>(opts.mapq);
    filter_user.threads = num_of_workers;
    filter_user.reads_file = reads_file;
    filter_user.buf_max_size = static_cast<size_t>(opts.mapping_pool * num_of_workers);
    filter_user.buf_0 = static_cast<MAPPING_INFO*>(malloc(sizeof(MAPPING_INFO) * filter_user.buf_max_size));
    filter_user.buf_1 = static_cast<MAPPING_INFO*>(malloc(sizeof(MAPPING_INFO) * filter_user.buf_max_size));
    if (!filter_user.buf_0 || !filter_user.buf_1)
    {
        time_error(-1, "Failed to create filter buffer.");
    }
    filter_user.activated_0 = true;
    filter_user.build_buf = filter_user.buf_0;
    filter_user.build_size = &filter_user.buf_0_size;
    //Prepare worker data.
    MAPPING_FILTER_WORKER* workers = new MAPPING_FILTER_WORKER[num_of_workers];
    if (!workers)
    {
        time_error(-1, "No enough memory for creating workers.");
    }
    assert(workers);
    size_t worker_start = 0;
    std::thread *worker_thread = new std::thread[num_of_workers];
    if (!worker_thread)
    {
        time_error(-1, "Failed to create worker threads.");
    }
    for (int32_t i = 0; i < num_of_workers; ++i)
    {
        //Clear the worker start.
        workers[i].start = false;
        //Configure the worker range.
        workers[i].buf_start = worker_start;
        worker_start += opts.mapping_pool;
        workers[i].buf_end = worker_start;
        //Initial the worker saving buffer.
        workers[i].mapping_buf_size = static_cast<size_t>(sizeof(HMR_MAPPING) * opts.mapping_pool);
        workers[i].mapping_buf = static_cast<char *>(malloc(workers[i].mapping_buf_size));
        if (!workers[i].mapping_buf)
        {
            time_error(-1, "No enough memory for the worker mapping buffer.");
        }
        workers[i].mapping_buf_offset = 0;
        //Start the worker thread.
        worker_thread[i] = std::thread(mapping_draft_filter_worker, std::ref(filter_user), std::ref(workers[i]), std::ref(contig_map));
    }
    //Prepare the draft user information.
    MAPPING_DRAFT_USER mapping_user{ contig_ids, invalid_id_set, filter_user, workers, contig_map };
    for (char* mapping_path : opts.mappings)
    {
        time_print("Loading reads from %s", mapping_path);
        //Build the reads mapping.
        hmr_mapping_read(mapping_path, MAPPING_PROC{ mapping_draft_n_contig, mapping_draft_contig, mapping_draft_read_align }, &mapping_user, num_of_workers);
        //Recover the mapping array.
        delete[] contig_map.contig_id_map;
        //Reset the contig map.
        contig_map.contig_id_map = NULL;
        contig_map.contig_idx = 0;
    }
    //Let the worker process the rest of the data.
    if ((*filter_user.build_size) != 0)
    {
        //Check whether the threads are running.
        if (filter_user.is_working)
        {
            //Wait for the filter to be complete.
            std::unique_lock<std::mutex> lock(filter_user.finished_mutex);
            filter_user.finished_cv.wait(lock, [&] { return !filter_user.is_working; });
        }
        //Start the last work.
        filter_user.proc_buf = filter_user.build_buf;
        filter_user.proc_size = filter_user.build_size;
        filter_user.is_working = true;
        filter_user.finished_counter = 0;
        //Update the working range.
        size_t worker_load = ((*filter_user.proc_size) + num_of_workers - 1) / num_of_workers;
        worker_start = 0;
        for (int32_t i = 0; i < num_of_workers; ++i)
        {
            //Start the worker.
            workers[i].buf_start = worker_start;
            worker_start += worker_load;
            workers[i].buf_end = hMin(worker_start, *filter_user.proc_size);
            //Start the current worker.
            workers[i].start = true;
            workers[i].start_cv.notify_one();
        }
        //Wait for worker complete their work.
        if (filter_user.is_working)
        {
            //Wait for the filter to be complete.
            std::unique_lock<std::mutex> lock(filter_user.finished_mutex);
            filter_user.finished_cv.wait(lock, [&] { return !filter_user.is_working; });
        }
    }
    //Exit all the workers.
    filter_user.exit = true;
    for (int32_t i = 0; i < num_of_workers; ++i)
    {
        workers[i].start = true;
        workers[i].start_cv.notify_one();
        //Wait for thread exit.
        worker_thread[i].join();
    }
    //Free the workers.
    delete[] worker_thread;
    free(filter_user.buf_0);
    free(filter_user.buf_1);
    //Reduce the worker result to a single variable.
    time_print("Contig edges built from %zu file(s).", opts.mappings.size());
    time_print("Merging edge map...");
    RAW_EDGE_MAP edges;
    for (int32_t i = 0; i < num_of_workers; ++i)
    {
        //Flush the entire buffer to the file.
        fwrite(workers[i].mapping_buf, workers[i].mapping_buf_offset, 1, reads_file);
        free(workers[i].mapping_buf);
        //Merge the edges together.
        for (const auto& edge_weight : workers[i].edges)
        {
            auto total_edges = edges.find(edge_weight.first);
            if (total_edges == edges.end())
            {
                edges.insert(edge_weight);
            }
            else
            {
                total_edges->second += edge_weight.second;
            }
        }
    }
    delete[] workers;
    fclose(reads_file);
    time_print("Edge data merged.");
    //Build the edge map.
    time_print("Calculating %zu edge weights...", edges.size());
    //Build the contig edge counters.
    HMR_EDGE_COUNTERS edge_counters;
    edge_counters.reserve(edges.size());
    double max_enzyme_counts_square = static_cast<double>(hSquare(max_enzyme_counts));
    for (const auto& edge_info : edges)
    {
        //Construct the edge counter.
        HMR_EDGE edge;
        edge.data = edge_info.first;
        int32_t pair_count = edge_info.second;
        double weight = max_enzyme_counts_square / static_cast<double>(contigs[edge.pos.start].enzyme_count) / static_cast<double>(contigs[edge.pos.end].enzyme_count) * pair_count;
        edge_counters.push_back(HMR_EDGE_INFO{ edge.pos.start, edge.pos.end, pair_count, 0, weight });
    }
    time_print("Done");
    //Dump the edge information into files.
    std::string path_edge = hmr_graph_path_edge(opts.output);
    time_print("Save contig edge information to %s", path_edge.data());
    hmr_graph_save_edge(path_edge.data(), edge_counters);
    time_print("Done");
    return 0;
}
;