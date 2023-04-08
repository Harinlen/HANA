#include <cassert>

#include "hmr_contig_graph.hpp"
#include "hmr_global.hpp"
#include "hmr_ui.hpp"

#include "draft_mapping.hpp"

bool position_in_range(int32_t pos, const ENZYME_RANGES& ranges)
{
    for (size_t i = 0; i < ranges.length; ++i)
    {
        const ENZYME_RANGE& r = ranges.ranges[i];
        //Check whether the position is in the range.
        if (r.start <= pos && pos <= r.end)
        {
            return true;
        }
    }
    //No position matched.
    return false;
}

inline int32_t contig_id(const MAPPING_CONTIG_MAP &map, int32_t refId)
{
    return (refId > -1 && refId < map.contig_idx) ? map.contig_id_map[refId] : -1;
}

void mapping_draft_filter_worker(MAPPING_FILTER_USER& filter, MAPPING_FILTER_WORKER& worker, const MAPPING_CONTIG_MAP& map)
{
    std::unique_lock<std::mutex> lock(worker.start_mutex);
    auto& edge_map = worker.edges;
    while (true)
    {
        //Wait for worker start.
        worker.start_cv.wait(lock, [&] { return worker.start; });
        if (filter.exit)
        {
            break;
        }
        //Loop of my area.
        for (size_t i = worker.buf_start; i < worker.buf_end; ++i)
        {
            MAPPING_INFO& mapping_info = filter.proc_buf[i];
            //Check whether the mapping info is valid, then check whether the position is in range.
            int32_t ref_index = contig_id(map, mapping_info.refID),
                next_ref_index = contig_id(map, mapping_info.next_refID);
            // 3852 stands for:
            // - read unmapped (0x4)
            // - mate unmapped (0x8)
            // - not primary alignment (0x100)
            // - read fails platform/vendor quality checks (0x200)
            // - read is PCR or optical duplicate (0x400)
            // - supplementary alignment (0x800)
            if (ref_index == -1 || next_ref_index == -1 //Reference index invalid detection.
                || mapping_info.mapq == 255 //Map quality is invalid.
                || mapping_info.mapq < filter.mapq //Check whether the mapping reaches the minimum quality
                || (mapping_info.flag & 3852) != 0 // Filtered source code.
                || ref_index == next_ref_index // We don't care about the pairs on the same contigs.
                || (!position_in_range(mapping_info.pos, filter.contig_ranges[ref_index]))) // Or the position is not in the position.
            {
                continue;
            }
            //Output the edge to the output buffer.
            if (worker.mapping_buf_size == worker.mapping_buf_offset)
            {
                {
                    std::unique_lock<std::mutex> lock(filter.reads_mutex);
                    //Flush the entire buffer to the file.
                    fwrite(worker.mapping_buf, worker.mapping_buf_size, 1, filter.reads_file);
                }
                //Reset the offset back to beginning.
                worker.mapping_buf_offset = 0;
            }
            //Save the edges to read buffer.
            HMR_MAPPING* output_read = reinterpret_cast<HMR_MAPPING*>(worker.mapping_buf + worker.mapping_buf_offset);
            (*output_read) = HMR_MAPPING{ ref_index, mapping_info.pos, next_ref_index, mapping_info.next_pos };
            worker.mapping_buf_offset += sizeof(HMR_MAPPING);
            //Construct the edge info.
            HMR_EDGE edge = hmr_graph_edge(ref_index, next_ref_index);
            //A new pair is found.
            auto edge_iter = edge_map.find(edge.data);
            if (edge_iter == edge_map.end())
            {
                //Insert one count to the data.
                edge_map.insert(std::make_pair(edge.data, 1));
            }
            else
            {
                ++(edge_iter->second);
            }
        }
        //Disable the start state.
        worker.start = false;
        //Increase the finish counter.
        {
            filter.counter_mutex.lock();
            ++filter.finished_counter;
            if (filter.finished_counter == filter.threads)
            {
                filter.is_working = false;
                filter.finished_cv.notify_all();
            }
            filter.counter_mutex.unlock();
        }
    }
}

void mapping_draft_n_contig(uint32_t n_ref, void* user)
{
    MAPPING_DRAFT_USER* mapping_user = reinterpret_cast<MAPPING_DRAFT_USER*>(user);
    //Allocate the contig id maps.
    mapping_user->map.contig_id_map = new int32_t[n_ref];
    mapping_user->map.contig_idx = 0;
}

void mapping_draft_contig(uint32_t name_length, char* name, uint32_t, void* user)
{
    MAPPING_DRAFT_USER* mapping_user = reinterpret_cast<MAPPING_DRAFT_USER*>(user);
    //Build the contig mapping index.
    std::string contig_name(name, name_length);
    auto name_finder = mapping_user->contig_ids.find(contig_name);
    int32_t contig_id = (name_finder == mapping_user->contig_ids.end()) ? -1 : name_finder->second;
    //Check the contig id.
    if (contig_id != -1)
    {
        //If the the contig id is in invalid set, then set the id to be -1.
        if (mapping_user->invalid_ids.find(contig_id) != mapping_user->invalid_ids.end())
        {
            contig_id = -1;
        }
    }
    //Record the contig id at the map.
    mapping_user->map.contig_id_map[mapping_user->map.contig_idx] = contig_id;
    ++mapping_user->map.contig_idx;
}

void mapping_draft_read_align(size_t, const MAPPING_INFO& mapping_info, void* user)
{
    MAPPING_DRAFT_USER* mapping_user = reinterpret_cast<MAPPING_DRAFT_USER*>(user);
    //Directly save the pair info to the buffer.
    MAPPING_FILTER_USER &filter = mapping_user->filter;
    if (*filter.build_size < filter.buf_max_size)
    {
        filter.build_buf[*filter.build_size] = mapping_info;
        (*filter.build_size)++;
        return;
    }
    //Check whether the threads are running.
    if (filter.is_working)
    {
        //Wait for the filter to be complete.
        std::unique_lock<std::mutex> lock(filter.finished_mutex);
        filter.finished_cv.wait(lock, [&] { return !filter.is_working; });
    }
    //Swap the build and current buf.
    if (filter.activated_0)
    {
        filter.activated_0 = false;
        filter.build_buf = filter.buf_1;
        filter.build_size = &filter.buf_1_size;
        filter.proc_buf = filter.buf_0;
        filter.proc_size = &filter.buf_0_size;
    }
    else
    {
        filter.activated_0 = true;
        filter.build_buf = filter.buf_0;
        filter.build_size = &filter.buf_0_size;
        filter.proc_buf = filter.buf_1;
        filter.proc_size = &filter.buf_1_size;
    }
    //Reset the build size.
    *filter.build_size = 0;
    filter.finished_counter = 0;
    //Start the threads.
    MAPPING_FILTER_WORKER* workers = mapping_user->workers;
    //Start all the workers.
    filter.is_working = true;
    for (int32_t i = 0; i < filter.threads; ++i)
    {
        //Start the worker.
        workers[i].start = true;
        workers[i].start_cv.notify_one();
    }
}
