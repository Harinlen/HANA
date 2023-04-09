#include <cstring>

#include "hmr_bgzf.hpp"
#include "hmr_bin_queue.hpp"
#include "hmr_ui.hpp"

#include "hmr_bam.hpp"

void hmr_bam_read(const char* filepath, BAM_MAPPING_PROC proc, void* user, int threads)
{
    //Open the .bam file as BGZF file.
    HMR_BGZF_HANDLER* bgzf_handler = hmr_bgzf_open(filepath, threads);
    //Fetch and check the magic number.
    auto buf = bgzf_handler->buffer;
    auto queue = bgzf_handler->queue;
    char* magic = hmr_bin_buf_fetch(buf, queue, 4);
    if (strncmp(magic, "BAM\1", 4))
    {
        time_error(1, "BAM header magic string incorrect.");
    }
    //Skip the header text.
    uint32_t l_text = hmr_bin_buf_fetch_uint32(buf, queue);
    char* text = hmr_bin_buf_fetch(buf, queue, l_text);
    //Fetch the n_ref.
    uint32_t n_ref = hmr_bin_buf_fetch_uint32(buf, queue);
    proc.proc_no_of_contig(n_ref, user);
    //Loop until all the reference are parsed.
    while (n_ref--)
    {
        //Format:
        // [name length] [name] [seq length]
        // name length include '\0'
        uint32_t l_name = hmr_bin_buf_fetch_uint32(buf, queue);
        char* name = hmr_bin_buf_fetch(buf, queue, l_name);
        uint32_t l_ref = hmr_bin_buf_fetch_uint32(buf, queue);
        proc.proc_contig(l_name - 1, name, l_ref, user);
    }
    //Fetch the rest of the data (align data).
    char* block_size_data = hmr_bin_buf_fetch(buf, queue, 4);
    size_t block_id = 0;
    while (block_size_data)
    {
        uint32_t block_size = *(reinterpret_cast<uint32_t*>(block_size_data));
        //Fetch the data of the fetch.
        char* block_data = hmr_bin_buf_fetch(buf, queue, block_size);
        //Call the process function.
        proc.proc_read_align(block_id, reinterpret_cast<BAM_BLOCK_HEADER*>(block_data), user);
        //Increase the block id.
        ++block_id;
        //Fetch the next block.
        block_size_data = hmr_bin_buf_fetch(buf, queue, 4);
    }
    //Close the BGZF file.
    hmr_bgzf_close(bgzf_handler);
}
