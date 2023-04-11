#include <cassert>
#include <cstring>
#include <thread>
#include <mutex>

#include "hmr_algorithm.hpp"
#include "hmr_bin_file.hpp"
#include "hmr_global.hpp"
#include "hmr_ui.hpp"
#include "hmr_path.hpp"

#include "hmr_contig_graph.hpp"

template<typename T>
struct GRAPH_LOAD_BUFFER
{
    T* buf;
    int32_t buf_size;
};

template<typename T>
struct GRAPH_BUF_LOADER
{
    GRAPH_LOAD_BUFFER<T> buf_0, buf_1;
    GRAPH_LOAD_BUFFER<T>* buf_loading, * buf_processing;
    int32_t buf_size;
    bool is_loaded, finished;
    std::mutex wait_processing_mutex, wait_loaded_mutex;
    std::condition_variable wait_processing_cv, wait_loaded_cv;
};

template<typename T>
void hmr_graph_buffer_init(GRAPH_LOAD_BUFFER<T>& buffer, int32_t buf_size)
{
    buffer.buf = static_cast<T *>(malloc(sizeof(T) * buf_size));
    if (!buffer.buf)
    {
        time_error(-1, "Failed to allocate memory for reads buffer.");
    }
    assert(buffer.buf);
    buffer.buf_size = buf_size;
}

template<typename T>
void hmr_graph_buffer_free(GRAPH_LOAD_BUFFER<T>& buffer)
{
    free(buffer.buf);
}

template<typename T>
void hmr_graph_buffer_loader(const char* filepath, GRAPH_BUF_LOADER<T>& loader, void(*size_proc)(uint64_t, void*), void *user)
{
    FILE* data_file;
    if (!bin_open(filepath, &data_file, "rb"))
    {
        time_error(-1, "Failed to open buffered data file %s", filepath);
    }
    //Check size loader.
    if (size_proc)
    {
        //Read the size from the file.
        uint64_t item_size;
        fread(&item_size, sizeof(uint64_t), 1, data_file);
        size_proc(item_size, user);
    }
    //Loop until all the data is loaded.
    std::unique_lock<std::mutex> processing_lock(loader.wait_processing_mutex);
    while (!loader.finished)
    {
        //Start to load data.
        loader.buf_loading->buf_size = static_cast<int32_t>(fread(loader.buf_loading->buf, sizeof(T), loader.buf_size, data_file));
        //Check whether we can still fill all the part of the buffer.
        loader.finished = loader.buf_loading->buf_size < loader.buf_size;
        //Set as loaded.
        loader.is_loaded = true;
        //Notify the main thread.
        loader.wait_loaded_cv.notify_one();
        //Wait for loaded flag to be off.
        if (!loader.finished)
        {
            loader.wait_processing_cv.wait(processing_lock, [&] {return !loader.is_loaded; });
        }
    }
    fclose(data_file);
}

template<typename T>
void hmr_graph_load_with_buffer(const char* filepath, int32_t buf_size,
    void(*size_proc)(uint64_t, void*),
    void(*proc)(T*, int32_t, void*), 
    void* user)
{
    //Allocate buffer for processing.
    GRAPH_BUF_LOADER<T> loader;
    hmr_graph_buffer_init(loader.buf_0, buf_size);
    hmr_graph_buffer_init(loader.buf_1, buf_size);
    loader.buf_loading = &loader.buf_0;
    loader.buf_processing = &loader.buf_1;
    loader.buf_size = buf_size;
    loader.is_loaded = false;
    loader.finished = false;
    //Loop for all the data is processed.
    std::thread loader_thread(hmr_graph_buffer_loader<T>, filepath, std::ref(loader), size_proc, user);
    while (!loader.finished)
    {
        //Wait for the loader loading a chunk.
        std::unique_lock<std::mutex> loaded_lock(loader.wait_loaded_mutex);
        loader.wait_loaded_cv.wait(loaded_lock, [&] {return loader.is_loaded; });
        //One chunk is loaded, we flip the buffer, let the thread to loaded next part.
        auto buf_temp = loader.buf_loading;
        loader.buf_loading = loader.buf_processing;
        loader.buf_processing = buf_temp;
        //Reset the loaded flag.
        loader.is_loaded = false;
        loader.wait_processing_cv.notify_one();
        //Processing the buffer.
        proc(loader.buf_processing->buf, loader.buf_processing->buf_size, user);
    }
    //Process the last part.
    if (loader.is_loaded)
    {
        proc(loader.buf_loading->buf, loader.buf_loading->buf_size, user);
    }
    //Close the thread.
    loader_thread.join();
    //Free the buffer.
    hmr_graph_buffer_free(loader.buf_0);
    hmr_graph_buffer_free(loader.buf_1);
}

std::string hmr_graph_path_contigs(const char* prefix)
{
    return std::string(prefix) + ".hmr_nodes";
}

std::string hmr_graph_path_reads(const char* prefix)
{
    return std::string(prefix) + ".hmr_reads";
}

std::string hmr_graph_path_nodes_invalid(const char* contig_path)
{
    return std::string(contig_path) + "_invalid";
}

std::string hmr_graph_path_edge(const char* prefix)
{
    return std::string(prefix) + ".hmr_edges";
}

std::string hmr_graph_path_contigs_invalid(const char* prefix)
{
    return std::string(prefix) + ".hmr_nodes_invalid";
}

std::string hmr_graph_path_allele_table(const char* prefix)
{
    return std::string(prefix) + ".hmr_allele_table";
}

std::string hmr_graph_path_cluster_name(const char* prefix, const int32_t index, const int32_t total)
{
    return std::string(prefix) + "_" + std::to_string(index) + "g" + std::to_string(total) + ".hmr_group";
}

void hmr_graph_load_contigs(const char* filepath, HMR_NODES& nodes, HMR_NODE_NAMES* names)
{
    //Load the file path.
    FILE* contig_file;
    if (!bin_open(filepath, &contig_file, "rb"))
    {
        time_error(-1, "Failed to open contig file %s", filepath);
    }
    //Read the information.
    int32_t contig_size;
    fread(&contig_size, sizeof(int32_t), 1, contig_file);
    nodes.resize(contig_size);
    fread(nodes.data(), sizeof(HMR_NODE), contig_size, contig_file);
    //Check whether we need to read the names or not.
    if (names)
    {
        HMR_NODE_NAMES& node_names = *names;
        node_names.resize(contig_size);
        for (int32_t i = 0; i < contig_size; ++i)
        {
            fread(&node_names[i].name_size, sizeof(int32_t), 1, contig_file);
            node_names[i].name = static_cast<char*>(malloc(node_names[i].name_size+1));
            assert(node_names[i].name);
            fread(node_names[i].name, sizeof(char), node_names[i].name_size, contig_file);
            node_names[i].name[node_names[i].name_size] = '\0';
        }
    }
    fclose(contig_file);
}

bool hmr_graph_save_contigs(const char* filepath, const HMR_CONTIGS& nodes)
{
    //Load the file path.
    FILE* contig_file;
    if (!bin_open(filepath, &contig_file, "wb"))
    {
        time_error(-1, "Failed to save contig file %s", filepath);
    }
    //Write the file of the size.
    int32_t contig_size = static_cast<int32_t>(nodes.contigs.size());
    fwrite(&contig_size, sizeof(int32_t), 1, contig_file);
    //Now write contig information.
    fwrite(nodes.contigs.data(), sizeof(HMR_NODE), nodes.contigs.size(), contig_file);
    for (int32_t i = 0; i < contig_size; ++i)
    {
        const HMR_NODE_NAME& contig_name = nodes.names[i];
        //Write the contig size.
        fwrite(&contig_name.name_size, sizeof(int32_t), 1, contig_file);
        fwrite(contig_name.name, sizeof(char), contig_name.name_size, contig_file);
    }
    //Close the file.
    fclose(contig_file);
    return true;
}

void hmr_graph_load_contig_ids(const char* filepath, HMR_CONTIG_ID_VEC& contig_ids)
{
    FILE* contig_ids_file;
    if (!bin_open(filepath, &contig_ids_file, "rb"))
    {
        time_error(-1, "Failed to load contig ids from %s", filepath);
    }
    //Read the number of contig ids.
    int32_t id_sizes;
    fread(&id_sizes, sizeof(int32_t), 1, contig_ids_file);
    contig_ids.resize(id_sizes);
    fread(contig_ids.data(), sizeof(int32_t), id_sizes, contig_ids_file);
    fclose(contig_ids_file);
}

bool hmr_graph_save_contig_ids(const char* filepath, const HMR_CONTIG_ID_VEC& contig_ids)
{
    FILE* contig_ids_file;
    if (!bin_open(filepath, &contig_ids_file, "wb"))
    {
        time_error(-1, "Failed to save contig ids to file %s", filepath);
        return false;
    }
    //Write the number of contig ids.
    int32_t id_sizes = static_cast<int32_t>(contig_ids.size());
    fwrite(&id_sizes, sizeof(int32_t), 1, contig_ids_file);
    fwrite(contig_ids.data(), sizeof(int32_t), id_sizes, contig_ids_file);
    fclose(contig_ids_file);
    return true;
}

void hmr_graph_load_edges(const char* filepath, int32_t buf_size, GRAPH_EDGE_SIZE_PROC size_proc, GRAPH_EDGE_PROC proc, void* user)
{
    hmr_graph_load_with_buffer(filepath, buf_size, size_proc, proc, user);
}

bool hmr_graph_save_edges(const char* filepath, const HMR_EDGE_COUNTERS& edges)
{
    FILE* edge_file;
    if (!bin_open(filepath, &edge_file, "wb"))
    {
        time_error(-1, "Failed to save edges to file %s", filepath);
    }
    //Write the number of edges.
    uint64_t edge_sizes = static_cast<uint64_t>(edges.size());
    fwrite(&edge_sizes, sizeof(uint64_t), 1, edge_file);
    //Write the edge data.
    fwrite(edges.data(), sizeof(HMR_EDGE_INFO), edge_sizes, edge_file);
    fclose(edge_file);
    return false;
}

bool hmr_graph_allele_conflict(const HMR_ALLELE_TABLE& allele_table, int32_t contig_id_a, int32_t contig_id_b)
{
    //Try to find allele record of contig a.
    const auto a_iter = allele_table.find(contig_id_a);
    //If there is no record for contig a, it means it is okay with everyone.
    if (a_iter == allele_table.end())
    {
        return false;
    }
    //If the contig b appears in a's record, it conflicts.
    return hmr_in_ordered_vector(contig_id_b, a_iter->second);
}

void hmr_graph_load_allele_table(const char* filepath, HMR_ALLELE_TABLE& allele_table)
{
    FILE* allele_table_file;
    if (!bin_open(filepath, &allele_table_file, "rb"))
    {
        time_error(-1, "Failed to open allele table file %s", filepath);
    }
    //Read the record counts.
    int32_t record_sizes;
    fread(&record_sizes, sizeof(int32_t), 1, allele_table_file);
    allele_table.reserve(record_sizes);
    for (int32_t i = 0; i < record_sizes; ++i)
    {
        //Read the contig id.
        int32_t contig_id, record_size;
        fread(&contig_id, sizeof(int32_t), 1, allele_table_file);
        fread(&record_size, sizeof(int32_t), 1, allele_table_file);
        //Read the conflict ids.
        HMR_CONTIG_ID_VEC conflict_ids;
        conflict_ids.resize(record_size);
        fread(conflict_ids.data(), sizeof(int32_t), record_size, allele_table_file);
        allele_table.insert(std::make_pair(contig_id, conflict_ids));
    }
}

bool hmr_graph_save_allele_table(const char* filepath, const HMR_ALLELE_TABLE& allele_table)
{
    FILE* allele_table_file;
    if (!bin_open(filepath, &allele_table_file, "wb"))
    {
        time_error(-1, "Failed to save allele table to file %s", filepath);
        return false;
    }
    //Write the record counts.
    int32_t record_sizes = static_cast<int32_t>(allele_table.size());
    fwrite(&record_sizes, sizeof(int32_t), 1, allele_table_file);
    for (const auto &record : allele_table)
    {
        fwrite(&record.first, sizeof(int32_t), 1, allele_table_file);
        int32_t record_size = static_cast<int32_t>(record.second.size());
        fwrite(&record_size, sizeof(int32_t), 1, allele_table_file);
        fwrite(record.second.data(), sizeof(int32_t), record_size, allele_table_file);
    }
    return true;
}

void hmr_graph_load_reads(const char* filepath, int32_t buf_size, HMR_READS_PROC proc, void* user)
{
    hmr_graph_load_with_buffer(filepath, buf_size, NULL, proc, user);
}
