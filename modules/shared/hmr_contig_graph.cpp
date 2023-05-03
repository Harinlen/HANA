#include <algorithm>
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
    std::mutex loaded_mutex, processed_mutex;
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
    //Get the total file size.
    fseek(data_file, 0L, SEEK_END);
#ifdef _MSC_VER
    size_t total_size = _ftelli64(data_file);
#else
    size_t total_size = ftello64(data_file);
#endif
    fseek(data_file, 0L, SEEK_SET);
    //For UI output.
    size_t report_size = (total_size + 9) / 10, report_pos = report_size;
    //Check size loader.
    if (size_proc)
    {
        //Read the size from the file.
        uint64_t item_size;
        fread(&item_size, sizeof(uint64_t), 1, data_file);
        size_proc(item_size, user);
    }
    //Loop until all the data is loaded.
    while (!loader.finished)
    {
        //Start to load data.
        loader.buf_loading->buf_size = static_cast<int32_t>(fread(loader.buf_loading->buf, sizeof(T), loader.buf_size, data_file));
        //Check whether we can still fill all the part of the buffer.
        loader.finished = loader.buf_loading->buf_size < loader.buf_size;
        //Check should we report the position.
#ifdef _MSC_VER
        size_t data_file_pos = _ftelli64(data_file);
#else
        size_t data_file_pos = ftello64(data_file);
#endif
        if (data_file_pos >= report_pos)
        {
            float percent = static_cast<float>(data_file_pos) / static_cast<float>(total_size) * 100.0f;
            time_print("File parsed %.1f%%", percent);
            report_pos += report_size;
        }
        //Set as loaded.
        loader.loaded_mutex.lock();
        loader.is_loaded = true;
        loader.loaded_mutex.unlock();
        //Wait for loaded flag to be off.
        if (!loader.finished)
        {
            //Wait for loaded data is processing.
            while(loader.is_loaded)
            {
                std::this_thread::sleep_for(std::chrono::nanoseconds(1));
            }
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
        while(!loader.is_loaded)
        {
            std::this_thread::sleep_for(std::chrono::nanoseconds(1));
        }
        //One chunk is loaded, we flip the buffer, let the thread to loaded next part.
        auto buf_temp = loader.buf_loading;
        loader.buf_loading = loader.buf_processing;
        loader.buf_processing = buf_temp;
        //Reset the loaded flag.
        loader.loaded_mutex.lock();
        loader.is_loaded = false;
        loader.loaded_mutex.unlock();
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

std::string hmr_graph_path_chromo_name(const char* seq_path)
{
    return std::string(path_basename(seq_path) + ".hmr_chromo");
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

inline void hmr_graph_read_contig_ids(HMR_CONTIG_ID_VEC& contig_ids, FILE *file)
{
    int32_t id_sizes;
    fread(&id_sizes, sizeof(int32_t), 1, file);
    contig_ids.resize(id_sizes);
    fread(contig_ids.data(), sizeof(int32_t), id_sizes, file);
}

void hmr_graph_load_contig_ids(const char* filepath, HMR_CONTIG_ID_VEC& contig_ids)
{
    FILE* contig_ids_file;
    if (!bin_open(filepath, &contig_ids_file, "rb"))
    {
        time_error(-1, "Failed to load contig ids from %s", filepath);
    }
    //Read the number of contig ids.
    hmr_graph_read_contig_ids(contig_ids, contig_ids_file);
    fclose(contig_ids_file);
}

inline void hmr_graph_write_contig_ids(const HMR_CONTIG_ID_VEC& contig_ids, FILE *file)
{
    //Write the number of contig ids.
    int32_t id_sizes = static_cast<int32_t>(contig_ids.size());
    fwrite(&id_sizes, sizeof(int32_t), 1, file);
    fwrite(contig_ids.data(), sizeof(int32_t), id_sizes, file);
}

bool hmr_graph_save_contig_ids(const char* filepath, const HMR_CONTIG_ID_VEC& contig_ids)
{
    FILE* contig_ids_file;
    if (!bin_open(filepath, &contig_ids_file, "wb"))
    {
        time_error(-1, "Failed to save contig ids to file %s", filepath);
        return false;
    }
    hmr_graph_write_contig_ids(contig_ids, contig_ids_file);
    fclose(contig_ids_file);
    return true;
}

void hmr_graph_load_contig_table(const char *filepath, HMR_CONTIG_ID_TABLE &contig_table)
{
    FILE* contig_table_file;
    if (!bin_open(filepath, &contig_table_file, "rb"))
    {
        time_error(-1, "Failed to load contig table from %s", filepath);
    }
    //Read the number of contig ids.
    int32_t num_of_rows;
    fread(&num_of_rows, sizeof(int32_t), 1, contig_table_file);
    contig_table.reserve(num_of_rows);
    for(int32_t i=0; i<num_of_rows; ++i)
    {
        HMR_CONTIG_ID_VEC row;
        hmr_graph_read_contig_ids(row, contig_table_file);
        contig_table.push_back(row);
    }
    fclose(contig_table_file);
}

bool hmr_graph_save_contig_table(const char *filepath, const HMR_CONTIG_ID_TABLE &contig_table)
{
    FILE* contig_table_file;
    if (!bin_open(filepath, &contig_table_file, "wb"))
    {
        time_error(-1, "Failed to save contig table to file %s", filepath);
        return false;
    }
    int32_t num_of_rows = static_cast<int32_t>(contig_table.size());
    fwrite(&num_of_rows, sizeof(int32_t), 1, contig_table_file);
    for(const auto &row: contig_table)
    {
        hmr_graph_write_contig_ids(row, contig_table_file);
    }
    fclose(contig_table_file);
    return true;
}

void hmr_graph_allele_insert_record(std::unordered_map<int32_t, HMR_CONTIG_ID_SET>& allele_set_map, int32_t id_a, int32_t id_b)
{
    auto iter_a = allele_set_map.find(id_a);
    if (iter_a == allele_set_map.end())
    {
        HMR_CONTIG_ID_SET a_set;
        a_set.insert(id_b);
        allele_set_map.insert(std::make_pair(id_a, a_set));
    }
    else
    {
        iter_a->second.insert(id_b);
    }
}

void hmr_graph_allele_map_init(HMR_ALLELE_MAP &allele_map, const HMR_CONTIG_ID_TABLE &allele_table)
{
    //Go through the allele map
    std::unordered_map<int32_t, HMR_CONTIG_ID_SET> allele_set_map;
    for(const auto &contig_ids: allele_table)
    {
        for (size_t i = 0; i < contig_ids.size() - 1; ++i)
        {
            int32_t id_i = contig_ids[i];
            for (size_t j = i + 1; j < contig_ids.size(); ++j)
            {
                int32_t id_j = contig_ids[j];
                hmr_graph_allele_insert_record(allele_set_map, id_i, id_j);
                hmr_graph_allele_insert_record(allele_set_map, id_j, id_i);
            }
        }
    }
    //Convert the record into set.
    for (const auto& iter : allele_set_map)
    {
        int32_t contig_id = iter.first;
        //Convert the set into vector.
        HMR_CONTIG_ID_VEC conflict_vec(iter.second.begin(), iter.second.end());
        std::sort(conflict_vec.begin(), conflict_vec.end());
        allele_map.insert(std::make_pair(contig_id, conflict_vec));
    }
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

void hmr_graph_load_reads(const char* filepath, int32_t buf_size, HMR_READS_PROC proc, void* user)
{
    hmr_graph_load_with_buffer(filepath, buf_size, NULL, proc, user);
}

void hmr_graph_load_chromosome(const char* filepath, CHROMOSOME_CONTIGS& seq)
{
    FILE* chromosome_file;
    if (!bin_open(filepath, &chromosome_file, "rb"))
    {
        time_error(-1, "Failed to load contig ids from %s", filepath);
    }
    //Read the number of contig ids.
    int32_t contig_sizes;
    fread(&contig_sizes, sizeof(int32_t), 1, chromosome_file);
    seq.resize(contig_sizes);
    fread(seq.data(), sizeof(HMR_DIRECTED_CONTIG), contig_sizes, chromosome_file);
    fclose(chromosome_file);
}

bool hmr_graph_save_chromosome(const char* filepath, const CHROMOSOME_CONTIGS& seq)
{
    FILE* chromosome_file;
    if (!bin_open(filepath, &chromosome_file, "wb"))
    {
        time_error(-1, "Failed to save contig ids to file %s", filepath);
        return false;
    }
    //Write the number of contig ids.
    int32_t id_sizes = static_cast<int32_t>(seq.size());
    fwrite(&id_sizes, sizeof(int32_t), 1, chromosome_file);
    fwrite(seq.data(), sizeof(HMR_DIRECTED_CONTIG), id_sizes, chromosome_file);
    fclose(chromosome_file);
    return true;

}
