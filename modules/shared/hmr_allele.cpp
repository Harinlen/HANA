#include "hmr_text_file.hpp"
#include "hmr_ui.hpp"

#include "hmr_allele.hpp"

void hmr_allele_table_load(const char* filepath, const HMR_CONTIGS& contig_table, ALLELE_PROC proc, void* user)
{
    //Build the contig name <-> index table.
    std::unordered_map<std::string, int32_t> contig_id_map;
    for (size_t i = 0; i < contig_table.size(); ++i)
    {
        contig_id_map.insert(std::make_pair(std::string(contig_table[i].name, contig_table[i].name_size), static_cast<int32_t>(i)));
    }
    //Read the allele table file.
    TEXT_LINE_HANDLE line_handle;
    if (!text_open_read_line(filepath, &line_handle))
    {
        time_error(-1, "Failed to read Allele table file %s", filepath);
    }
    //Loop and parse the line.
    char* line = NULL;
    size_t line_length = 0;
    ssize_t line_size = 0;
    while ((line_size = (line_handle.parser(&line, &line_length, &line_handle.buf, line_handle.file_handle))) != -1)
    {
        //Line format:
        // Chr01 \t 23399 \t tig00005935 \t tig00033851 \t tig00062802 \ttig00060970 \t
        // Do not trust the last one alway follow a \t.
        //Split the line into several parts.
        std::list<size_t> spliter_pos;
        for (size_t i = 0; i < line_length; ++i)
        {
            if (line[i] == '\t')
            {
                spliter_pos.push_back(i);
            }
        }
        //We are caring about names after the third column.
        if (spliter_pos.size() < 3)
        {
            continue;
        }
        //Pop the first item
        spliter_pos.pop_front();
        size_t name_start = spliter_pos.front() + 1;
        spliter_pos.pop_front();
        //Extract the name value.
        CONTIG_ID_VECTOR contig_ids;
        contig_ids.reserve(spliter_pos.size()+1);
        while (!spliter_pos.empty())
        {
            //Extract the contig name.
            size_t name_stop = spliter_pos.front();
            std::string contig_name(line + name_start, name_stop - name_start);
            //Find the contig index.
            auto iter = contig_id_map.find(contig_name);
            if (iter == contig_id_map.end())
            {
                time_error(-1, "Unknown contig name '%s'", contig_name.c_str());
            }
            contig_ids.push_back(iter->second);
            //Go to next column.
            name_start = name_stop + 1;
            spliter_pos.pop_front();
        }
        //Check whether we reach the end of the line.
        if (name_start < line_length - 1)
        {
            //Treat the start to the end of line is the target.
            std::string contig_name(line + name_start, line_length - name_start);
            auto iter = contig_id_map.find(contig_name);
            if (iter == contig_id_map.end())
            {
                time_error(-1, "Unknown contig name '%s'", contig_name.c_str());
            }
            contig_ids.push_back(iter->second);
        }
        //Call the processing function.
        proc(contig_ids, user);
    }
    //Close the file.
    text_close_read_line(&line_handle);
}

void hmr_allele_table_insert(HMR_ALLELE_TABLE* table, int32_t a, int32_t b)
{
    auto iter_a = table->find(a);
    if (iter_a == table->end())
    {
        //Build the set for a.
        CONTIG_ID_SET set_a;
        set_a.insert(b);
        table->insert(std::make_pair(a, set_a));
        return;
    }
    //Update the data.
    iter_a->second.insert(b);
}

void hmr_allele_table_set_build(const CONTIG_ID_VECTOR &conflict, void* user)
{
    HMR_ALLELE_TABLE* table = reinterpret_cast<HMR_ALLELE_TABLE*>(user);
    //Insert the value to table.
    size_t j_max = conflict.size(), i_max = j_max - 1;
    for (size_t i = 0; i < i_max; ++i)
    {
        for (size_t j = i + 1; j < j_max; ++j)
        {
            hmr_allele_table_insert(table, i, j);
            hmr_allele_table_insert(table, j, i);
        }
    }
}

void hmr_allele_table_conflict_set(const char* filepath, const HMR_CONTIGS& contig_table, HMR_ALLELE_TABLE& table)
{
    hmr_allele_table_load(filepath, contig_table, hmr_allele_table_set_build, &table);
}
