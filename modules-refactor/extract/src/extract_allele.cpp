#include "hmr_global.hpp"
#include "hmr_text_file.hpp"
#include "hmr_ui.hpp"

#include "extract_allele.hpp"

std::list<size_t> allele_find_tab(char* line, size_t len)
{
    std::list<size_t> stops;
    for (size_t i = 0; i < len; ++i)
    {
        if (line[i] == '\t')
        {
            stops.push_back(i);
        }
    }
    return stops;
}

void allele_insert_record(HMR_ALLELE_TABLE& allele_table, int32_t id_a, int32_t id_b)
{
    auto iter_a = allele_table.find(id_a);
    if (iter_a == allele_table.end())
    {
        HMR_CONTIG_ID_SET a_set;
        a_set.insert(id_b);
        allele_table.insert(std::make_pair(id_a, a_set));
    }
    else
    {
        iter_a->second.insert(id_b);
    }
}

HMR_ALLELE_TABLE extract_allele_table(const char* filepath, CONTIG_INDEX_MAP* index_map)
{
    //Open the text file.
    TEXT_LINE_HANDLE line_handle;
    if (!text_open_read_line(filepath, &line_handle))
    {
        time_error(-1, "Failed to read allele table file %s", filepath);
    }
    //Read the allele table line by line.
    ssize_t line_size = 0;
    char* line = NULL;
    size_t len = 0;
    HMR_ALLELE_TABLE allele_table;
    while ((line_size = (line_handle.parser(&line, &len, &line_handle.buf, line_handle.file_handle))) != -1)
    {
        //Finding all the '\t'.
        std::list<size_t> column_stops = allele_find_tab(line, static_cast<size_t>(line_size));
        //We only care about the lines with 4 or more columns.
        if (column_stops.size() < 3)
        {
            continue;
        }
        //Try to read the column indices.
        HMR_CONTIG_ID_CHAIN contig_id_chain;
        //Ignore the first column (Chromosome label).
        column_stops.pop_front();
        //Record the second column start.
        size_t column_start = column_stops.front() + 1;
        column_stops.pop_front();
        //Extract the rest of the column data.
        while (!column_stops.empty())
        {
            //Get the contig name.
            std::string contig_name(line + column_start, column_stops.front() - column_start);
            int32_t contig_id = extract_contig_index_get(*index_map, contig_name);
            if (contig_id != -1)
            {
                contig_id_chain.push_back(contig_id);
            }
            column_start = column_stops.front() + 1;
            column_stops.pop_front();
        }
        //Check whether the column start reachs the end.
        if (column_start != line_size - 1)
        {
            std::string contig_name(line + column_start, line_size - column_start);
            int32_t contig_id = extract_contig_index_get(*index_map, contig_name);
            if (contig_id != -1)
            {
                contig_id_chain.push_back(contig_id);
            }
        }
        HMR_CONTIG_ID_VEC contig_ids;
        hMoveListToVector(contig_id_chain, contig_ids);
        for (size_t i = 0; i < contig_ids.size() - 1; ++i)
        {
            int32_t id_i = contig_ids[i];
            for (size_t j = i + 1; j < contig_ids.size(); ++j)
            {
                int32_t id_j = contig_ids[j];
                allele_insert_record(allele_table, id_i, id_j);
                allele_insert_record(allele_table, id_j, id_i);
            }
        }
    }
    //Close the file.
    text_close_read_line(&line_handle);
    return allele_table;
}
