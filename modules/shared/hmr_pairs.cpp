#include <vector>

#include "hmr_text_file.hpp"
#include "hmr_ui.hpp"

#include "hmr_pairs.hpp"

void hmr_pairs_fill_tabs(char *line, size_t line_size, size_t *tab_stops, int &tab_stop_count)
{
    //Find 7 tabs from the line.
    tab_stop_count = 0;
    for(size_t i=0; i<line_size && tab_stop_count < 7; ++i)
    {
        if(line[i] == '\t')
        {
            tab_stops[tab_stop_count] = i;
            line[i] = '\0';
            ++tab_stop_count;
        }
    }
}

void hmr_pairs_read(const char *filepath, PAIR_PROC proc, void *user)
{
    //Read the file path as a text file.
    TEXT_LINE_HANDLE line_handle;
    if(!text_open_read_line(filepath, &line_handle))
    {
        time_error(-1, "Failed to read pairs file %s", filepath);
    }
    //Loop and detect line.
    char *line = NULL;
    size_t len = 0, line_length = 0;
    ssize_t line_size = 0;
    //Line parsing variables.
    size_t tab_stops[7];
    int tab_stop_count;
    int32_t line_pos_0, line_pos_1;
    while((line_size = (line_handle.parser(&line, &len, &line_handle.buf, line_handle.file_handle))) != -1)
    {
        line_length = static_cast<size_t>(line_size);
        //We are caring about column 1, 2, 3, 4 and 7. (starts with 0)
        printf("%c%c%c%c\t%u\t%zu\n", line[0], line[1], line[2], line[3], line[0], line_length);
        //Ignore the # start comment lines.
        if(line_length < 1 || line[0] == '#')
        {
            continue;
        }
        //Find out all the tab line from the line.
        hmr_pairs_fill_tabs(line, line_length, tab_stops, tab_stop_count);
        if(tab_stop_count < 7)
        {
            continue;
        }
        //Increase the position to be the start position.
        ++tab_stops[0]; ++tab_stops[2];
        line_pos_0 = atoi(line + tab_stops[1] + 1);
        line_pos_1 = atoi(line + tab_stops[3] + 1);
        proc(line + tab_stops[0], tab_stops[1] - tab_stops[0],
                line + tab_stops[2], tab_stops[3] - tab_stops[2],
                line_pos_0, line_pos_1, line + tab_stops[6] + 1, user);
    }
}
