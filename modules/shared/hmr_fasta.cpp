#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstdint>
#include <cstring>

#include "hmr_text_file.hpp"
#include "hmr_bin_file.hpp"
#include "hmr_seq.hpp"
#include "hmr_ui.hpp"
#include "hmr_char.hpp"

#include "hmr_fasta.hpp"

typedef struct FASTA_PARSE
{
    char *seq_name, *seq_data;
    size_t seq_name_len, seq_data_len;
    FASTA_PROC user_parser;
    void *user;
} FASTA_PARSE;

void fasta_name_trim(char* name, size_t *name_len)
{
    for (size_t i = 0, i_max = *name_len; i < i_max; ++i)
    {
        if (is_space(name[i]))
        {
            name[i] = '\0';
            *name_len = i;
            return;
        }
    }
}

inline void fasta_yield_line(char *line, size_t line_size, FASTA_PARSE &args, int32_t &index)
{
    if(line_size == 0 || line == NULL)
    {
        return;
    }
    //Check whether the line start is with '>'.
    if(line[0] == '>')
    {
        //Check the last is empty or not.
        if(args.seq_name != NULL)
        {
            if(args.seq_data != NULL)
            {
                //Upper the sequence.
                hmr_seq_upper(args.seq_data, args.seq_data_len);
                //Yield the line parser.
                args.user_parser(index, args.seq_name, args.seq_name_len, args.seq_data, args.seq_data_len, args.user);
                ++index;
            }
            else
            {
                //Clear the sequence name.
                free(args.seq_name);
            }
        }
        //Set the sequence name.
        args.seq_name = static_cast<char *>(malloc(line_size));
        assert(NULL != args.seq_name);
        args.seq_name_len = line_size-1;
        memcpy(args.seq_name, line+1, args.seq_name_len);
        args.seq_name[args.seq_name_len] = '\0';
        //Trim the name to the first spacing.
        fasta_name_trim(args.seq_name, &args.seq_name_len);
        //Reset the sequence.
        args.seq_data = NULL;
        args.seq_data_len = 0;
    }
    else
    {
        //Check the sequence.
        size_t seq_extend_len = args.seq_data_len + line_size;
        char *seq_data = static_cast<char*>(realloc(args.seq_data, seq_extend_len + 1));
        assert(seq_data);
        args.seq_data = seq_data;
        //Copy the data.
        memcpy(args.seq_data + args.seq_data_len, line, line_size);
        args.seq_data_len = seq_extend_len;
        args.seq_data[args.seq_data_len] = '\0';
    }
}

void hmr_fasta_read(const char *filepath, FASTA_PROC parser, void *user)
{
    FASTA_PARSE args;
    TEXT_LINE_HANDLE line_handle;
    if(!text_open_read_line(filepath, &line_handle))
    {
        time_error(-1, "Failed to read FASTA file %s", filepath);
    }
    //Loop and detect line.
    char *line = NULL;
    size_t len = 0, line_length = 0;
    int32_t index = 0;
    ssize_t line_size = 0;
    args.seq_name = NULL;
    args.seq_data = NULL;
    args.seq_name_len = 0;
    args.seq_data_len = 0;
    args.user_parser = parser;
    args.user = user;
    while((line_size = (line_handle.parser(&line, &len, &line_handle.buf, line_handle.file_handle))) != -1)
    {
        //Trimmed line.
        line_length = static_cast<size_t>(line_size);
        trimmed_right(line, line_length);
        //Yield the line, call the function.
        fasta_yield_line(line, line_length, args, index);
    }
    //At the end of the line, yield the last result.
    hmr_seq_upper(args.seq_data, args.seq_data_len);
    parser(index, args.seq_name, args.seq_name_len, args.seq_data, args.seq_data_len, user);
    //Close the file.
    text_close_read_line(&line_handle);
}
