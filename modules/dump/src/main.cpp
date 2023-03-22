#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <map>
#include <string>

#include "hmr_args.hpp"
#include "hmr_ui.hpp"
#include "hmr_mapping.hpp"
#include "hmr_bin_file.hpp"
#include "hmr_contig_graph_type.hpp"
#include "hmr_contig_graph.hpp"
#include "hmr_text_file.hpp"

#include "args_dump.hpp"

extern HMR_ARGS opts;

typedef int (*DumpOpeartion)(int, char*[]);

typedef struct
{
    std::string description;
    DumpOpeartion op;
} DUMP_OP;

constexpr int DUMP_RECORD_MAX = 1048575;
constexpr size_t DUMP_BUFFER_SIZE = DUMP_RECORD_MAX * sizeof(HMR_MAPPING);

typedef struct
{
    FILE* dump_file;
    int buffer_offset;
    HMR_MAPPING* buffer;
} BAM_DUMP_USER;

void dump_bam_mapping_draft_n_contig(uint32_t n_ref, void* user) {}
void dump_bam_mapping_draft_contig(uint32_t name_length, char* name, uint32_t length, void* user) {}
void dump_bam_mapping_draft_read_align(size_t id, const MAPPING_INFO& mapping_info, void* user)
{
    //Convert the user to const user.
    BAM_DUMP_USER* dump_user = reinterpret_cast<BAM_DUMP_USER*>(user);
    //Write the index, pos and quality to the dump file.
    if (dump_user->buffer_offset == DUMP_RECORD_MAX)
    {
        fwrite(dump_user->buffer, DUMP_BUFFER_SIZE, 1, dump_user->dump_file);
        dump_user->buffer_offset = 0;
    }
    //Write the text line buffer.
    dump_user->buffer[dump_user->buffer_offset++] = HMR_MAPPING{ mapping_info.refID, mapping_info.pos, mapping_info.next_refID, mapping_info.next_pos };
}

int dump_bam(int argc, char* argv[])
{
    //Prase the arguments.
    parse_arguments(argc, argv);
    //Loop in all the bam file and dump to output directory.
    time_print("Dumping BAM file %s", opts.bam_file);
    BAM_DUMP_USER dump_user {NULL, 0, NULL};
    dump_user.buffer = static_cast<HMR_MAPPING*>(malloc(DUMP_BUFFER_SIZE));
    assert(dump_user.buffer);
    if (!bin_open(opts.output, &dump_user.dump_file, "w"))
    {
        time_error(-1, "Failed to open the output file %s", opts.output);
    }
    hmr_mapping_read(opts.bam_file, MAPPING_PROC{ dump_bam_mapping_draft_n_contig, dump_bam_mapping_draft_contig, dump_bam_mapping_draft_read_align }, &dump_user, opts.threads);
    if (dump_user.buffer_offset)
    {
        fwrite(dump_user.buffer, dump_user.buffer_offset * sizeof(HMR_MAPPING), 1, dump_user.dump_file);
    }
    fclose(dump_user.dump_file);
    return 0;
}

typedef struct DUMP_EDGE_USER
{
    FILE* text_out;
} DUMP_EDGE_USER;

void dump_edge_size(uint64_t size, void* user) {}
void dump_edge_proc(const HMR_EDGE_INFO& edge_info, void* user)
{
    DUMP_EDGE_USER* edge_user = reinterpret_cast<DUMP_EDGE_USER*>(user);
    fprintf(edge_user->text_out, "%-14d  %-14d %-14d %lf\n", edge_info.edge.pos.start, edge_info.edge.pos.end, edge_info.pairs, edge_info.weights);
}

int dump_edge(int argc, char* argv[])
{
    //Prase the arguments.
    parse_arguments(argc, argv);
    //Loop and print the data into text file.
    time_print("Dumping edge info from %s", opts.edge_file);
    DUMP_EDGE_USER user {NULL};
    if (!text_open_write(opts.output, &user.text_out))
    {
        time_error(-1, "Failed to open output file to dump the data.");
    }
    //Read the edge info.
    hmr_graph_load_edge(opts.edge_file, dump_edge_size, dump_edge_proc, &user);
    return 0;
}

std::map<std::string, DUMP_OP> dump_ops = {
    {"bam", DUMP_OP {"Dump the information of the bam file", &dump_bam}},
    {"edge", DUMP_OP {"Dump the information of the HMR edge file", &dump_edge}},
};

void exit_command_help(int exitCode)
{
    printf("usage: <command> [<options>]\n");
    printf("command arguments:\n");
    for (auto const& op_iter : dump_ops)
    {
        printf("    %-10s\t%s\n", op_iter.first.c_str(), op_iter.second.description.c_str());
    }
    exit(exitCode);
}

void lower_command(char* source)
{
    int source_length = strlen(source);
    for (int i = 0; i < source_length; ++i)
    {
        if (source[i] >= 'A' && source[i] <= 'Z')
        {
            source[i] = 'a' + source[i] - 'A';
        }
    }
}

int main(int argc, char* argv[])
{
    //Decide the operation based on the first argument.
    if (argc < 2)
    {
        exit_command_help(0);
    }
    lower_command(argv[1]);
    std::string command = std::string(argv[1]);
    auto command_iter = dump_ops.find(command);
    if (command_iter == dump_ops.end())
    {
        printf("Unknown command '%s'.", command.c_str());
        exit_command_help(1);
    }
    //Call the function of the command.
    return command_iter->second.op(argc - 1, argv + 1);
}