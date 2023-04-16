#include <unordered_map>

#include "hmr_path.hpp"
#include "hmr_text_file.hpp"
#include "hmr_contig_graph.hpp"
#include "hmr_args.hpp"
#include "hmr_ui.hpp"
#include "hmr_bam.hpp"

#include "args_dump.hpp"

typedef int (*DUMP_PROC)(int, char* []);

extern HMR_ARGS opts;

typedef struct DUMP_BAM
{
    FILE* dump_file;
    std::string* names;
    int32_t i;
} DUMP_BAM;

void dump_bam_n_contig(uint32_t num_of_contigs, void* user)
{
    DUMP_BAM* d = static_cast<DUMP_BAM*>(user);
    d->names = new std::string[num_of_contigs];
    d->i = 0;
}
void dump_bam_contig(uint32_t name_length, char* name, uint32_t, void* user)
{
    DUMP_BAM* d = static_cast<DUMP_BAM*>(user);
    d->names[d->i] = std::string(name, name_length);
    ++d->i;
}
void dump_bam_read_align(size_t, const BAM_BLOCK_HEADER* bam_block, void* user)
{
    DUMP_BAM* d = static_cast<DUMP_BAM*>(user);
    fprintf(d->dump_file, "%s\t%d\t%s\t%d\n", d->names[bam_block->refID].c_str(), bam_block->pos, d->names[bam_block->next_refID].c_str(), bam_block->next_pos);
}

int dump_bam(int argc, char* argv[])
{
    if (argc < 2)
    {
        printf("Usage: bam [bam file path] -o [dump file path]\n");
        exit(-1);
    }
    //Read the arguments.
    parse_arguments(argc - 1, argv + 1);
    //Get the bam file path.
    const char* filepath = argv[1];
    if (!path_can_read(filepath))
    {
        printf("Failed to open bam file %s\n", filepath);
        exit(-1);
    }
    std::string output_path = std::string(filepath) + ".txt";
    FILE* dump_file;
    if (!text_open_write(output_path.c_str(), &dump_file))
    {
        time_error(-1, "Failed to open output file.");
        exit(-1);
    }
    //Start parsing the bam file.
    DUMP_BAM user;
    user.dump_file = dump_file;
    hmr_bam_read(filepath, BAM_MAPPING_PROC{ dump_bam_n_contig, dump_bam_contig, dump_bam_read_align }, &user, 1);
    fclose(dump_file);
    return 0;
}

typedef struct DUMP_READS
{
    FILE *file;
    HMR_CONTIGS node_infos;
} DUMP_READS;

void dump_reads_proc(HMR_MAPPING* mapping, int32_t buf_size, void *user)
{
    DUMP_READS *d = static_cast<DUMP_READS*>(user);
    for(int32_t i=0; i<buf_size; ++i)
    {
        fprintf(d->file, "%s\t%d\t%s\t%d\n",
                d->node_infos.names[mapping[i].refID].name, mapping[i].pos,
                d->node_infos.names[mapping[i].next_refID].name, mapping[i].next_pos);
    }
}

int dump_reads(int argc, char* argv[])
{
    if (argc < 2)
    {
        printf("Usage: group [group file path] -n [nodes file path]\n");
    }
    //Read the arguments.
    parse_arguments(argc - 1, argv + 1);
    //Get the bam file path.
    const char* filepath = argv[1];
    if (!path_can_read(filepath))
    {
        printf("Failed to open group file %s\n", filepath);
        exit(-1);
    }
    if (!opts.nodes)
    {
        printf("Please provide the node file path.\n");
        exit(-1);
    }
    if (!path_can_read(opts.nodes))
    {
        printf("Failed to open node file %s\n", opts.nodes);
        exit(-1);
    }
    std::string output_path = std::string(filepath) + ".txt";
    FILE* dump_file;
    if (!text_open_write(output_path.c_str(), &dump_file))
    {
        time_error(-1, "Failed to open output file.");
    }
    //Load the contigs.
    DUMP_READS user;
    user.file = dump_file;
    hmr_graph_load_contigs(opts.nodes, user.node_infos.contigs, &user.node_infos.names);
    //Load the reads.
    hmr_graph_load_reads(filepath, 512 << 10, dump_reads_proc, &user);
    fclose(dump_file);
    return 0;
}

void dump_edges_size_proc(uint64_t, void*)
{
}

void dump_edges_proc(HMR_EDGE_INFO* edges, int32_t edge_size, void* user)
{
    DUMP_READS *d = static_cast<DUMP_READS*>(user);
    for(int32_t i=0; i<edge_size; ++i)
    {
        fprintf(d->file, "%s\t%s\t%lu\t%lf\n",
                d->node_infos.names[edges[i].start].name, d->node_infos.names[edges[i].end].name,
                edges[i].pairs, edges[i].weights);
    }
}

int dump_edges(int argc, char* argv[])
{
    if (argc < 2)
    {
        printf("Usage: edges [edges file path] -n [nodes file path]\n");
    }
    //Read the arguments.
    parse_arguments(argc - 1, argv + 1);
    //Get the bam file path.
    const char* filepath = argv[1];
    if (!path_can_read(filepath))
    {
        printf("Failed to open group file %s\n", filepath);
        exit(-1);
    }
    if (!opts.nodes)
    {
        printf("Please provide the node file path.\n");
        exit(-1);
    }
    if (!path_can_read(opts.nodes))
    {
        printf("Failed to open node file %s\n", opts.nodes);
        exit(-1);
    }
    std::string output_path = std::string(filepath) + ".txt";
    FILE* dump_file;
    if (!text_open_write(output_path.c_str(), &dump_file))
    {
        time_error(-1, "Failed to open output file.");
    }
    //Load the contigs.
    DUMP_READS user;
    user.file = dump_file;
    hmr_graph_load_contigs(opts.nodes, user.node_infos.contigs, &user.node_infos.names);
    //Load the edges.
    hmr_graph_load_edges(filepath, 512 << 10, dump_edges_size_proc, dump_edges_proc, &user);
    fclose(dump_file);
    return 0;
}

int dump_nodes(int argc, char* argv[])
{
    if (argc < 1)
    {
        printf("Usage: nodes [nodes file path]\n");
        exit(-1);
    }
    //Get the bam file path.
    const char* filepath = argv[1];
    if (!path_can_read(filepath))
    {
        printf("Failed to open nodes file %s\n", filepath);
        exit(-1);
    }
    std::string output_path = std::string(filepath) + ".txt";
    FILE* dump_file;
    if (!text_open_write(output_path.c_str(), &dump_file))
    {
        time_error(-1, "Failed to open output file.");
        exit(-1);
    }
    //Start parsing the bam file.
    HMR_CONTIGS node_infos;
    hmr_graph_load_contigs(filepath, node_infos.contigs, &node_infos.names);
    //Write the contig name and length.
    fprintf(dump_file, "#\tContig\tRECounts\tLength\n");
    for (size_t i = 0; i < node_infos.contigs.size(); ++i)
    {
        fprintf(dump_file, "%zu\t%s\t%d\t%d\n", i, node_infos.names[i].name, node_infos.contigs[i].enzyme_count, node_infos.contigs[i].length);
    }
    fclose(dump_file);
    return 0;
}

int dump_group(int argc, char* argv[])
{
    if (argc < 2)
    {
        printf("Usage: group [group file path] -n [nodes file path]\n");
    }
    //Read the arguments.
    parse_arguments(argc - 1, argv + 1);
    //Get the bam file path.
    const char* filepath = argv[1];
    if (!path_can_read(filepath))
    {
        printf("Failed to open group file %s\n", filepath);
        exit(-1);
    }
    if (!opts.nodes)
    {
        printf("Please provide the node file path.\n");
        exit(-1);
    }
    if (!path_can_read(opts.nodes))
    {
        printf("Failed to open node file %s\n", opts.nodes);
        exit(-1);
    }
    std::string output_path = std::string(filepath) + ".txt";
    FILE* dump_file;
    if (!text_open_write(output_path.c_str(), &dump_file))
    {
        time_error(-1, "Failed to open output file.");
    }
    //Load the contigs.
    HMR_CONTIGS node_infos;
    hmr_graph_load_contigs(opts.nodes, node_infos.contigs, &node_infos.names);
    //Load the group.
    HMR_CONTIG_ID_VEC contig_ids;
    hmr_graph_load_contig_ids(filepath, contig_ids);
    fprintf(dump_file, "#\tContig\tLength\n");
    for (int32_t contig_id : contig_ids)
    {
        fprintf(dump_file, "%d\t%s\t%d\n", contig_id, node_infos.names[contig_id].name, node_infos.contigs[contig_id].length);
    }
    fclose(dump_file);
    return 0;
}

int dump_allele_table(int argc, char* argv[])
{
    if (argc < 2)
    {
        printf("Usage: allele [allele table file path] -n [nodes file path]\n");
    }
    //Read the arguments.
    parse_arguments(argc - 1, argv + 1);
    //Get the bam file path.
    const char* filepath = argv[1];
    if (!path_can_read(filepath))
    {
        printf("Failed to open group file %s\n", filepath);
        exit(-1);
    }
    if (!opts.nodes)
    {
        printf("Please provide the node file path.\n");
        exit(-1);
    }
    if (!path_can_read(opts.nodes))
    {
        printf("Failed to open node file %s\n", opts.nodes);
        exit(-1);
    }
    std::string output_path = std::string(filepath) + ".txt";
    FILE* dump_file;
    if (!text_open_write(output_path.c_str(), &dump_file))
    {
        time_error(-1, "Failed to open output file.");
    }
    //Load the contigs.
    HMR_CONTIGS node_infos;
    hmr_graph_load_contigs(opts.nodes, node_infos.contigs, &node_infos.names);
    //Load the allele table.
    HMR_CONTIG_ID_TABLE allele_table;
    hmr_graph_load_contig_table(filepath, allele_table);
    for (const HMR_CONTIG_ID_VEC &record : allele_table)
    {
        for (const int32_t cid: record)
        {
            fprintf(dump_file, "%s\t", node_infos.names[cid].name);
        }
        fprintf(dump_file, "\n");
    }
    fclose(dump_file);
    return 0;
}



std::unordered_map<std::string, DUMP_PROC> dump_proc_map = {
    {"bam", dump_bam},
    {"nodes", dump_nodes},
    {"group", dump_group},
    {"reads", dump_reads},
    {"edges", dump_edges},
    {"allele", dump_allele_table},
};

void help_exit()
{
    printf("usage: op [param]\n");
    printf("Support operations:\n");
    for (const auto &iter : dump_proc_map)
    {
        printf("    %s\n", iter.first.c_str());
    }
    exit(-1);
}

int main(int argc, char* argv[])
{
    //Check the operations.
    if (argc < 2)
    {
        help_exit();
    }
    std::string op = argv[1];
    const auto iter = dump_proc_map.find(op);
    if (iter == dump_proc_map.end())
    {
        printf("Unknown operation '%s'", op.c_str());
        help_exit();
    }
    //Or else, call the target function.
    return iter->second(argc - 1, argv + 1);
}
