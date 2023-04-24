#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <list>
#include <vector>

#include "hmr_args.hpp"
#include "hmr_contig_graph.hpp"
#include "hmr_fasta.hpp"
#include "hmr_path.hpp"
#include "hmr_global.hpp"
#include "hmr_text_file.hpp"
#include "hmr_ui.hpp"

#include "args_build.hpp"

extern HMR_ARGS opts;

typedef struct CONTIG_SEQ
{
    char* name, * seq;
    size_t name_size, seq_size;
} CONTIG_SEQ;

typedef std::vector<CONTIG_SEQ> CONTIG_DICT;
typedef std::list<CONTIG_SEQ> CONTIG_DICT_LIST;


static char reverse_bp[26] = {
    'T', // 'A',
    'B',
    'G', //'C',
    'D',
    'E',
    'F',
    'C', //'G',
    'H',
    'I',
    'J',
    'K',
    'L',
    'M',
    'N',
    'O',
    'P',
    'Q',
    'R',
    'S',
    'A', //'T',
    'U',
    'V',
    'W',
    'X',
    'Y',
    'Z',
};

void build_fasta_loader(int32_t seq_index, char *name, size_t name_len, char *seq_data, size_t seq_data_len, void *user)
{
    CONTIG_DICT_LIST *contigs = reinterpret_cast<CONTIG_DICT_LIST*>(user);
    //Construct the contig building dictionary.
    contigs->push_back(CONTIG_SEQ{ name, seq_data, name_len, seq_data_len });
}

typedef struct CHROMOSOME_BUILD
{
    CONTIG_DICT contigs;
    FILE *output_fasta, *output_agp;
} CHROMOSOME_BUILD;

void build_chromosome(const HMR_DIRECTED_CONTIG &contig_info, CHROMOSOME_BUILD& builder)
{
    CONTIG_DICT& contigs = builder.contigs;
    //Extract the sequence from the directory.
    CONTIG_SEQ &contig = contigs[contig_info.id];
    //Check the direction.
    if(contig_info.direction)
    {
        //It is reversed.
        char *rev_buffer = static_cast<char *>(malloc(contig.seq_size + 1)),
                *seq = contig.seq;
        assert(rev_buffer);
        size_t i_max = contig.seq_size;
        for(size_t i=0; i<i_max; ++i)
        {
            size_t rev_idx = i_max - i - 1;
            if(seq[rev_idx] >= 'A' && seq[rev_idx] <= 'Z')
            {
                rev_buffer[i] = reverse_bp[seq[rev_idx] - 'A'];
            }
        }
        fwrite(rev_buffer, i_max, 1, builder.output_fasta);
        free(rev_buffer);
    }
    else
    {
        //Just simply write the sequence.
        fwrite(contig.seq, contig.seq_size, 1, builder.output_fasta);
    }
}

size_t build_agp(size_t offset, std::string title, const HMR_DIRECTED_CONTIG& contig_info, CHROMOSOME_BUILD& builder, size_t *counter)
{
    //Extract the sequence from the directory.
    const CONTIG_SEQ& contig = builder.contigs[contig_info.id];
    size_t next_offset = offset + contig.seq_size;
    fprintf(builder.output_agp, "%s\t%zu\t%zu\t%zu\tW\t%s\t1\t%zu\t%c\n", title.c_str(), offset, next_offset - 1, *counter, contig.name, contig.seq_size, contig_info.direction ? '-' : '+');
    ++(*counter);
    return next_offset;
}

int main(int argc, char *argv[])
{
    //Parse the arguments.
    parse_arguments(argc, argv);
    //Check the arguments are meet the requirements.
    if (!opts.fasta) { help_exit(-1, "Missing FASTA file path."); }
    if (!path_can_read(opts.fasta)) { time_error(-1, "Cannot read FASTA file %s", opts.fasta); }
    if (opts.chromosomes.empty()) { help_exit(-1, "Missing chromosome sequence file paths."); }
    if (!opts.output) { help_exit(-1, "Missing output fasta file path."); }
    //Load the FASTA and cache all the sequence.
    CHROMOSOME_BUILD builder;
    {
        CONTIG_DICT_LIST contig_list;
        time_print("Constructing FASTA sequence index from %s", opts.fasta);
        hmr_fasta_read(opts.fasta, build_fasta_loader, &contig_list);
        hMoveListToVector(contig_list, builder.contigs);
        time_print("%zu sequences loaded.", contig_list.size());
    }
    //Open the output file to write data.
    assert(opts.output);
    std::string fasta_path = std::string(opts.output) + "_build.fasta";
    std::string agp_path = std::string(opts.output) + "_build.agp";
    time_print("Opening output FASTA file %s", fasta_path.c_str());
    if(!text_open_write(fasta_path.c_str(), &builder.output_fasta))
    {
        time_error(-1, "Failed to open the output FASTA file.");
    }
    time_print("Opening output AGP file %s", agp_path.c_str());
    if (!text_open_write(agp_path.c_str(), &builder.output_agp))
    {
        time_error(-1, "Failed to open the output AGP file.");
    }
    for(size_t i=0; i<opts.chromosomes.size(); ++i)
    {
        char *chromo_path = opts.chromosomes[i];
        time_print("Constructing chromosome %zu from %s", i+1, chromo_path);
        //Write the name of the chromosome.
        std::string group_name = "Group_" + std::to_string(i + 1);
        fprintf(builder.output_fasta, ">%s\n", group_name.c_str());
        CHROMOSOME_CONTIGS seq;
        hmr_graph_load_chromosome(chromo_path, seq);
        if (seq.empty())
        {
            continue;
        }
        for (const auto &contig_info: seq)
        {
            build_chromosome(contig_info, builder);
        }
        fprintf(builder.output_fasta, "\n");
        //Build AGP file.
        size_t offset = 1, contig_counter = 1;
        //Output the first contig.
        offset = build_agp(offset, group_name, seq[0], builder, &contig_counter);
        for (size_t info_id = 1; info_id < seq.size(); ++info_id)
        {
            //Write a 100 gap.
            size_t next_offset = offset + 100;
            fprintf(builder.output_agp, "%s\t%zu\t%zu\t%zu\tU\t100\tcontig\tyes\tmap\n", group_name.c_str(), offset, next_offset - 1, contig_counter);
            ++contig_counter;
            offset = build_agp(next_offset, group_name, seq[info_id], builder, &contig_counter);
        }
    }
    fclose(builder.output_fasta);
    fclose(builder.output_agp);
    time_print("Build complete.");
    return 0;
}
