#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <unordered_map>

#include "hmr_args.hpp"
#include "hmr_contig_graph.hpp"
#include "hmr_fasta.hpp"
#include "hmr_path.hpp"
#include "hmr_text_file.hpp"
#include "hmr_ui.hpp"

#include "args_build.hpp"

extern HMR_ARGS opts;

typedef struct CONTIG_SEQ
{
    char *seq;
    size_t seq_size;
} CONTIG_SEQ;

typedef std::unordered_map<int32_t, CONTIG_SEQ> CONTIG_DICT;
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

void build_fasta_loader(int32_t seq_index, char *, size_t, char *seq_data, size_t seq_data_len, void *user)
{
    CONTIG_DICT *contigs = reinterpret_cast<CONTIG_DICT *>(user);
    //Construct the contig building dictionary.
    contigs->insert(std::make_pair(seq_index, CONTIG_SEQ{ seq_data, seq_data_len }));
}

typedef struct CHROMOSOME_BUILD
{
    CONTIG_DICT contigs;
    FILE *output_fasta;
} CHROMOSOME_BUILD;

void build_chromosome(const HMR_DIRECTED_CONTIG &contig_info, CHROMOSOME_BUILD& builder)
{
    CONTIG_DICT &contigs = builder.contigs;
    //Extract the sequence from the directory.
    auto seq_iter = contigs.find(contig_info.id);
    if(seq_iter == contigs.end())
    {
        return;
    }
    CONTIG_SEQ contig = seq_iter->second;
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
    time_print("Constructing FASTA sequence index from %s", opts.fasta);
    hmr_fasta_read(opts.fasta, build_fasta_loader, &builder.contigs);
    time_print("%zu sequences loaded.", builder.contigs.size());
    //Open the output file to write data.
    time_print("Opening output file %s", opts.output);
    if(!text_open_write(opts.output, &builder.output_fasta))
    {
        time_error(-1, "Failed to open the output file.");
    }
    for(size_t i=0; i<opts.chromosomes.size(); ++i)
    {
        char *chromo_path = opts.chromosomes[i];
        time_print("Constructing chromosome %zu from %s", i+1, chromo_path);
        //Write the name of the chromosome.
        fprintf(builder.output_fasta, ">Chromosome_%zu\n", i+1);
        CHROMOSOME_CONTIGS seq;
        hmr_graph_load_chromosome(chromo_path, seq);
        for (const auto &contig_info: seq)
        {
            build_chromosome(contig_info, builder);
        }
        fprintf(builder.output_fasta, "\n");
    }
    fclose(builder.output_fasta);
    time_print("Build complete.");
    return 0;
}
