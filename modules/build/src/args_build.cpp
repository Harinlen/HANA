#include <cstdlib>

#include "args_build.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-f", "--fasta"}, "FASTA", "Contig FASTA file (.fasta/.fasta.gz)", LAMBDA_PARSE_ARG {opts.fasta = arg[0]; }},
    { {"-c", "--chromosome"}, "CHROMO 1, CHROMO 2...", "Hi-C reads mapping files (.bam/.hmr_mapping)", LAMBDA_PARSE_ARG { opts.chromosomes = arg; }},
    { {"-o", "--output"}, "OUTPUT", "Output build file prefix", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
};
