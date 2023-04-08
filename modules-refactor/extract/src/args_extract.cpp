#include <cstdlib>

#include "args_extract.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-f", "--fasta"}, "FASTA", "Contig FASTA file (.fasta/.fasta.gz)", LAMBDA_PARSE_ARG {opts.fasta = arg[0]; }},
    { {"-m", "--mapping"}, "MAPPING 1, MAPPING 2...", "Hi-C reads mapping files (.bam/.pairs)", LAMBDA_PARSE_ARG { opts.mappings = arg; }},
    { {"-o", "--output"}, "OUTPUT", "Output file prefix", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
    { {"-e", "--enzyme"}, "ENZYME", "Enzyme to find in the sequence", LAMBDA_PARSE_ARG {opts.enzyme = arg[0]; }},
    { {"-q", "--mapq"}, "MAPQ", "MAPQ of mapping lower bound (default: 40)", LAMBDA_PARSE_ARG {opts.mapq = atoi(arg[0]); }},
    { {"-t", "--threads"}, "THREADS", "Number of threads (default: 1)", LAMBDA_PARSE_ARG { opts.threads = atoi(arg[0]); }},
};
