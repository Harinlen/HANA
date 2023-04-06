#include <cstdlib>

#include "args_orientation.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-n", "--nodes"}, "NDOES", "HMR contig node file (.hmr_contig)", LAMBDA_PARSE_ARG {opts.nodes = arg[0]; }},
    { {"-r", "--reads"}, "READS", "HMR paired-reads file (.hmr_reads)", LAMBDA_PARSE_ARG { opts.reads = arg[0];}},
    { {"-s", "--seq"}, "SEQ", "HMR sorted contig sequence file (.hmr_seq)", LAMBDA_PARSE_ARG { opts.seq = arg[0];}},
    { {"-o", "--output"}, "OUTPUT", "Output chromosome sequence file (.hmr_chromo)", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
    { {"-b", "--buffer-size"}, "BUFFER_SIZE", "Buffered reads count (default: 65535)", LAMBDA_PARSE_ARG { opts.buffer_size = atoi(arg[0]);}},
};
