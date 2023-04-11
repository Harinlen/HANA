#include <cstdlib>

#include "args_orientation.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-n", "--nodes"}, "NDOES", "HMR contig node file (.hmr_contig)", LAMBDA_PARSE_ARG {opts.nodes = arg[0]; }},
    { {"-r", "--reads"}, "READS", "HMR paired-reads file (.hmr_reads)", LAMBDA_PARSE_ARG { opts.reads = arg[0];}},
    { {"-s", "--seq"}, "SEQ 1, SEQ 2...", "HMR sorted contig sequence file (.hmr_seq)", LAMBDA_PARSE_ARG { opts.seq = arg;}},
    { {"-b", "--buffer-size"}, "BUFFER_SIZE", "HMR paired-reads buffer size (unit: K, default: 512)", LAMBDA_PARSE_ARG {opts.read_buffer_size = atoi(arg[0]); }},
};
