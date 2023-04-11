#include <cstdlib>

#include "args_ordering.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-n", "--nodes"}, "NDOES", "HMR contig node file (.hmr_contig)", LAMBDA_PARSE_ARG {opts.nodes = arg[0]; }},
    { {"-e", "--edge"}, "EDGE", "HMR edge file (.hmr_edge)", LAMBDA_PARSE_ARG { opts.edge = arg[0];}},
    { {"-g", "--group"}, "GROUP", "HMR contig group file (.hmr_group)", LAMBDA_PARSE_ARG { opts.group = arg[0];}},
    { {"-o", "--output"}, "OUTPUT", "Output ordered contig sequence file (.hmr_seq)", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
    { {"-t", "--threads"}, "THREAS", "Number of threads (default: 1)", LAMBDA_PARSE_ARG { opts.threads = atoi(arg[0]); }},
    { {"-b", "--buffer-size"}, "BUFFER_SIZE", "HMR edge buffer size (unit: K, default: 512)", LAMBDA_PARSE_ARG {opts.read_buffer_size = atoi(arg[0]); }},
    { {"-s", "--seed"}, "SEED", "Fixed random seed, 0 for no special seed (default: 0)", LAMBDA_PARSE_ARG {opts.seed = atoll(arg[0]); }},
    { {"--mutapb"}, "MUTATION", "Mutation probability (default: 0.2)", LAMBDA_PARSE_ARG {opts.mutapb = atof(arg[0]); }},
    { {"--ngen"}, "NUM_OF_GENERATION", "Number of generations for convergence (default: 5000)", LAMBDA_PARSE_ARG {opts.ngen = atoi(arg[0]); }},
    { {"--max-gen"}, "MAX_GENERATION", "Limits of trial generations (default: 1000000)", LAMBDA_PARSE_ARG {opts.max_gen = atoll(arg[0]); }},
    { {"--npop"}, "NUM_OF_POP", "Candidate sequences size (default: 100)", LAMBDA_PARSE_ARG {opts.npop = atoi(arg[0]); }},
};
