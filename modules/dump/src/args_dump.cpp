#include <cstdlib>

#include "args_dump.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-b", "--bam"}, "BAM FILE", "BAM file to be handled", LAMBDA_PARSE_ARG {opts.bam_file = arg[0]; }},
    { {"-e", "--edge"}, "EDGE FILE", "HMR edge file to be handled (*.hmr_edge)", LAMBDA_PARSE_ARG {opts.edge_file = arg[0]; }},
    { {"-o", "--output"}, "OUTPUT", "Output file path", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
    { {"-t", "--threads"}, "THREAS", "Number of threads (default: 1)", LAMBDA_PARSE_ARG { opts.threads = atoi(arg[0]); }},
};
