#include <cstdlib>

#include "args_ordering.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-n", "--nodes"}, "NDOES", "HMR contig node file (.hmr_contig)", LAMBDA_PARSE_ARG {opts.nodes = arg[0]; }},
    { {"-e", "--edge"}, "EDGE", "HMR edge file (.hmr_edge)", LAMBDA_PARSE_ARG { opts.edge = arg[0];}},
    { {"-g", "--group"}, "GROUP", "HMR contig group file (.hmr_group)", LAMBDA_PARSE_ARG { opts.group = arg[0];}},
    { {"-o", "--output"}, "OUTPUT", "Output partition file (.hmr_partition)", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
};
