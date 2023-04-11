#include <cstdlib>

#include "args_dump.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-n", "--nodes"}, "NODES", "HMR Nodes file (.hmr_nodes)", LAMBDA_PARSE_ARG {opts.nodes = arg[0]; }},
    { {"-o", "--output"}, "OUTPUT", "Output path", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
};
