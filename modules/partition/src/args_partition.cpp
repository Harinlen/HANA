#include <cstdlib>

#include "args_partition.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-n", "--nodes"}, "NODES", "HMR contig node file (.hmr_nodes)", LAMBDA_PARSE_ARG {opts.nodes = arg[0]; }},
    { {"-e", "--edges"}, "EDGES", "HMR edges file (.hmr_edges)", LAMBDA_PARSE_ARG { opts.edges = arg[0];}},
    { {"-a", "--allele"}, "ALLELE_TABLE", "HMR allele table file (.hmr_allele_table)", LAMBDA_PARSE_ARG { opts.allele = arg[0]; }},
    { {"-g", "--group"}, "GROUP", "Number of groups to be partitioned", LAMBDA_PARSE_ARG {opts.groups = atoi(arg[0]); }},
    { {"-o", "--output"}, "OUTPUT", "Output partition file (.hmr_partition)", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
    { {"-b", "--buffer-size"}, "BUFFER_SIZE", "HMR edge buffer size (unit: K, default: 512)", LAMBDA_PARSE_ARG {opts.read_buffer_size = atoi(arg[0]); }},
    { {"--non-informative-ratio"}, "NON_INFO_RATIO", "Skipped contigs recover cutoff (default: 3)", LAMBDA_PARSE_ARG {opts.non_informative_ratio = atoi(arg[0]); }},
};
