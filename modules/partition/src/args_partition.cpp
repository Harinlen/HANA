#include <cstdlib>

#include "args_partition.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-n", "--nodes"}, "NDOES", "HMR contig node file (.hmr_contig)", LAMBDA_PARSE_ARG {opts.nodes = arg[0]; }},
    { {"-e", "--edge"}, "EDGE", "HMR edge file (.hmr_edge)", LAMBDA_PARSE_ARG { opts.edge = arg[0];}},
    { {"-g", "--group"}, "GROUP", "Number of homologous chromosomes groups", LAMBDA_PARSE_ARG {opts.groups = atoi(arg[0]); }},
    { {"-x", "--table"}, "ALLELE_TABLE", "Allele contig table (.hmr_allele/.ctg.table)", LAMBDA_PARSE_ARG {opts.allele_table = arg[0]; }},
    { {"-o", "--output"}, "OUTPUT", "Output partition file (.hmr_partition)", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
    { {"-t", "--threads"}, "THREAS", "Number of threads (default: 1)", LAMBDA_PARSE_ARG { opts.threads = atoi(arg[0]); }},
};
