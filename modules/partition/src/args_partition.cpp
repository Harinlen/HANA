#include <cstdlib>

#include "args_partition.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-n", "--nodes"}, "NDOES", "HMR contig node file (.hmr_contig)", LAMBDA_PARSE_ARG {opts.nodes = arg[0]; }},
    { {"-e", "--edge"}, "EDGE", "HMR edge file (.hmr_edge)", LAMBDA_PARSE_ARG { opts.edge = arg[0];}},
    { {"-g", "--group"}, "GROUP", "Number of homologous chromosomes groups", LAMBDA_PARSE_ARG {opts.groups = atoi(arg[0]); }},
    { {"-r", "--re"}, "MIN_RE", "Number of minimum REs (default: 25)", LAMBDA_PARSE_ARG {opts.min_re = atoi(arg[0]); }},
    { {"--max-link-density"}, "MAX_LINK_DENSITY", "Density threshold of repetitive contig (default: 2)", LAMBDA_PARSE_ARG {opts.max_link_density = atoi(arg[0]); }},
    { {"-a", "--table"}, "ALLELE_TABLE", "Allele contig table (.hmr_allele/.ctg.table)", LAMBDA_PARSE_ARG {opts.allele_table = arg[0]; }},
    { {"-o", "--output"}, "OUTPUT", "Output partition file (.hmr_partition)", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
};
