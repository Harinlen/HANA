#include <cstdlib>

#include "args_draft.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-n", "--nodes"}, "NODES", "HMR contig node file (.hmr_nodes)", LAMBDA_PARSE_ARG {opts.nodes = arg[0]; }},
    { {"-r", "--reads"}, "READS", "HMR paired-reads file (.hmr_reads)", LAMBDA_PARSE_ARG { opts.reads = arg[0]; }},
    { {"-a", "--allele-table"}, "ALLELE_TABLE", "Allele contig table (.hmr_allele_table)", LAMBDA_PARSE_ARG {opts.allele_table = arg[0]; }},
    { {"-o", "--output"}, "OUTPUT", "Output graph prefix", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
    { {"-b", "--buffer-size"}, "BUFFER_SIZE", "HMR paired-reads buffer size (unit: K, default: 512)", LAMBDA_PARSE_ARG {opts.read_buffer_size = atoi(arg[0]); }},
    { {"--min-links"}, "MIN_LINKS", "Minimum number of links for contig pair (default: 3)", LAMBDA_PARSE_ARG {opts.min_links = atoi(arg[0]); }},
    { {"--min-re"}, "MIN_RE", "Minimum number of RE sites in a contig (default: 10)", LAMBDA_PARSE_ARG {opts.min_re = atoi(arg[0]); }},
    { {"--max-link-density"}, "MAX_DENSITY", "Maximum allowed link density (default: 2)", LAMBDA_PARSE_ARG {opts.max_density = atof(arg[0]); }},
};
