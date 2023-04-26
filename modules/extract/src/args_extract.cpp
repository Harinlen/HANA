#include <cstdlib>

#include "args_extract.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-f", "--fasta"}, "FASTA", "Contig FASTA file (.fasta/.fasta.gz)", LAMBDA_PARSE_ARG {opts.fasta = arg[0]; }},
    { {"-a", "--allele"}, "ALLELE", "Contig allele table file (.ctg.table)", LAMBDA_PARSE_ARG {opts.allele = arg[0]; }},
    { {"-m", "--mapping"}, "MAPPING 1, MAPPING 2...", "Hi-C reads mapping files (.bam/.pairs)", LAMBDA_PARSE_ARG { opts.mappings = arg; }},
    { {"-o", "--output"}, "OUTPUT", "Output file prefix", LAMBDA_PARSE_ARG {opts.output = arg[0]; }},
    { {"-e", "--enzyme"}, "ENZYME", "Enzyme to find in the sequence", LAMBDA_PARSE_ARG {opts.enzyme = arg[0]; }},
    { {"-t", "--threads"}, "THREADS", "Number of threads (default: 1)", LAMBDA_PARSE_ARG { opts.threads = atoi(arg[0]); }},
    { {"-r", "--range"}, "ENZYME_RANGE", "Hi-C pairs valid distance around enzyme (default: 1000)", LAMBDA_PARSE_ARG {opts.range = atoi(arg[0]) >> 1; }},
    { {"-q", "--mapq"}, "MAPQ", "BAM minimum mapping quality (default: 40)", LAMBDA_PARSE_ARG {opts.mapq = atoi(arg[0]); }},
    { {"--search-buffer"}, "SEARCH_BUF_SIZE", "FASTA searching buffer size per thread (default: 32)", LAMBDA_PARSE_ARG { opts.fasta_pool = atoi(arg[0]); }},
    { {"--pairs-read-len"}, "PAIRS_READ_LEN", "The length used for pairs file reads (default: 150)", LAMBDA_PARSE_ARG { opts.pairs_read_len = atoi(arg[0]); }},
    { {"--mapping-buffer"}, "MAP_BUF_SIZE", "Mapping parse buffer size (unit: K, default: 512)", LAMBDA_PARSE_ARG { opts.mapping_pool = atoi(arg[0]); }},
    { {"--no-flag"}, "", "Skip the flag checking", LAMBDA_PARSE_ARG { (void)arg; opts.skip_flag = true; }},
    { {"--no-range"}, "", "Skip the enzyme range checking", LAMBDA_PARSE_ARG { (void)arg; opts.skip_range = true; }},
};
