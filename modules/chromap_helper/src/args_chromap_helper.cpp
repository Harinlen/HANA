#include <cstdlib>

#include "args_chromap_helper.hpp"

#include "hmr_args_types.hpp"

HMR_ARGS opts;

HMR_ARG_PARSER args_parser = {
    { {"-q", "--fastq"}, "FASTQ", "FASTQ reads file (.fastq.gz/.fastq)", LAMBDA_PARSE_ARG { opts.reads = arg[0]; }},
};
