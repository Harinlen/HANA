#ifndef HMR_ARGS_TYPES_H
#define HMR_ARGS_TYPES_H

#include <string>
#include <set>
#include <vector>

typedef struct HMR_ARGS HMR_ARGS;

#define PARSE_ARGS  (const std::vector<char *> &arg)

typedef void(*ARG_PARSE)PARSE_ARGS;

#define PARSE_ARG(x)        void (x)PARSE_ARGS
#define LAMBDA_PARSE_ARG    []PARSE_ARGS

typedef struct HMR_ARG_INFO
{
    std::set<std::string> opts;
    std::string val_name;
    std::string description;
    ARG_PARSE set_val;
} HMR_ARG_INFO;

typedef std::vector<HMR_ARG_INFO> HMR_ARG_PARSER;

#endif // HMR_ARGS_TYPES_H
