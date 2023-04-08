#include <cstdlib>
#include <cstring>
#include <cstdarg>
#include <unordered_map>
#include <algorithm>

#include "hmr_args_types.hpp"

#include "hmr_args.hpp"

extern HMR_ARGS opts;
extern HMR_ARG_PARSER args_parser;

typedef struct PARSE_RULE
{
    bool has_param;
    ARG_PARSE parser;
} PARSE_RULE;

typedef std::unordered_map<std::string, PARSE_RULE> PARSE_MAP;

void help_exit(int exit_code, const char* msg, ...)
{
    //Print the option message.
    if (msg != NULL)
    {
        //Generate the exit message.
        va_list args;
        va_start(args, msg);
        char buffer[1025];
        vsnprintf(buffer, 1024, msg, args);
        va_end(args);
        //Print the message.
        printf("%s\n", buffer);
    }
    //Print the help information.
    printf("usage: [-h]");
    for (const HMR_ARG_INFO& info : args_parser)
    {
        std::string min_length_opt;
        for (const auto& opt : info.opts)
        {
            if (min_length_opt.empty() || opt.size() < min_length_opt.size())
            {
                min_length_opt = opt;
            }
        }
        printf(" %s", min_length_opt.data());
        if (!info.val_name.empty())
        {
            printf(" %s", info.val_name.data());
        }
    }
    printf("\noptional arguments:\n"
        "  -h, --help            Show this help message and exit\n");
    for (const HMR_ARG_INFO& info : args_parser)
    {
        const char* val_name = info.val_name.c_str();
        //Prepare the line of option line.
        char arg_buf[1024], arg_item[1024];
        arg_buf[0] = '\0';
        bool not_first = false;
        std::vector<std::string> arg_names(info.opts.begin(), info.opts.end());
        std::sort(arg_names.begin(), arg_names.end(),
            [](std::string& lhs, std::string& rhs) {
                return lhs.length() < rhs.length();
            });
        //Check whether we have option.
        if (info.val_name.empty())
        {
            //No other otpions.
            for (const auto& i : arg_names)
            {
                //Check the format.
                if (not_first)
                {
#ifdef _MSC_VER
                    sprintf_s(arg_item, 1023, ", %s", i.c_str());
#else
                    sprintf(arg_item, ", %s", i.c_str());
#endif
                }
                else
                {
#ifdef _MSC_VER
                    sprintf_s(arg_item, 1023, "%s", i.c_str());
#else
                    sprintf(arg_item, "%s", i.c_str());
#endif
                }
                not_first = true;
#ifdef _MSC_VER
                strcat_s(arg_buf, 1024, arg_item);
#else
                strcat(arg_buf, arg_item);
#endif
            }
        }
        else
        {
            for (const auto& i : arg_names)
            {
                if (not_first)
                {
#ifdef _MSC_VER
                    sprintf_s(arg_item, 1023, ", %s %s", i.c_str(), val_name);
#else
                    sprintf(arg_item, ", %s %s", i.c_str(), val_name);
#endif
                }
                else
                {
#ifdef _MSC_VER
                    sprintf_s(arg_item, 1023, "%s %s", i.c_str(), val_name);
#else
                    sprintf(arg_item, "%s %s", i.c_str(), val_name);
#endif
                }
                not_first = true;
#ifdef _MSC_VER
                strcat_s(arg_buf, 1024, arg_item);
#else
                strcat(arg_buf, arg_item);
#endif
            }
        }
        printf("  %s", arg_buf);
        size_t arg_len = strlen(arg_buf);
        if (arg_len < 21)
        {
            //Fill to 22 spacing.
            printf("%*s", static_cast<int>(20 - arg_len), "");
        }
        else
        {
            printf("\n%*s  ", 20, "");
        }
        printf("  %s\n", info.description.c_str());
    }
    exit(exit_code);
}

void parse_arguments(int argc, char* argv[])
{
    //Check argument count.
    if (argc == 1)
    {
        help_exit(0, NULL);
    }
    //Build the parsing map.
    PARSE_MAP parse_map;
    //All the option name is add to the dictionary.
    for (const HMR_ARG_INFO& info : args_parser)
    {
        for (const std::string& opt_name : info.opts)
        {
            parse_map.insert(std::make_pair(opt_name, PARSE_RULE{ !info.val_name.empty(), info.set_val }));
        }
    }
    //Construct the parsing table.
    int arg_id = 1;
    while (arg_id < argc)
    {
        const char* opt = argv[arg_id];
        //Check the prefix of the argc.
        size_t arg_size = strlen(opt);
        if (arg_size < 1 || opt[0] != '-')
        {
            //Unknown option detected.
            help_exit(-1, "Invalid option %s", opt);
        }
        //Check for help.
        std::string opt_str = std::string(opt);
        if (opt_str == "-h" || opt_str == "--help")
        {
            help_exit(0, NULL);
        }
        //Try to find the arguments.
        auto opt_finder = parse_map.find(opt_str);
        if (opt_finder == parse_map.end())
        {
            help_exit(-1, "Unknown option %s", opt);
        }
        //Find all the following parameters.
        std::vector<char*> args;
        if (opt_finder->second.has_param)
        {
            while (arg_id + 1 < argc)
            {
                char* next_arg = argv[arg_id + 1];
                size_t next_arg_size = strlen(next_arg);
                //If the next arguments start with '-', then we stop searching.
                if (next_arg_size < 1 || next_arg[0] == '-')
                {
                    break;
                }
                args.push_back(next_arg);
                ++arg_id;
            }
        }
        //Passing the arguments to the parser.
        opt_finder->second.parser(args);
        //Keep parsing.
        ++arg_id;
    }
}
