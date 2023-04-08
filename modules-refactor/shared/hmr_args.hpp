#ifndef HMR_ARGS_H
#define HMR_ARGS_H

/* Display the help information and exit the program */
void help_exit(int exit_code, const char* msg, ...);

/* Parse the CPP arguments */
void parse_arguments(int argc, char* argv[]);

#endif // HMR_ARGS_H