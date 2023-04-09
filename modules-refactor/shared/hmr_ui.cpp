#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cstdarg>
#include <mutex>
#ifdef _MSC_VER
#include <Windows.h>
#endif

#include "hmr_ui.hpp"

#ifdef _MSC_VER
HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
#define sprintf     sprintf_s
std::mutex ui_mutex;
#endif

void time_print(const char *str, ...)
{
    time_t now = time(0);
    va_list args;
    va_start(args, str);
    char buffer[1025];
    vsnprintf(buffer, 1024, str, args);
    va_end(args);
#ifdef _MSC_VER
    struct tm tstruct;
    localtime_s(&tstruct, &now);
    {
        std::unique_lock<std::mutex> lock(ui_mutex);
        SetConsoleTextAttribute(hStdOut, 2);
        printf("[%02d:%02d:%02d]",
            tstruct.tm_hour, tstruct.tm_min, tstruct.tm_sec);
        SetConsoleTextAttribute(hStdOut, 7);
        printf(" %s\n", buffer);
    }
#else
    struct tm *tstruct = localtime(&now);
    printf("\033[32m[%02d:%02d:%02d]\033[0m %s\n",
           tstruct->tm_hour, tstruct->tm_min, tstruct->tm_sec, buffer);
#endif
}

void time_error(int exitCode, const char *str, ...)
{
    //Print the data.
    va_list args;
    va_start(args, str);
    char buffer[1025];
    vsnprintf(buffer, 1024, str, args);
    va_end(args);
    time_print("%s", buffer);
    //Exit the program.
    exit(exitCode);
}
