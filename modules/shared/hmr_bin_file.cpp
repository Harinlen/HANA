#include "hmr_bin_file.hpp"

bool bin_open(const char* filepath, FILE** file, char const* mode)
{
#ifdef _MSC_VER
    FILE* bin_file = NULL;
    if (fopen_s(&bin_file, filepath, mode) != 0)
    {
        return false;
    }
#else
    FILE* bin_file = fopen(filepath, mode);
    if (bin_file == NULL)
    {
        return false;
    }
#endif
    *file = bin_file;
    return true;
}
