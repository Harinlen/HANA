#include <cstring>

#include "hmr_bin_file.hpp"

#include "hmr_path.hpp"

void path_split(const char* filepath, size_t length, std::string& filepath_base, std::string& filepath_suffix)
{
    //When the length is 0, auto detect the length of the path.
    if (length == 0)
    {
        length = strlen(filepath);
    }
    std::string filepath_str(filepath, length);
    //Find the dot from the last one.
    size_t dot_pos = filepath_str.find_last_of('.');
    if (dot_pos == std::string::npos)
    {
        //Failed to find the position, the base is the full filepath.
        filepath_base = filepath_str;
        filepath_suffix = std::string();
        return;
    }
    //Set the result.
    filepath_base = filepath_str.substr(0, dot_pos);
    filepath_suffix = filepath_str.substr(dot_pos);
}

bool path_ends_with(const char* filepath, const char* suffix)
{
    size_t path_len = strlen(filepath);
    size_t suffix_len = strlen(suffix);
    if (path_len < suffix_len)
    {
        return false;
    }
    //Compare the path suffix len.
    for (const char* f_p = filepath + (path_len - suffix_len), *s_p = suffix, *s_max = suffix + suffix_len; s_p < s_max; )
    {
        if (*f_p != *s_p)
        {
            return false;
        }
        ++f_p;
        ++s_p;
    }
    return true;
}

std::string path_suffix(const char *filepath, size_t length)
{
    //Find the dot from the last one.
    std::string filepath_base, filepath_suffix;
    //Split the filepath.
    path_split(filepath, length, filepath_base, filepath_suffix);
    //Get the suffix of the file.
    return filepath_suffix;
}

std::string path_basename_core(const char* filepath, size_t length)
{
    //Find the dot from the last one.
    std::string filepath_base, filepath_suffix;
    //Split the filepath.
    path_split(filepath, length, filepath_base, filepath_suffix);
    //Get the basename of the file.
    return filepath_base;
}

std::string path_basename(const char* filepath, size_t length)
{
    //If the suffix is gz, remove the gz.
    if (path_suffix(filepath, length) == ".gz")
    {
        return path_basename_core(filepath, strlen(filepath) - 3);
    }
    //Just get the basename.
    return path_basename_core(filepath, length);
}

bool path_can_read(const char* filepath)
{
    //Try to open the file to have a test.
    FILE* file_test;
    bool result = bin_open(filepath, &file_test, "rb");
    //If we open it successfully, close it right now.
    if (result)
    {
        fclose(file_test);
    }
    return result;
}
