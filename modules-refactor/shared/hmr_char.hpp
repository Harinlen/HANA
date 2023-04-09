#ifndef HMR_CHAR_H
#define HMR_CHAR_H

#include <string>

/* Check whether the input char is a spacing char */
inline int is_space(char c)
{
    return c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r' || c == ' ';
}

/* Trim the string right to the none spacing char */
inline void trimmed_right(char* buf, size_t& size)
{
    while (size > 0 && is_space(buf[size - 1]))
    {
        buf[--size] = '\0';
    }
}

#endif // HMR_CHAR_H