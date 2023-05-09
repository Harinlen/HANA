#include "hmr_seq.hpp"

void hmr_seq_upper(char *seq, size_t seq_len)
{
    //Convert the original char in upper case letter.
    for (char* s = seq, *e = seq + seq_len; s < e; ++s)
    {
        if ((*s) >= 'a' && (*s) <= 'z')
        {
            (*s) = 'A' + ((*s) - 'a');
        }
    }
}

bool hmr_seq_valid(const char* seq, size_t seq_len)
{
    //Seq must be upper case, and should be one of A, T, G, C, U.
    for (const char* s = seq, *e = seq + seq_len; s < e; ++s)
    {
        if ((*s) != 'A' && (*s) != 'T' && (*s) != 'C' && (*s) != 'G' && (*s) != 'U')
        {
            return false;
        }
    }
    return true;
}
