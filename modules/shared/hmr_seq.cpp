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
