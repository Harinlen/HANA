#include "hmr_ui.hpp"
#include "hmr_seq.hpp"

#include "hmr_enzyme.hpp"

void hmr_enzyme_formalize(char* enzyme, const size_t& enzyme_size)
{
    //Convert the original char in upper case letter.
    hmr_seq_upper(enzyme, enzyme_size);
    //Check whether the alias dictionary is built.
    if (!hmr_seq_valid(enzyme, enzyme_size))
    {
        time_error(-1, "Invalid restriction sites '%s'.", enzyme);
    }
}