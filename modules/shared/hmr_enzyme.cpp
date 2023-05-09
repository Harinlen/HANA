#include <cstring>
#include <set>

#include "hmr_ui.hpp"
#include "hmr_seq.hpp"

#include "hmr_enzyme.hpp"

//void hmr_enzyme_formalize(char* enzyme, const size_t& enzyme_size)
//{
//    //Convert the original char in upper case letter.
//    hmr_seq_upper(enzyme, enzyme_size);
//    //Check whether the alias dictionary is built.
//    if (!hmr_seq_valid(enzyme, enzyme_size))
//    {
//        time_error(-1, "Invalid restriction sites '%s'.", enzyme);
//    }
//}

ENZYME_VEC hmr_render_enzyme(const std::string &raw, size_t pos, const std::string &target)
{
    ENZYME_VEC enzymes;
    enzymes.reserve(target.size());
    for(size_t i=0; i<target.size(); ++i)
    {
        std::string result = raw;
        result[pos] = target[i];
        enzymes.push_back(result);
    }
    return enzymes;
}

ENZYME_VEC hmr_expand_enzyme(const std::string &enzyme)
{
    const char *enzyme_data = enzyme.data();
    for(size_t i=0; i<enzyme.size(); ++i)
    {
        if(enzyme_data[i] == 'R') { return hmr_render_enzyme(enzyme, i, "AG"); }
        if(enzyme_data[i] == 'Y') { return hmr_render_enzyme(enzyme, i, "CT"); }
        if(enzyme_data[i] == 'S') { return hmr_render_enzyme(enzyme, i, "CG"); }
        if(enzyme_data[i] == 'W') { return hmr_render_enzyme(enzyme, i, "AT"); }
        if(enzyme_data[i] == 'K') { return hmr_render_enzyme(enzyme, i, "GT"); }
        if(enzyme_data[i] == 'M') { return hmr_render_enzyme(enzyme, i, "AC"); }
        if(enzyme_data[i] == 'B') { return hmr_render_enzyme(enzyme, i, "CGT"); }
        if(enzyme_data[i] == 'D') { return hmr_render_enzyme(enzyme, i, "AGT"); }
        if(enzyme_data[i] == 'H') { return hmr_render_enzyme(enzyme, i, "ACT"); }
        if(enzyme_data[i] == 'V') { return hmr_render_enzyme(enzyme, i, "ACG"); }
        if(enzyme_data[i] == 'N') { return hmr_render_enzyme(enzyme, i, "ACGT"); }
    }
    return ENZYME_VEC();
}

ENZYME_VEC hmr_enzyme_formalize(std::vector<char *> &enzymes)
{
    //Convert the original char in upper case letter
    for(char *enzyme_str: enzymes)
    {
        hmr_seq_upper(enzyme_str, strlen(enzyme_str));
    }
    //Construct the enzyme vector.
    ENZYME_VEC enzyme_vec;
    for(char *enzyme_str: enzymes)
    {
        enzyme_vec.push_back(std::string(enzyme_str));
    }
    //Expand the enzyme sites.
    bool vec_is_changed = true;
    while(vec_is_changed)
    {
        //Set the vector changed flag to false.
        vec_is_changed = false;
        //Check each string.
        for(size_t i=0; i<enzyme_vec.size(); ++i)
        {
            //Check whether we can expand the enzyme.
            ENZYME_VEC result = hmr_expand_enzyme(enzyme_vec[i]);
            if(!result.empty())
            {
//                for(const auto &k: enzyme_vec)
//                {
//                    printf("%s\n", k.c_str());
//                }
//                printf("++++++++++\n");
//                for(const auto &k: result)
//                {
//                    printf("%s\n", k.c_str());
//                }
//                printf("----------\n");
                //Reset the flag.
                vec_is_changed = true;
                //Reconstruct the enzyme list.
                enzyme_vec.reserve(enzyme_vec.size() + result.size());
                enzyme_vec.erase(enzyme_vec.begin() + i);
//                for(const auto &k: enzyme_vec)
//                {
//                    printf("%s\n", k.c_str());
//                }
//                printf("**********\n");
                enzyme_vec.insert(enzyme_vec.begin() + i, result.begin(), result.end());
//                for(const auto &k: enzyme_vec)
//                {
//                    printf("%s\n", k.c_str());
//                }
//                exit(1);
                break;
            }
        }
    }
    //Make sure each enzyme only appears once.
    {
//        std::set<std::string> enzyme_set(enzyme_vec.begin(), enzyme_vec.end());
//        enzyme_vec = ENZYME_VEC(enzyme_set.begin(), enzyme_set.end());
    }
    //Check enzyme validation.
    for (const std::string &enzyme_str: enzyme_vec)
    {
        if (!hmr_seq_valid(enzyme_str.data(), enzyme_str.size()))
        {
            time_error(-1, "Invalid restriction sites '%s'.", enzyme_str.data());
        }
    }
    return enzyme_vec;
}
