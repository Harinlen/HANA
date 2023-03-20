#include <algorithm>
#include <cctype>
#include <cstring>
#include <string>
#include <unordered_map>

#include "hmr_ui.hpp"
#include "hmr_bin_file.hpp"
#include "hmr_seq.hpp"

#include "hmr_enzyme.hpp"

typedef struct ENZYME_SEQ
{
    std::vector<std::string> names;
    std::string sequence;
} ENZYME_SEQ;

static std::vector<ENZYME_SEQ> known_enzymes = {
    { {"TaqI", "Taq1"}, "TCGA"},
    { {"HaeIII", "Hae3"}, "GGCC"},
    { {"DpnI", "Dpn1", "DpnII", "Dpn2", "MboI", "Mbo1", "Sau3AI"}, "GATC"},
    { {"AluI", "Alu1"}, "AGCT"},
    { {"NlaIII", "Nla3", "FaeI", "Fae1", "FatI", "Fat1", "Hin1II", "Hsp92II"}, "CATG"},
    { {"HpaII", "Hpa2"}, "CCGG"},
    { {"FokI", "Fok1"}, "GGATG"},
    { {"AaaI", "Aaa1"}, "CGGCG"},
    { {"HgaI", "Hga1"}, "GACGC"},
    { {"BglII", "Bgl2"}, "AGATCT"},
    { {"EcoRV", "EcoR5"}, "GATATC"},
    { {"EcoRI", "EcoR1"}, "GAATTC"},
    { {"BamHI", "BamH1"}, "GGATCC"},
    { {"HindIII", "Hind3"}, "AAGCTT"},
    { {"KpnI", "Kpn1"}, "GGTACC"},
    { {"XbaI", "Xba1"}, "TCTAGA"},
    { {"XhoI", "Xho1"}, "CTCGAG"},
    { {"SacI", "Sac1"}, "GAGCTC"},
    { {"PstI", "Pst1"}, "CTGCAG"},
    { {"SmaI", "Sma1"}, "CCCGGG"},
    { {"PvuII", "Pvu2"}, "CAGCTG"},
    { {"SalI", "Sal1"}, "GTCGAC"},
    { {"ScaI", "Sca1"}, "AGTACT"},
    { {"SpeI", "Spe1"}, "ACTAGT"},
    { {"SphI", "Sph1"}, "GCATGC"},
    { {"StuI", "Stu1"}, "AGGCCT"},
    { {"NdeI", "Nde1"}, "CATATG"},
    { {"NotI", "Not1"}, "GCGGCCGC"},
};

static std::unordered_map<std::string, std::string> known_enzyme_alias;

void hmr_enzyme_formalize(char* enzyme, const char** nuc_seq, int* nuc_seq_size)
{
    size_t length = strlen(enzyme);
    //Convert the original char in upper case letter.
    hmr_seq_upper(enzyme, length);
    //Check whether the alias dictionary is built.
    if (known_enzyme_alias.empty())
    {
        //Build the dictionary.
        for (ENZYME_SEQ& enzyme_info : known_enzymes)
        {
            for (auto& enzyme_name : enzyme_info.names)
            {
                //Convert the name into upper case.
                std::transform(enzyme_name.begin(), enzyme_name.end(), enzyme_name.begin(), 
                    [](char c) {return std::toupper(c); });
                known_enzyme_alias.insert(std::make_pair(enzyme_name, enzyme_info.sequence));
            }
        }
    }
    //Check whether it is an known alias name.
    auto known_finder = known_enzyme_alias.find(std::string(enzyme));
    if (known_finder != known_enzyme_alias.end())
    {
        (*nuc_seq) = known_finder->second.data();
        (*nuc_seq_size) = static_cast<int>(strlen(*nuc_seq));
        return;
    }
    //Or else we have to check the enzyme is valid or not.
    for (char* s = enzyme, *e = enzyme + length; s < e; ++s)
    {
        //Check invalid nuc.
        if ((*s) != 'A' && (*s) != 'C' && (*s) != 'T' && (*s) != 'G')
        {
            time_error(-1, "Invalid nucleotide base '%c' found in enzyme '%s'.", *s, enzyme);
        }
    }
    *nuc_seq = enzyme;
    *nuc_seq_size = static_cast<int>(length);
}
