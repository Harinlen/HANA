#ifndef HMR_BIN_ENZYME_H
#define HMR_BIN_ENZYME_H

#include <vector>
#include <string>

typedef std::vector<std::string> ENZYME_VEC;

ENZYME_VEC hmr_enzyme_formalize(std::vector<char*> &enzyme);

#endif // HMR_BIN_ENZYME
