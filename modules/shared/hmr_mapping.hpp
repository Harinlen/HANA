#ifndef HMR_MAPPING_H
#define HMR_MAPPING_H

#include <cstdlib>

#include "hmr_mapping_type.hpp"

void hmr_mapping_read(const char* filepath, MAPPING_PROC proc, void *user, int threads);

#endif // HMR_MAPPING_H