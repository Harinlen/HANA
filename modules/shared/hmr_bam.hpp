#ifndef HMR_BAM_H
#define HMR_BAM_H

#include "hmr_mapping_type.hpp"

void hmr_bam_read(const char *filepath, MAPPING_PROC proc, void* user, int threads);

#endif // HMR_BAM_H