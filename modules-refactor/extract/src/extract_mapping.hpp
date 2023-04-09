#ifndef EXTRACT_MAPPING_H
#define EXTRACT_MAPPING_H

#include "extract_mapping_type.hpp"

void extract_mapping_file(const char *filepath, CONTIG_INDEX_MAP* index_map, FILE* reads_file, CONTIG_ENZYME_RANGES* contig_enzyme_ranges, uint8_t mapq, int32_t buffer_size, int32_t threads);

#endif // EXTRACT_MAPPING_H