#ifndef EXTRACT_MAPPING_H
#define EXTRACT_MAPPING_H

#include "extract_mapping_type.hpp"

constexpr auto CHECK_FLAG_RANGE = 1 << 0;
constexpr auto CHECK_FLAG_FLAG = 1 << 1;

void extract_mapping_file(const char* filepath, CONTIG_INDEX_MAP* index_map, FILE* reads_file, CONTIG_ENZYME_RANGES* contig_enzyme_ranges, uint16_t check_flag, uint8_t mapq, int32_t thread_buffer_size, int32_t threads);

#endif // EXTRACT_MAPPING_H