#ifndef EXTRACT_INDEX_MAP_H
#define EXTRACT_INDEX_MAP_H

#include <unordered_map>

typedef std::unordered_map<std::string, int32_t> CONTIG_INDEX_MAP;

int32_t extract_contig_index_get(const CONTIG_INDEX_MAP& contig_name_map, const std::string& contig_name);

#endif // EXTRACT_INDEX_MAP_H