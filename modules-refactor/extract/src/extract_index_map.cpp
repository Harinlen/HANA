#include "extract_index_map.hpp"

int32_t extract_contig_index_get(const CONTIG_INDEX_MAP& contig_name_map, const std::string& contig_name)
{
    const auto& contig_name_iter = contig_name_map.find(contig_name);
    return contig_name_iter == contig_name_map.cend() ? -1 : contig_name_iter->second;
}