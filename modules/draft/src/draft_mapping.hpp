#ifndef DRAFT_MAPPING_H
#define DRAFT_MAPPING_H

#include "hmr_mapping_type.hpp"

#include "draft_map_type.hpp"

void mapping_draft_filter_worker(MAPPING_FILTER_USER& filter, MAPPING_FILTER_WORKER& worker, const MAPPING_CONTIG_MAP& map);

void mapping_draft_n_contig(uint32_t n_ref, void* user);
void mapping_draft_contig(uint32_t name_length, char* name, uint32_t, void* user);
void mapping_draft_read_align(size_t, const MAPPING_INFO& mapping_info, void* user);

#endif // DRAFT_MAPPING_H