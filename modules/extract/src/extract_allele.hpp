#ifndef EXTRACT_ALLELE_H
#define EXTRACT_ALLELE_H

#include "hmr_contig_graph_type.hpp"

#include "extract_index_map.hpp"

HMR_CONTIG_ID_TABLE extract_allele_table(const char *filepath, CONTIG_INDEX_MAP *index_map);

#endif // EXTRACT_ALLELE_H
