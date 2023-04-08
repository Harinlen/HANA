#ifndef HMR_ALLELE_H
#define HMR_ALLELE_H

#include "hmr_allele_type.hpp"

typedef void (*ALLELE_PROC)(const CONTIG_ID_VECTOR &, void *);
void hmr_allele_table_load(const char *filepath, const HMR_CONTIGS &contig_table, ALLELE_PROC proc, void *user);
void hmr_allele_table_conflict_set(const char* filepath, const HMR_CONTIGS& contig_table, HMR_ALLELE_TABLE &table);

#endif // HMR_ALLELE_H