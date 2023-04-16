#ifndef EXTRACT_FASTA_TYPE_H
#define EXTRACT_FASTA_TYPE_H

#include <vector>
#include <cstdint>

typedef struct ENZYME_RANGE
{
    int32_t start, end;
} ENZYME_RANGE;

typedef std::vector<ENZYME_RANGE> CONTIG_ENZYME_RANGE;
typedef std::vector<CONTIG_ENZYME_RANGE> CONTIG_ENZYME_RANGES;

#endif // EXTRACT_FASTA_TYPE_H
