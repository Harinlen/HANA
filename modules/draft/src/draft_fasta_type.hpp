#ifndef DRAFT_FASTA_TYPE_H
#define DRAFT_FASTA_TYPE_H

#include <cstdint>

typedef struct ENZYME_RANGE
{
    int32_t start, end;
} ENZYME_RANGE;

typedef struct ENZYME_RANGES
{
    ENZYME_RANGE* ranges;
    size_t length, counter;
} ENZYME_RANGES;

#endif // DRAFT_FASTA_TYPE_H