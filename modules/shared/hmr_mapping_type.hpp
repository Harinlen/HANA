#ifndef HMR_MAPPING_TYPE_H
#define HMR_MAPPING_TYPE_H

#include <cstdint>

/* Critical mapping information from the BAM file */
typedef struct MAPPING_INFO
{
    int32_t refID;
    int32_t pos;
    int32_t next_refID;
    int32_t next_pos;
    uint16_t flag;
    uint8_t mapq;
} MAPPING_INFO;

typedef void (*MAPPING_N_CONTIG)(uint32_t, void*);
typedef void (*MAPPING_CONTIG)(uint32_t, char*, uint32_t, void*);
typedef void (*MAPPING_READ_ALIGN)(size_t, const MAPPING_INFO &, void*);

typedef struct MAPPING_PROC
{
    MAPPING_N_CONTIG proc_no_of_contig;
    MAPPING_CONTIG proc_contig;
    MAPPING_READ_ALIGN proc_read_align;
} MAPPING_PROC;

#endif // HMR_MAPPING_TYPE_H