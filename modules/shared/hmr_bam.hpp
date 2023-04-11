#ifndef HMR_BAM_H
#define HMR_BAM_H

#include <cstdint>
#include <cstdlib>

/* Critical mapping information from the BAM file */
typedef struct BAM_BLOCK_HEADER
{
    int32_t refID;
    int32_t pos;
    uint8_t l_read_name;
    uint8_t mapq;
    uint16_t bin;
    uint16_t n_cigar_op;
    uint16_t flag;
    uint32_t l_seq;
    int32_t next_refID;
    int32_t next_pos;
    int32_t tlen;
    void* data;
} BAM_BLOCK_HEADER;

typedef void (*BAM_N_CONTIG)(uint32_t num_of_contigs, void* user);
typedef void (*BAM_CONTIG)(uint32_t name_length, char* name, uint32_t , void* user);
typedef void (*BAM_READ_ALIGN)(size_t block_id, const BAM_BLOCK_HEADER *bam_block, void* user);

typedef struct BAM_MAPPING_PROC
{
    BAM_N_CONTIG proc_no_of_contig;
    BAM_CONTIG proc_contig;
    BAM_READ_ALIGN proc_read_align;
} BAM_MAPPING_PROC;

void hmr_bam_read(const char *filepath, BAM_MAPPING_PROC proc, void* user, int threads);

#endif // HMR_BAM_H
