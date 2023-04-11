#ifndef HMR_FASTA_H
#define HMR_FASTA_H

#include <cstdint>

/* FASTA file processing function type */
/* Please notice, `name` and `seq` needs to be free by the PROC function */
typedef void (*FASTA_PROC)(int32_t index, char *name, size_t name_size, char *seq, size_t seq_size, void *user);

/* FASTA file parser */
void hmr_fasta_read(const char *filepath, FASTA_PROC parser, void *user);

#endif // HMR_FASTA_H
