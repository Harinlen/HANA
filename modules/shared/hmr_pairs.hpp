#ifndef HMR_PAIRS_H
#define HMR_PAIRS_H

#include <cstdlib>

typedef void (*PAIR_PROC)(const char *ref, size_t ref_len, const char *next_ref, size_t next_ref_len,
                          int32_t pos, int32_t next_pos, const char *types, void *user);
void hmr_pairs_read(const char *filepath, PAIR_PROC proc, void *user);

#endif // HMR_PAIRS_H
