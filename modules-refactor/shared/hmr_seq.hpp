#ifndef HMR_SEQ_H
#define HMR_SEQ_H

#include <cstdint>
#include <cstddef>

void hmr_seq_upper(char *seq, size_t seq_len);
bool hmr_seq_valid(char *seq, size_t seq_len);

#endif // HMR_SEQ_H
