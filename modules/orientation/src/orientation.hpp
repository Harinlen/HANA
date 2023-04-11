#ifndef ORIENTATION_H
#define ORIENTATION_H

#include "hmr_contig_graph_type.hpp"

constexpr auto DIRECTION_POSITIVE = 0;
constexpr auto DIRECTION_NEGATIVE = 1;

typedef struct ORIENTATION_SEQUENCE
{
    HMR_CONTIG_ID_VEC contig_id_seq;
    std::unordered_map<int32_t, int32_t> contig_id_map;
    double* contig_size;
    double* cost_buffer;
    double* cost_matrix[2];
} ORIENTATION_SEQUENCE;

typedef struct ORIENTATION_INFO
{
    std::vector<ORIENTATION_SEQUENCE> sequences;
    std::vector<int32_t> belongs;
} ORIENTATION_INFO;

void orientation_init(std::vector<char*> seq_paths, const HMR_NODES& nodes, ORIENTATION_INFO& info);
void orientation_calc_gradient(HMR_MAPPING* mapping, int32_t buf_size, void* user);
CHROMOSOME_CONTIGS orientation_extract(const ORIENTATION_SEQUENCE& sequence);

#endif // ORIENTATION_H