#ifndef ORDERING_DESCENT_H
#define ORDERING_DESCENT_H

#include "ordering_type.hpp"

ORDERING_TIG* ordering_init_alloc(int32_t contig_size);

void ordering_init(ORDERING_INFO& info, const CONTIG_ID_VECTOR& start_order);
void ordering_optimize_phase(int32_t phase, int npop, int ngen, double mutapb,
                             ORDERING_INFO &info, std::mt19937_64& rng, int threads);

#endif // ORDERING_DESCENT_H
