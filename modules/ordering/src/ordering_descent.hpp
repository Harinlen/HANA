#ifndef ORDERING_DESCENT_H
#define ORDERING_DESCENT_H

#include <random>

#include "ordering_type.hpp"

void ordering_descent_init(const HMR_CONTIG_ID_VEC& contig_group, const HMR_NODES& contigs, ORDERING_INFO& info);
HMR_CONTIG_ID_VEC ordering_descent_optimize(int32_t phase, int npop, int ngen, uint64_t maxgen, double mutapb, ORDERING_INFO& info, std::mt19937_64& rng, int threads);

#endif // ORDERING_DESCENT_H
