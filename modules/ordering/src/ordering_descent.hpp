#ifndef ORDERING_DESCENT_H
#define ORDERING_DESCENT_H

#include "ordering_type.hpp"

struct ORDERING_TIG;

ORDERING_TIG *ordering_init(ORDERING_INFO &info, int32_t *order = NULL);

void ordering_optimize_phase(int32_t phase,
                             int npop,
                             int ngen,
                             double mutapb,
                             ORDERING_INFO &info);

#endif // ORDERING_DESCENT_H
