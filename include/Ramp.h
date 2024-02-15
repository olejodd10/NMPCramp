#pragma once

#include <stddef.h>

#include "Types.h"
#include "IterableSet.h"
#include "IndexedVectors.h"

#define RAMP_ERROR_INFEASIBLE -1
#define RAMP_ERROR_RANK_2_UPDATE -2

void ramp_init(size_t n_H, void (*get_column_m4)(size_t, real_t*));

void ramp_cleanup(void);

void ramp_enable_infeasibility_error(real_t min, real_t max);

void ramp_disable_infeasibility_error(void);

int ramp_solve(size_t n_H, size_t n_a, iterable_set_t *a_set, indexed_vectors_t *invq, real_t y[n_H]);
