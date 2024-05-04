#pragma once

#include <stddef.h>

#include "Types.h"

#define RAMP_ERROR_INFEASIBLE -1
#define RAMP_ERROR_RANK_2_UPDATE -2

#define RAMP_HOTSTART_NONE 0
#define RAMP_HOTSTART_M4_UNCHANGED 1
#define RAMP_HOTSTART_M4_CHANGED 2

void ramp_init(size_t n_H, size_t n_a, void (*get_column_M4)(size_t, real_t*));

void ramp_cleanup(void);

void ramp_enable_infeasibility_error(real_t min, real_t max);

void ramp_disable_infeasibility_error(void);

int ramp_solve(size_t n_H, size_t n_a, int hotstart_variant, const real_t m[n_H], real_t y[n_H]);
