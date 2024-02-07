#pragma once

#include <stddef.h>

#include "Types.h"

void mmc_trajectory_shift(size_t n_u, size_t N, const real_t u_opt[N][n_u], real_t u[N][n_u]);
