#pragma once

#include <stddef.h>

#include "Types.h"

void mmc_mpc_init(
        size_t n_x,
        size_t n_u,
        size_t N,

        real_t q1,
        real_t q2,

        const real_t x_min[n_x],
        const real_t x_max[n_x],
        real_t n_sm,
        real_t insertion_index_deviation_allowance,
        const real_t u_min[n_u],
        const real_t u_max[n_u]
);

void mmc_mpc_cleanup(void);

int mmc_mpc_solve(size_t n_x, size_t n_u, size_t N, 
        const real_t x1_ref[N], real_t x2_ref, 
        const real_t A[N][n_x][n_x], const real_t B[N][n_x][n_u], const real_t d[N][n_x], const real_t x0[n_x], 
        real_t x[N][n_x], real_t u[N][n_u]);
