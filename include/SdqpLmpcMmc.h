#pragma once

#include <stddef.h>

#include "Types.h"

void sdqp_lmpc_mmc_init(
        size_t n_x,
        size_t n_u,
        size_t N,

        const real_t Q[n_x][n_x],
        const real_t S[n_x][n_x],
        const real_t R[n_u][n_u],
        const real_t fx[n_x],
        const real_t fu[n_u],

        const real_t x_min[n_x],
        const real_t x_max[n_x],
        real_t n_sm,
        const real_t u_min[n_u],
        const real_t u_max[n_u]
);

void sdqp_lmpc_mmc_cleanup(void);

int sdqp_lmpc_mmc_solve(size_t n_x, size_t n_u, size_t N, 
        const real_t A[N][n_x][n_x], const real_t B[N][n_x][n_u], const real_t d[N][n_x], const real_t x0[n_x], 
        real_t x[N][n_x], real_t u[N][n_u]);
