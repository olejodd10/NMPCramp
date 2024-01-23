#pragma once

#include <stddef.h>

#include "Types.h"

void sdqp_lmpc_constant_init(
        size_t n_x,
        size_t n_u,
        size_t n_y,
        size_t n_t,
        size_t N,

        const real_t Q[n_x][n_x],
        const real_t S[n_x][n_x],
        const real_t R[n_u][n_u],
        const real_t fx[n_x],
        const real_t fu[n_u],

        const real_t A[n_x][n_x],
        const real_t B[n_x][n_u],
        const real_t C[n_y][n_x],

        real_t y_min,
        real_t y_max,
        const real_t Lt[n_t][n_x],
        const real_t lt[n_t],
        real_t u_min,
        real_t u_max); 

void sdqp_lmpc_constant_cleanup(void);

int sdqp_lmpc_constant_solve(
        size_t n_x, size_t n_u, size_t N, 
        const real_t x0[n_x], 
        real_t x[n_x*N], real_t u[n_u*N]);
