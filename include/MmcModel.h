#pragma once

#include <stddef.h>

#include "Types.h"

#define N_X 4
#define N_U 2

void mmc_model_get_init(real_t R, real_t Rc, real_t L, real_t Lc, real_t C, real_t Ts, real_t n_sm, size_t N, real_t A[N][N_X][N_X], real_t B[N][N_X][N_U], real_t d[N][N_X]);

void mmc_model_get(size_t N, const real_t x[N][N_X], const real_t u[N][N_U], const real_t vf[N], real_t Vdc, real_t A[N][N_X][N_X], real_t B[N][N_X][N_U], real_t d[N][N_X]);
