#pragma once

#include "Types.h"
#include "MmcModel.h"

void mmc_s_function_start(int N, real_t R, real_t Rc, real_t L, real_t Lc, real_t C, real_t Ts, real_t freq, real_t n_sm, real_t q1, real_t q2, const real_t x_min[N_X], const real_t x_max[N_X], real_t insertion_index_deviation_allowance, const real_t u_min[N_U], const real_t u_max[N_U]);

void mmc_s_function_terminate(void);

int mmc_s_function(int N, real_t phase_Vf, real_t phase_Iv, real_t amp_Vf, real_t amp_Iv, real_t Vdc, real_t Icir_ref, const real_t x0[N_X], real_t x[N][N_X], real_t u[N][N_U]);
