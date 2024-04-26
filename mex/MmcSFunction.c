#include "MmcSFunction.h"

#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "Types.h"
#include "Utils.h"
#include "MmcTrajectory.h"
#include "MmcModel.h"
#include "SdqpLmpcMmc.h"

#define MMC_DELAY 6 // Delay for Vf and Iv_ref
#define TS1_FACTOR 5
#define LINEARIZATIONS 1

// N_X and N_U are imported from MmcModel
#define V_SIGMA_TRAJ_0 300.0
static const size_t X_TRAJ_0[N_X] = {0.0, 0.0, V_SIGMA_TRAJ_0, V_SIGMA_TRAJ_0};

static real_t *m_A;
static real_t *m_B;
static real_t *m_d;

static real_t *m_x_traj;
static real_t *m_u_traj;

static real_t *m_Vf;
static real_t *m_Iv_ref;

static real_t m_omega;

void mmc_s_function_start(int N, real_t R, real_t Rc, real_t L, real_t Lc, real_t C, real_t Ts, real_t freq, real_t n_sm, real_t q1, real_t q2, const real_t x_min[N_X], const real_t x_max[N_X], real_t insertion_index_deviation_allowance, const real_t u_min[N_U], const real_t u_max[N_U]) {
    m_A = (real_t*)malloc(sizeof(real_t)*N*N_X*N_X);
    m_B = (real_t*)malloc(sizeof(real_t)*N*N_X*N_U);
    m_d = (real_t*)malloc(sizeof(real_t)*N*N_X);

    m_x_traj = (real_t*)malloc(sizeof(real_t)*N*N_X);
    m_u_traj = (real_t*)malloc(sizeof(real_t)*N*N_U);

    m_Vf = (real_t*)malloc(sizeof(real_t)*N);
    m_Iv_ref = (real_t*)malloc(sizeof(real_t)*N);

    m_omega = 2.0*M_PI*freq*Ts;

    mmc_model_get_rk4_init(R, Rc, L, Lc, C, Ts, TS1_FACTOR*Ts, n_sm, (size_t)N);    
    sdqp_lmpc_mmc_init(N_X, N_U, (size_t)N, q1, q2, x_min, x_max, n_sm, insertion_index_deviation_allowance, u_min, u_max);

    // Initialize trajectories
    for (int i = 0; i < N; ++i) {
        m_x_traj[i*N_X + 0] = X_TRAJ_0[0];
        m_x_traj[i*N_X + 1] = X_TRAJ_0[1];
        m_x_traj[i*N_X + 2] = X_TRAJ_0[2];
        m_x_traj[i*N_X + 3] = X_TRAJ_0[3];
    }
    for (int i = 0; i < N*N_U; ++i) {
        m_u_traj[i] = n_sm/2.0;
    }
}

void mmc_s_function_terminate(void) {
    sdqp_lmpc_mmc_cleanup();

    free(m_A);
    free(m_B);
    free(m_d);

    free(m_x_traj);
    free(m_u_traj);

    free(m_Vf);
    free(m_Iv_ref);
}

int mmc_s_function(int N, real_t phase_Vf, real_t phase_Iv, real_t amp_Vf, real_t amp_Iv, real_t Vdc, real_t Icir_ref, const real_t x0[N_X], real_t x[N][N_X], real_t u[N][N_U]) {
    // Adjust trajectories
    memcpy(m_x_traj, x0, sizeof(real_t)*N_X);

    // Predict Vf and Iv_ref
    real_t phase = m_omega*MMC_DELAY;
    for (int i = 0; i < N; ++i) {
        m_Vf[i] = amp_Vf*sin(phase + phase_Vf);
        m_Iv_ref[i] = amp_Iv*sin(phase + phase_Iv);
        phase += (i == 0 ? m_omega : (real_t)TS1_FACTOR*m_omega);
    }

    for (int i = 0; i < LINEARIZATIONS; ++i) {
        // Get new model
        mmc_model_get_rk4((size_t)N, CAST_CONST_2D_VLA(m_x_traj, N_X), CAST_CONST_2D_VLA(m_u_traj, N_U), m_Vf, Vdc, CAST_3D_VLA(m_A, N_X, N_X), CAST_3D_VLA(m_B, N_X, N_U), CAST_2D_VLA(m_d, N_X));

        // Solve
        int err = sdqp_lmpc_mmc_solve(N_X, N_U, (size_t)N, m_Iv_ref, Icir_ref, CAST_CONST_3D_VLA(m_A, N_X, N_X), CAST_CONST_3D_VLA(m_B, N_X, N_U), CAST_CONST_2D_VLA(m_d, N_X), x0, x, u);

        // Store solutions for trajectory adjustments in future time steps
        memcpy(&m_x_traj[N_X], x, sizeof(real_t)*(N-1)*N_X);
        memcpy(m_u_traj, u, sizeof(real_t)*N*N_U);

        if (err) {
            return err;
        }
    }

    return 0;
}
