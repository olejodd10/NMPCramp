#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Types.h"
#include "LinAlg.h"
#include "MmcTrajectory.h"
#include "Extrapolate.h"
#include "MmcModel.h"
#include "SdqpLmpcMmc.h"
#include "Csv.h"
#include "Simulate.h"
#include "Utils.h"
#include "Timer.h"

// System parameters
#define R 10.0e-3
#define L 1.5e-3
#define C 20.0e-3

#define Rc 0.0
#define Lc 0.0

#define N_SM 18.0

#define FREQ 50.0
#define PHASE_0 0.0

#define Vf_amp 100.0
#define Vf_phase 0.0 // Phase relative to PHASE_0
#define Vdc 300.0

// Costs and references
#define q1 1.0
#define q2 0.3

#define P_ref 7.5e3
#define Idc_ref P_ref/Vdc
#define Icir_ref Idc_ref/3.0

#define Iv_ref_amp 2.0*P_ref/(3.0*Vf_amp)
#define Iv_phase M_PI // Phase relative to PHASE_0

// Model and discretization
#define N_X 4
#define N_U 2

#define Ts 70.0e-6

// Constraints
#define INSERTION_INDEX_DEVIATION_ALLOWANCE 2.0

#define u1_min 0.0
#define u2_min 0.0

static const real_t U_MIN[N_U] = {u1_min, u2_min};

#define u1_max N_SM
#define u2_max N_SM

static const real_t U_MAX[N_U] = {u1_max, u2_max};

// Initial conditions
#define Iv_0 (Iv_ref_amp*sin(PHASE_0 + Iv_phase))
#define Icir_0 Icir_ref
#define Vsigma_u_0 Vdc
#define Vsigma_l_0 Vdc


static int simulate_mmc(const char* output_dir, size_t N, size_t simulation_timesteps) {
    // Allocation for output
    real_t *A = (real_t*)malloc(sizeof(real_t)*N*N_X*N_X);
    real_t *B = (real_t*)malloc(sizeof(real_t)*N*N_X*N_U);
    real_t *d = (real_t*)malloc(sizeof(real_t)*N*N_X);

    real_t *x = (real_t*)malloc(sizeof(real_t)*N*N_X);
    real_t *u = (real_t*)malloc(sizeof(real_t)*N*N_U);

    real_t *Vf = (real_t*)malloc(sizeof(real_t)*N);
    real_t *Iv_ref = (real_t*)malloc(sizeof(real_t)*N);

    real_t *xout = (real_t*)malloc(sizeof(real_t)*(simulation_timesteps+1)*N_X);
    real_t *uout = (real_t*)malloc(sizeof(real_t)*simulation_timesteps*N_U);
    real_t *tout = (real_t*)malloc(sizeof(real_t)*simulation_timesteps);

    // Initialize model
    mmc_model_get_fe_init(R, Rc, L, Lc, C, Ts, N_SM, N, CAST_3D_VLA(A, N_X, N_X), CAST_3D_VLA(B, N_X, N_U), CAST_2D_VLA(d, N_X));

    // Initialize solver
    sdqp_lmpc_mmc_init(N_X, N_U, N, q1, q2, N_SM, INSERTION_INDEX_DEVIATION_ALLOWANCE, U_MIN, U_MAX);

    // Initialize state
    xout[0] = Iv_0;
    xout[1] = Icir_0;
    xout[2] = Vsigma_u_0;
    xout[3] = Vsigma_l_0;

    // Initialize trajectory
    real_t omega = 2.0*M_PI*FREQ*Ts;
    for (size_t i = 0; i < N; ++i) {
        x[i*N_X + 0] = 0.0;
        x[i*N_X + 1] = Icir_ref;
        x[i*N_X + 2] = Vdc;
        x[i*N_X + 3] = Vdc;
    }
    for (size_t i = 0; i < N; ++i) {
        u[i*N_U + 0] = N_SM/2.0;
        u[i*N_U + 1] = N_SM/2.0;
    }

    long long simulation_time_ns = 0;
    for (size_t i = 0; i < simulation_timesteps; ++i) {
        timer_reset();
        // Predict trajectory
        memcpy(x, &xout[i*N_X], sizeof(real_t)*N_X); // We prefer simulation/measurement value to MPC value
        mmc_trajectory_shift(N_U, N, CAST_CONST_2D_VLA(u, N_U), CAST_2D_VLA(u, N_U));
        // Extrapolate references and disturbances
        extrapolate_sine(N, Vf_amp, FREQ, Ts, PHASE_0 + Vf_phase + 2.0*M_PI*FREQ*Ts*(real_t)i, 0.0, Vf);
        extrapolate_sine(N, Iv_ref_amp, FREQ, Ts, PHASE_0 + Iv_phase + 2.0*M_PI*FREQ*Ts*(real_t)i, 0.0, Iv_ref);
        // Get linearized discrete model
        mmc_model_get_fe(N, CAST_CONST_2D_VLA(x, N_X), CAST_CONST_2D_VLA(u, N_U), 
                Vf, Vdc, 
                CAST_3D_VLA(A, N_X, N_X), CAST_3D_VLA(B, N_X, N_U), CAST_2D_VLA(d, N_X));
        // Solve QP
        int err = sdqp_lmpc_mmc_solve(N_X, N_U, N, Iv_ref, Icir_ref, 
                CAST_CONST_3D_VLA(A, N_X, N_X), CAST_CONST_3D_VLA(B, N_X, N_U), CAST_CONST_2D_VLA(d, N_X), &xout[i*N_X], 
                CAST_2D_VLA(x, N_X), CAST_2D_VLA(u, N_U));
        if (err) {
            printf("Error while solving: %d\n", err);
            return 1;
        }
        long timestep_elapsed_ns = timer_elapsed_ns();
        tout[i] = ((real_t)timestep_elapsed_ns)/1000.0;
        simulation_time_ns += (long long)timestep_elapsed_ns;
        // Simulate using linearized discrete model
        simulate_affine(N_X, N_U, CAST_CONST_2D_VLA(A, N_X), &xout[i*N_X], CAST_CONST_2D_VLA(B, N_U), u, d, &xout[(i+1)*N_X]);
        memcpy(&uout[i*N_U], u, sizeof(real_t)*N_U);
    }

    // Timer output
    printf("%ld timesteps with horizon %ld finished in %lld ms\n", simulation_timesteps, N, simulation_time_ns/1000000);

    // Save output
    char path[128];
    sprintf(path, "%s/xoutN%ld.csv", output_dir, N);
    if (csv_save_matrix(path, simulation_timesteps+1, N_X, CAST_CONST_2D_VLA(xout, N_X)) < 0) {
        printf("Error while saving xout.\n");
        return 1;
    }
    sprintf(path, "%s/uoutN%ld.csv", output_dir, N);
    if (csv_save_matrix(path, simulation_timesteps, N_U, CAST_CONST_2D_VLA(uout, N_U)) < 0) {
        printf("Error while saving uout.\n");
        return 1;
    }
    sprintf(path, "%s/toutN%ld.csv", output_dir, N);
    if (csv_save_vector(path, simulation_timesteps, tout) < 0) {
        printf("Error while saving tout.\n");
        return 1;
    }

    //Cleanup
    sdqp_lmpc_mmc_cleanup();

    free(A);
    free(B);
    free(d);

    free(x);
    free(u);

    free(Vf);
    free(Iv_ref);

    free(xout);
    free(uout);
    free(tout);

    return 0;
}


int main(int argc, char *argv[]) {
    if (argc < 3) {
        return 1;
    } else {
        const char *output_dir = argv[1];
        size_t N = (size_t)atoi(argv[2]);
        size_t simulation_timesteps = (size_t)atoi(argv[3]);
        return simulate_mmc(output_dir, N, simulation_timesteps);
    }
    return 1;
}
