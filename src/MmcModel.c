#include "MmcModel.h"

#include <stddef.h>
#include <string.h>

#include "LinAlg.h"
#include "Types.h"

// (2.14) linearized and discretized using forward Euler.
// The matrices below match (2.14) except that they are multiplied with Ts.

static real_t A_00 = 0.0;
static real_t A_11 = 0.0;

static real_t B0_02 = 0.0;
static real_t B0_12 = 0.0;
static real_t B0_20 = 0.0;
static real_t B0_21 = 0.0;

static real_t B1_03 = 0.0;
static real_t B1_13 = 0.0;
static real_t B1_30 = 0.0;
static real_t B1_31 = 0.0;

static real_t d_0 = 0.0;
static real_t d_1 = 0.0;

void mmc_model_get_init(real_t R, real_t Rc, real_t L, real_t Lc, real_t C, real_t Ts, real_t n_sm, size_t N, real_t A[N][N_X][N_X], real_t B[N][N_X][N_U], real_t d[N][N_X]) {
    A_00 = -Ts*(R+2.0*Rc)/(L+2.0*Lc);
    A_11 = -Ts*R/L;

    B0_02 = Ts/((L+2.0*Lc)*n_sm);
    B0_12 = -Ts/(2*n_sm*L);
    B0_20 = -Ts/(2.0*C);
    B0_21 = Ts/C;

    B1_03 = -B0_02;
    B1_13 = B0_12;
    B1_30 = -B0_20;
    B1_31 = B0_21;

    d_0 = 2.0*Ts/(L+2.0*Lc);
    d_1 = Ts/(2.0*L);

    // Sets matrices to zero in unused elements, and initializes elements that rely on parameters but not trajectories
    memset(A, 0, sizeof(real_t)*N*N_X*N_X);
    memset(B, 0, sizeof(real_t)*N*N_X*N_U);
    memset(d, 0, sizeof(real_t)*N*N_X);
    
    for (size_t i = 0; i < N; ++i) {
        A[i][0][0] = 1.0 + A_00;
        A[i][1][1] = 1.0 + A_11;
        A[i][2][2] = 1.0;
        A[i][3][3] = 1.0;
    }
}

// Assumes A, B and d are zero in unused elements
void mmc_model_get(size_t N, const real_t x[N][N_X], const real_t u[N][N_U], const real_t vf[N], real_t Vdc, real_t A[N][N_X][N_X], real_t B[N][N_X][N_U], real_t d[N][N_X]) {
    for (size_t i = 0; i < N; ++i) {
        A[i][0][2] = B0_02*u[i][0];
        A[i][0][3] = B1_03*u[i][1];
        A[i][1][2] = B0_12*u[i][0];
        A[i][1][3] = B1_13*u[i][1];
        A[i][2][0] = B0_20*u[i][0];
        A[i][2][1] = B0_21*u[i][0];
        A[i][3][0] = B1_30*u[i][1];
        A[i][3][1] = B1_31*u[i][1];

        B[i][0][0] = B0_02*x[i][2];
        B[i][0][1] = B1_03*x[i][3];
        B[i][1][0] = B0_12*x[i][2];
        B[i][1][1] = B1_13*x[i][3];
        B[i][2][0] = B0_20*x[i][0] + B0_21*x[i][1];
        B[i][3][1] = B1_30*x[i][0] + B1_31*x[i][1];

        d[i][0] = d_0*vf[i] - linalg_vector_inner_product(N_U, B[i][0], u[i]);
        d[i][1] = d_1*Vdc - linalg_vector_inner_product(N_U, B[i][1], u[i]);    
        d[i][2] = - linalg_vector_inner_product(N_U, B[i][2], u[i]);
        d[i][3] = - linalg_vector_inner_product(N_U, B[i][3], u[i]);
    }
}
