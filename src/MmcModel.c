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

static real_t m_Ts0 = 0.0;
static real_t m_Ts1 = 0.0;

static real_t m_Al_T[N_X][N_X];
static real_t m_Bl_T[N_U][N_X];
static real_t m_dl[N_X];

static real_t m_temp1[N_X][N_X];
static real_t m_temp2[N_X][N_X];
static real_t m_temp3[N_X][N_X];

void mmc_model_get_fe_init(real_t R, real_t Rc, real_t L, real_t Lc, real_t C, real_t Ts, real_t n_sm, size_t N, real_t A[N][N_X][N_X], real_t B[N][N_X][N_U], real_t d[N][N_X]) {
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
void mmc_model_get_fe(size_t N, const real_t x[N][N_X], const real_t u[N][N_U], const real_t vf[N], real_t Vdc, real_t A[N][N_X][N_X], real_t B[N][N_X][N_U], real_t d[N][N_X]) {
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

void mmc_model_get_rk4_init(real_t R, real_t Rc, real_t L, real_t Lc, real_t C, real_t Ts0, real_t Ts1, real_t n_sm, size_t N) {
    m_Ts0 = Ts0;
    m_Ts1 = Ts1;

    A_00 = -1.0*(R+2.0*Rc)/(L+2.0*Lc);
    A_11 = -1.0*R/L;

    B0_02 = 1.0/((L+2.0*Lc)*n_sm);
    B0_12 = -1.0/(2*n_sm*L);
    B0_20 = -1.0/(2.0*C);
    B0_21 = 1.0/C;

    B1_03 = -B0_02;
    B1_13 = B0_12;
    B1_30 = -B0_20;
    B1_31 = B0_21;

    d_0 = 2.0/(L+2.0*Lc);
    d_1 = 1.0/(2.0*L);

    memset(m_Al_T, 0, sizeof(real_t)*N_X*N_X);
    memset(m_Bl_T, 0, sizeof(real_t)*N_X*N_U);
    memset(m_dl, 0, sizeof(real_t)*N_X);

    m_Al_T[0][0] = A_00;
    m_Al_T[1][1] = A_11;
}


static void linearize(size_t N, const real_t x[N_X], const real_t u[N_U], real_t Vf, real_t Vdc, real_t Al_T[N_X][N_X], real_t Bl_T[N_U][N_X], real_t dl[N_X]) {
        Al_T[2][0] = B0_02*u[0];
        Al_T[3][0] = B1_03*u[1];
        Al_T[2][1] = B0_12*u[0];
        Al_T[3][1] = B1_13*u[1];
        Al_T[0][2] = B0_20*u[0];
        Al_T[1][2] = B0_21*u[0];
        Al_T[0][3] = B1_30*u[1];
        Al_T[1][3] = B1_31*u[1];

        Bl_T[0][0] = B0_02*x[2];
        Bl_T[1][0] = B1_03*x[3];
        Bl_T[0][1] = B0_12*x[2];
        Bl_T[1][1] = B1_13*x[3];
        Bl_T[0][2] = B0_20*x[0] + B0_21*x[1];
        Bl_T[1][3] = B1_30*x[0] + B1_31*x[1];

        dl[0] = d_0*Vf;
        dl[1] = d_1*Vdc;
        dl[2] = 0.0;
        dl[3] = 0.0;
        for (size_t i = 0; i < N_U; ++i) {
            linalg_vector_add_scaled(N_X, dl, Bl_T[i], -u[i], dl);
        }
}

static void discretize_rk4(const real_t x[N_X], const real_t u[N_U], real_t Ts, const real_t Al_T[N_X][N_X], const real_t Bl_T[N_U][N_X], const real_t dl[N_X], real_t A[N_X][N_X], real_t B[N_X][N_U], real_t d[N_X]) {
    memset(m_temp1, 0, sizeof(real_t)*N_X*N_X);
    memset(m_temp3, 0, sizeof(real_t)*N_X*N_X);
    for (size_t i = 0; i < N_X; ++i) {
        m_temp1[i][i] = Ts;
        m_temp3[i][i] = Ts;
    }

    for (size_t k = 1; k < 4; ++k) { // The upper limit is the order of the RK method
        if (k % 2 == 1) {
            for (size_t i = 0; i < N_X; ++i) {
                for (size_t j = 0; j < N_X; ++j) {
                    m_temp2[i][j] = (Ts/(real_t)(k+1))*linalg_vector_inner_product(N_X, m_temp1[i], Al_T[j]);
                }
            }
            for (size_t i = 0; i < N_X; ++i) {
                linalg_vector_add(N_X, m_temp3[i], m_temp2[i], m_temp3[i]);
            }
        } else {
            for (size_t i = 0; i < N_X; ++i) {
                for (size_t j = 0; j < N_X; ++j) {
                    m_temp1[i][j] = (Ts/(real_t)(k+1))*linalg_vector_inner_product(N_X, m_temp2[i], Al_T[j]);
                }
            }
            for (size_t i = 0; i < N_X; ++i) {
                linalg_vector_add(N_X, m_temp3[i], m_temp1[i], m_temp3[i]);
            }
        }
    }

    linalg_matrix_product(N_X, N_X, N_X, m_temp3, Al_T, A);
    for (size_t i = 0; i < N_X; ++i) {
        A[i][i] += 1.0;
    }

    linalg_matrix_product(N_X, N_X, N_U, m_temp3, Bl_T, B);

    linalg_matrix_vector_product(N_X, N_X, m_temp3, dl, d);
}

void mmc_model_get_rk4(size_t N, const real_t x[N][N_X], const real_t u[N][N_U], const real_t vf[N], real_t Vdc, real_t A[N][N_X][N_X], real_t B[N][N_X][N_U], real_t d[N][N_X]) {
    for (size_t i = 0; i < N; ++i) {
        real_t Ts = i == 0 ? m_Ts0 : m_Ts1;
        linearize(N, x[i], u[i], vf[i], Vdc, m_Al_T, m_Bl_T, m_dl);
        discretize_rk4(x[i], u[i], Ts, m_Al_T, m_Bl_T, m_dl, A[i], B[i], d[i]);
    }
}
