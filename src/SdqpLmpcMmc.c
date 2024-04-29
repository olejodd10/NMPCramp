#include "SdqpLmpcMmc.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "LinAlg.h"
#include "Ramp.h"
#include "Types.h"
#include "Utils.h"

static size_t m_n_x;
static size_t m_n_u;
static size_t m_N;

static size_t m_n_z;
static size_t m_n_M;
static size_t m_n_H;
static size_t m_n_a;

static real_t m_q1;
static real_t m_q2;

static const real_t *m_A;
static const real_t *m_B;

static const real_t *m_x_min;
static const real_t *m_x_max;
static real_t m_n_sm;
static real_t m_insertion_index_deviation_allowance;
static const real_t *m_u_min;
static const real_t *m_u_max;

static real_t *m_y;

static real_t *m_temp1; // Use x instead? Has length n_M and x is unused in initialize_y
static real_t *m_temp2;

// in-place version, v = (I-Ahat)^-1 * v
static void multiply_inv_eye_sub_Ahat_inplace(size_t n_x, size_t N, const real_t A[N][n_x][n_x], real_t v[N][n_x]) {
    for (size_t i = 1; i < N; ++i) {
        for (size_t j = 0; j < n_x; ++j) {
            v[i][j] += linalg_vector_inner_product(n_x, A[i][j], v[i-1]);
        }
    }
}

static void multiply_inv_eye_sub_Ahat_T_inplace(size_t n_x, size_t N, const real_t A[N][n_x][n_x], real_t v[N][n_x]) {
    for (size_t i = N-1; i >= 1; --i) { 
        for (size_t j = 0; j < n_x; ++j) {
            linalg_vector_add_scaled(n_x, v[i-1], A[i][j], v[i][j], v[i-1]);
        }
    }
}

static void get_column_M4a(
        size_t local_index, 

        size_t n_x,
        size_t n_u,
        size_t N,

        size_t n_M,
        size_t n_a,

        real_t q1,
        real_t q2,

        const real_t A[N][n_x][n_x],
        const real_t B[N][n_x][n_u],

        real_t *column_M4) 
{
    // Top block/first n_a elements first
    // Don't confuse MH2T P MH2 with MH2T P MH

    // Get column of Bhat
    memset(m_temp1, 0, sizeof(real_t)*n_M);
    for (size_t i = 0; i < n_x; ++i) {
        m_temp1[(local_index/n_u)*n_x+i] = B[local_index/n_u][i][local_index % n_u];
    }

    // Multiply with (...)^-1
    multiply_inv_eye_sub_Ahat_inplace(n_x, N, A, CAST_2D_VLA(m_temp1, n_x));

    // Multiply with Qhat
    for (size_t i = 0; i < N; ++i) {
        m_temp2[i*n_x + 0] = q1*m_temp1[i*n_x + 0];
        m_temp2[i*n_x + 1] = q2*m_temp1[i*n_x + 1];
        m_temp2[i*n_x + 2] = 0.0;
        m_temp2[i*n_x + 3] = 0.0;
    }

    // Multiply with (...)^-T
    multiply_inv_eye_sub_Ahat_T_inplace(n_x, N, A, CAST_2D_VLA(m_temp2, n_x));

    // Multiply with Bhat_T, write result to first n_a elements of column_M4
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < n_x; ++j) {
            linalg_vector_add_scaled(n_u, &column_M4[i*n_u], B[i][j], -m_temp2[i*n_x + j], &column_M4[i*n_u]); // Note the negation of the scaling factor
        }
    }

    // Add (subtract) Rhat column
    // R is zero for MMC

    // Add element on the diagonal, stemming from I
    column_M4[local_index] += 1.0;

    // Bottom block/last n_b elements
    // Don't confuse -HbMH2 with -HbMH
    // Note that m_temp1 already contains (...)^-1*Bhat(:,index), So we really just need to multiply with -Hb

    // Top two blocks
    for (size_t i = 0; i < n_M; ++i) {
        column_M4[n_a + i] = m_temp1[i]; // Note sign
        column_M4[n_a + n_M + i] = -m_temp1[i]; // Note sign
    }

    // M2 blocks. Remember the negation
    column_M4[n_a + 2*n_M + local_index/n_u] = 1.0; // Note the signs here!
    column_M4[n_a + 2*n_M + N + local_index/n_u] = -1.0;

    // The I-block of -HbMH2. Note the double negation
    column_M4[n_a + 2*(n_M+N) + local_index] += 1.0;
}

static void get_column_M4b(
        size_t local_index, 

        size_t n_x,
        size_t n_u,
        size_t N,

        size_t n_M,
        size_t n_a,

        const real_t A[N][n_x][n_x],
        const real_t B[N][n_x][n_u],

        real_t *column_M4) 
{
    // Top block/first n_a elements first
    if (local_index < 2*n_M) {
        memset(m_temp1, 0, sizeof(real_t)*n_M);
        if (local_index < n_M) {
            m_temp1[local_index] = -1.0;
        } else {
            m_temp1[local_index - n_M] = 1.0;
        } 
        multiply_inv_eye_sub_Ahat_T_inplace(n_x, N, A, CAST_2D_VLA(m_temp1, n_x));
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < n_x; ++j) {
                linalg_vector_add_scaled(n_u, &column_M4[i*n_u], B[i][j], m_temp1[i*n_x + j], &column_M4[i*n_u]);
            }
        }
    } else if (local_index < 2*n_M + N) {
        for (size_t i = 0; i < n_u; ++i) {
            column_M4[(local_index - 2*n_M)*n_u + i] = -1.0;
        }
    } else if (local_index < 2*(n_M + N)) {
        for (size_t i = 0; i < n_u; ++i) {
            column_M4[(local_index - 2*n_M - N)*n_u + i] = 1.0;
        }
    } else {
        column_M4[local_index - 2*(n_M + N)] = -1.0;
    }

    // Bottom block/last n_b elements
    // Add I-element
    column_M4[n_a + local_index] += 1.0;
}

static void get_column_M4(size_t index, real_t *column_M4) {
    memset(column_M4, 0, sizeof(real_t)*m_n_H);
    if (index < m_n_a) {
        get_column_M4a(index,
            m_n_x, m_n_u, m_N,
            m_n_M, m_n_a,
            m_q1, m_q2,
            CAST_CONST_3D_VLA(m_A, m_n_x, m_n_x), CAST_CONST_3D_VLA(m_B, m_n_x, m_n_u),
            column_M4);
    } else {
        get_column_M4b(index - m_n_a,
            m_n_x, m_n_u, m_N,
            m_n_M, m_n_a,
            CAST_CONST_3D_VLA(m_A, m_n_x, m_n_x), CAST_CONST_3D_VLA(m_B, m_n_x, m_n_u),
            column_M4);
    }
}

static void initialize_y(
        size_t n_x,
        size_t n_u,
        size_t N,

        size_t n_M,
        size_t n_H, 
        size_t n_a, 

        real_t q1,
        real_t q2,
        const real_t x1_ref[N],
        real_t x2_ref,

        const real_t A[N][n_x][n_x],
        const real_t B[N][n_x][n_u],
        const real_t d[N][n_x],

        const real_t x_min[n_x],
        const real_t x_max[n_x],
        real_t n_sm,
        real_t insertion_index_deviation_allowance,
        const real_t u_min[n_u],
        const real_t u_max[n_u],

        const real_t x0[n_x],

        real_t y[n_H]) 
{
    // MH[0 ha] (result n_z, but last n_a are just ha) (actually multiplication with MH2 (not transposed))
    for (size_t i = 0; i < N; ++i) {
        // Compute Bhat*[0 ha] 
        linalg_matrix_vector_product(n_x, n_u, B[i], u_max, &m_temp2[i*n_x]);
    }
    // Add d
    for (size_t i = 0; i < N; ++i) {
        linalg_vector_add(n_x, &m_temp2[i*n_x], d[i], &m_temp2[i*n_x]);
    }
    // Now, before multiplying with (...)^-1, is a perfect time to add A0x0 (result n_x) to get MH[A0x0+d ha] easily
    for (size_t i = 0; i < n_x; ++i) {
        m_temp2[i] += linalg_vector_inner_product(n_x, A[0][i], x0);
    }
    // Finish computing MH[A0x0+d ha]
    multiply_inv_eye_sub_Ahat_inplace(n_x, N, A, CAST_2D_VLA(m_temp2, n_x));

    // Multiply out of place with Hb, store result in last n_b of y
    for (size_t i = 0; i < n_M; ++i) {
        y[n_a + i] = -m_temp2[i];
    }
    for (size_t i = 0; i < n_M; ++i) {
        y[n_a + n_M + i] = m_temp2[i];
    }
    real_t temp = 0.0;
    for (size_t i = 0; i < n_u; ++i) {
        temp += u_max[i];
    }
    for (size_t i = 0; i < N; ++i) {
        y[n_a + 2*n_M + i] = -temp;
        y[n_a + 2*n_M + N + i] = temp;
    }
    for (size_t i = 0; i < n_a; ++i) {
        y[n_a + 2*(n_M + N) + i] = -u_max[i % n_u];
    }

    // Subtract hb from last n_b of y
    // Note the signs
    for (size_t i = 0; i < n_M; ++i) {
        y[n_a + i] += x_min[i % n_x];
    }
    for (size_t i = 0; i < n_M; ++i) {
        y[n_a + n_M + i] -= x_max[i % n_x];
    }
    for (size_t i = 0; i < N; ++i) {
        y[n_a + 2*n_M + i] -= insertion_index_deviation_allowance - n_sm;
    }
    for (size_t i = 0; i < N; ++i) {
        y[n_a + 2*n_M + N + i] -= insertion_index_deviation_allowance + n_sm;
    }
    for (size_t i = 0; i < n_u*N; ++i) {
        y[n_a + 2*(n_M + N) + i] += u_min[i % n_u];
    }

    // Multiply MH[A0x0+d ha] with P, then add f (result n_z, but n_a of it can go directly to y because of I in MH2T)
    for (size_t i = 0; i < N; ++i) {
        m_temp1[i*n_x + 0] = q1*(m_temp2[i*n_x + 0] - x1_ref[i]);
        m_temp1[i*n_x + 1] = q2*(m_temp2[i*n_x + 1] - x2_ref);
        m_temp1[i*n_x + 2] = 0.0;
        m_temp1[i*n_x + 3] = 0.0;
    }
    // R and f_u are zero for MMC, so latter parts vanish

    // Multiply MH2T (result n_a, stored in first n_a of y)
    multiply_inv_eye_sub_Ahat_T_inplace(n_x, N, A, CAST_2D_VLA(m_temp1, n_x));
    memset(y, 0, sizeof(real_t)*n_a); // Needed for transpose product
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < n_x; ++j) {
            linalg_vector_add_scaled(n_u, &y[i*n_u], B[i][j], m_temp1[i*n_x + j], &y[i*n_u]);
        }
    }
    
    // "Optional":
    // If starting from non-empty active set and invq is correct, compute y = invq*(M3x0+m5)
}

static void compute_x_u(size_t n_x, size_t n_u, size_t N, size_t n_H, 
        const real_t A[N][n_x][n_x], const real_t B[N][n_x][n_u], const real_t d[N][n_x], const real_t x0[n_x], const real_t u_max[n_u], const real_t y[n_H],
        real_t x[N][n_x], real_t u[N][n_u])
{
    // Find ha-ra, store in last n_a elements of z (u)
    for (size_t i = 0; i < N*n_u; ++i) { // n_a == N*n_u
        u[i/n_u][i%n_u] = y[i] > 0.0 ? u_max[i % n_u] - y[i] : u_max[i % n_u]; // y_a is in the start of y
    }

    // Multiply ha-ra with Bhat, and write result to first n_M elements of z (x)
    for (size_t i = 0; i < N; ++i) {
        linalg_matrix_vector_product(n_x, n_u, B[i], u[i], x[i]);
    }
    // Add d
    for (size_t i = 0; i < N; ++i) {
        linalg_vector_add(n_x, x[i], d[i], x[i]);
    }
    // Add A0x0 to first n_x elements
    for (size_t i = 0; i < n_x; ++i) {
        x[0][i] += linalg_vector_inner_product(n_x, A[0][i], x0);
    }
    
    // Multiply x with (I-Ahat)^-1 and write result to x
    multiply_inv_eye_sub_Ahat_inplace(n_x, N, A, x);
}

void sdqp_lmpc_mmc_init(
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
) {
	m_n_x = n_x;
	m_n_u = n_u;
	m_N = N;

	m_n_z = (n_x+n_u)*N;
    m_n_M = n_x*N;
    m_n_H = 2*N*(n_x+1+n_u);
    m_n_a = n_u*N;

    m_q1 = q1;
    m_q2 = q2;

    m_x_min = (real_t*)x_min;
    m_x_max = (real_t*)x_max;
    m_n_sm = n_sm;
    m_insertion_index_deviation_allowance = insertion_index_deviation_allowance;
    m_u_min = (real_t*)u_min;
    m_u_max = (real_t*)u_max;

    m_y = (real_t*)malloc(sizeof(real_t)*m_n_H);

    m_temp1 = (real_t*)malloc(sizeof(real_t)*m_n_M);
    m_temp2 = (real_t*)malloc(sizeof(real_t)*m_n_M);

    ramp_init(m_n_H, m_n_a, get_column_M4);
    ramp_enable_infeasibility_error(1e-12, 1e12);
}

void sdqp_lmpc_mmc_cleanup(void) {
	free(m_y);

	free(m_temp1);
	free(m_temp2);

    ramp_cleanup();
}

int sdqp_lmpc_mmc_solve(size_t n_x, size_t n_u, size_t N, const real_t x1_ref[N], real_t x2_ref, const real_t A[N][n_x][n_x], const real_t B[N][n_x][n_u], const real_t d[N][n_x], const real_t x0[n_x], real_t x[N][n_x], real_t u[N][n_u]) {
    // initialize y
    m_A = (real_t*)A;
    m_B = (real_t*)B;
    initialize_y(n_x, n_u, N, m_n_M, m_n_H, m_n_a,
        m_q1, m_q2, x1_ref, x2_ref,
		A, B, d,
        m_x_min, m_x_max, m_n_sm, m_insertion_index_deviation_allowance, m_u_min, m_u_max, 
        x0, m_y);
    int err = ramp_solve(m_n_H, m_n_a, RAMP_HOTSTART_NONE, m_y);
    if (err) {
        return err;
    }
    compute_x_u(n_x, n_u, N, m_n_H, 
            A, B, d,
            x0, m_u_max, m_y, 
            x, u);
    return 0;
}
