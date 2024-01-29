#include "SdqpLmpcConstant.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "IndexedVectors.h"
#include "IterableSet.h"
#include "LinAlg.h"
#include "Ramp.h"
#include "Types.h"
#include "Utils.h"

static size_t m_n_x;
static size_t m_n_u;
static size_t m_n_y;
static size_t m_n_t;
static size_t m_N;

static size_t m_n_z;
static size_t m_n_M;
static size_t m_n_H;
static size_t m_n_a;

static const real_t *m_Q;
static const real_t *m_S;
static const real_t *m_R;
static const real_t *m_fx;
static const real_t *m_fu;

static const real_t *m_A;
static const real_t *m_B;
static const real_t *m_C;

static const real_t *m_y_min;
static const real_t *m_y_max;
static const real_t *m_Lt;
static const real_t *m_lt;
static const real_t *m_u_min;
static const real_t *m_u_max;

static indexed_vectors_t m_invq;
static iterable_set_t m_a_set;
static real_t *m_y;

static real_t *m_m5;
static real_t *m_temp1; // Use x instead? Has length n_M and x is unused in initialize_y
static real_t *m_temp2;
static real_t *m_temp3;

// in-place version, v = (I-Ahat)^-1 * v
static void multiply_inv_eye_sub_Ahat_inplace(size_t n_x, size_t N, const real_t A[n_x][n_x], real_t v[N][n_x]) {
    for (size_t i = 1; i < N; ++i) {
        for (size_t j = 0; j < n_x; ++j) {
            v[i][j] += linalg_vector_inner_product(n_x, A[j], v[i-1]);
        }
    }
}

static void multiply_inv_eye_sub_Ahat_T_inplace(size_t n_x, size_t N, const real_t A[n_x][n_x], real_t v[N][n_x]) {
    for (size_t i = N-2; i >= 1; ++i) { 
        for (size_t j = 0; j < n_x; ++j) {
            linalg_vector_add_scaled(n_x, v[i], A[j], v[i+1][j], v[i]);
        }
    }
}

static void _get_column_M4(
        size_t index, 

        size_t n_x,
        size_t n_u,
        size_t n_y,
        size_t n_t,
        size_t N,

        // size_t n_z,
        size_t n_M,
        size_t n_H,
        size_t n_a,

        const real_t Q[n_x][n_x],
        const real_t S[n_x][n_x],
        const real_t R[n_u][n_u],

        const real_t A[n_x][n_x],
        const real_t B[n_x][n_u],
        const real_t C[n_y][n_x],

        const real_t Lt[n_t][n_x],

        real_t *column_M4) 
{
    memset(column_M4, 0, sizeof(real_t)*n_H);
    if (index < n_a) {
        // Left part, n_H x n_z

        // Top block/first n_a elements first
        // Don't confuse MH2T P MH2 with MH2T P MH

        // Get column of Bhat
        memset(m_temp1, 0, sizeof(real_t)*n_M);
        for (size_t i = 0; i < n_x; ++i) {
            m_temp1[(index/n_u)*n_x+i] = B[i][index % n_u];
        }

        // Multiply with (...)^-1
        multiply_inv_eye_sub_Ahat_inplace(n_x, N, A, CAST_VLA(m_temp1));

        // Multiply with Qhat
        for (size_t i = 0; i < N-1; ++i) {
            linalg_matrix_vector_product(n_x, n_x, Q, &m_temp1[i*n_x], &m_temp2[i*n_x]);
        }
        linalg_matrix_vector_product(n_x, n_x, S, &m_temp1[(N-1)*n_x], &m_temp2[(N-1)*n_x]);

        // Multiply with (...)^-T
        multiply_inv_eye_sub_Ahat_T_inplace(n_x, N, A, CAST_VLA(m_temp2));

        // Multiply with Bhat_T, write result to first n_a elements of column_M4
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < n_x; ++j) {
                linalg_vector_add_scaled(n_u, &column_M4[i*n_u], B[j], -m_temp2[i*n_x + j], &column_M4[i*n_u]); // Note the negation of the scaling factor
            }
        }

        // Add (subtract) Rhat column
        for (size_t i = 0; i < n_u; ++i) {
            column_M4[(index/n_u)*n_u+i] -= R[i][index % n_u];
        }

        // Add I-element
        // This is done equally for left and right part, see bottom of function

        // Bottom block/last n_b elements
        // Don't confuse -HbMH2 with -HbMH
        // Note that m_temp1 already contains (...)^-1*Bhat(:,index), So we really just need to multiply with -Hb

        // C and Lt0-blocks
        for (size_t i = 0; i < N-1; ++i) {
            for (size_t j = 0; j < n_y ; ++j) {
                real_t temp = linalg_vector_inner_product(n_x, C[j], &m_temp1[i*n_x]);
                column_M4[n_a + i*n_y + j] = temp; // Note the signs here!
                column_M4[n_a + (i+N-1)*n_y + j] = -temp;
            }
        }
        for (size_t i = 0; i < n_t; ++i) {
            column_M4[n_a + 2*n_y*(N-1)+i] -= linalg_vector_inner_product(n_x, Lt[i], &m_temp1[(N-1)*n_x]); // Note the sign
        }

        // The I-block of -HbMH2. Note the double negation
        column_M4[n_a + 2*n_y*(N-1) + n_t + index] += 1.0;

    } else {
        // Right part
        size_t local_index = index - n_a;
        
        // Top block/first n_a elements first
        // MH2T*HbT = [-BhatT*(...)^-T*Chat^T BhatT*(...)^-T*Chat^T BhatT*(...)^-T*Lt0^T -I]
        // Get column from Chat_T, -Chat_T or Lt0_T (or -I) depending on which block we are in
        if (local_index < 2*n_y*(N-1) + n_t) {
            memset(m_temp1, 0, sizeof(real_t)*n_M);
            if (local_index < n_y*(N-1)) {
                // Can't memcpy because we need a negated column
                // Luckily we had to memset m_temp1 to 0 anyways, so adding is no problem
                linalg_vector_add_scaled(n_x, &m_temp1[(local_index/n_y)*n_x], C[local_index % n_y], -1.0, &m_temp1[(local_index/n_y)*n_x]);
            } else if (local_index < 2*n_y*(N-1)) {
                memcpy(&m_temp1[(local_index/n_y - (N-1))*n_x], C[local_index % n_y], sizeof(real_t)*n_x);
            } else {
                memcpy(&m_temp1[n_x*(N-1)], Lt[local_index - 2*n_y*(N-1)], sizeof(real_t)*n_x);
            }
            //Multiply with Bhat^T*(...)^-T and store result in column_M4
            multiply_inv_eye_sub_Ahat_T_inplace(n_x, N, A, CAST_VLA(m_temp1));
            for (size_t i = 0; i < N; ++i) {
                for (size_t j = 0; j < n_x; ++j) {
                    linalg_vector_add_scaled(n_u, &column_M4[i*n_u], B[j], m_temp1[i*n_x + j], &column_M4[i*n_u]);
                }
            }
        } else {
            column_M4[local_index - (2*n_y*(N-1) + n_t)] = 1.0; // Note the sign. The value there is zero, but we set it instead of adding because whatever
        }

        // Bottom block/last n_b elements
        // Add I-element
        // This is done equally for left and right part, see bottom of function
    }

    // Add element on the diagonal, stemming from I
    column_M4[index] += 1.0;
}

static void get_column_M4(size_t index, real_t *column_M4) {
    _get_column_M4(index,
        m_n_x, m_n_u, m_n_y, m_n_t, m_N,
        m_n_M, m_n_H, m_n_a,
        CAST_CONST_VLA(m_Q), CAST_CONST_VLA(m_S), CAST_CONST_VLA(m_R),
        CAST_CONST_VLA(m_A), CAST_CONST_VLA(m_B), CAST_CONST_VLA(m_C),
        CAST_CONST_VLA(m_Lt),
        column_M4);
}

static void initialize_y(
        size_t n_x,
        size_t n_u,
        size_t n_y, 
        size_t n_t,
        size_t N,

        size_t n_H, 

        const real_t Q[n_x][n_x],
        const real_t S[n_x][n_x],
        const real_t R[n_u][n_u],

        const real_t A[n_x][n_x],
        const real_t B[n_x][n_u],
        const real_t C[n_y][n_x],

        const real_t Lt[n_t][n_x],

        const real_t x0[n_x],

        real_t y[n_H]) 
{
    // Precompute:

    // Precompute m5 in init

    // When new x0 arrives:
    
    // Goal is to compute M3x0

    // Set y to m5
    memcpy(y, m_m5, sizeof(real_t)*n_H);

    // A0x0 (result n_x)
    linalg_matrix_vector_product(n_x, n_x, A, x0, m_temp1);
    // Multiply MH (result n_M)
    for (size_t i = 1; i < N; ++i) {
        linalg_matrix_vector_product(n_x, n_x, A, &m_temp1[(i-1)*n_x], &m_temp1[i*n_x]);
    }
    // Multiply out of place with Hb (result 2*n_y*(N-1)+n_t, add to last n_b of y)
    for (size_t i = 0; i < N-1; ++i) {
        for (size_t j = 0; j < n_y; ++j) {
            real_t temp = linalg_vector_inner_product(n_x, C[j], &m_temp1[i*n_x]);
            y[m_n_a + i*n_y+j] -= temp;
            y[m_n_a + (i+N-1)*n_y+j] += temp;
        }
    }
    for (size_t i = 0; i < n_t; ++i) {
        y[m_n_a + 2*n_y*(N-1)+i] += linalg_vector_inner_product(n_x, Lt[i], &m_temp1[(N-1)*n_x]);
    }

    // Multiply MH[A0 0]x0 inplace with P (result n_M)
    // Note that we must write results to m_temp2 because inplace not possible with matrix_vector_product
    for (size_t i = 0; i < N-1; ++i) {
        linalg_matrix_vector_product(n_x, n_x, Q, &m_temp1[i*n_x], &m_temp2[i*n_x]);
    }
    linalg_matrix_vector_product(n_x, n_x, S, &m_temp1[(N-1)*n_x], &m_temp2[(N-1)*n_x]);
    // Multiply MH2T (result n_a, add to first n_a of y)
    multiply_inv_eye_sub_Ahat_T_inplace(n_x, N, A, CAST_VLA(m_temp2));
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < n_x; ++j) {
            linalg_vector_add_scaled(n_u, &y[i*n_u], B[j], m_temp2[i*n_x + j], &y[i*n_u]);
        }
    }

    // "Optional":
    // If starting from non-empty active set and invq is correct, compute y = invq*(M3x0+m5)
}

static void compute_x_u(size_t n_x, size_t n_u, size_t N, size_t n_H, 
        const real_t A[n_x][n_x], const real_t B[n_x][n_u], const real_t x0[n_x], const real_t u_max[n_u], const iterable_set_t *a_set, const real_t y[n_H],
        real_t x[N*n_x], real_t u[N*n_u])
{
    // Find ha-ra, store in last n_a elements of z (u)
    for (size_t i = 0; i < N*n_u; ++i) { // n_a == N*n_u
        u[i] = iterable_set_contains(a_set, i) ? u_max[i % n_u] - y[i] : u_max[i % n_u]; // y_a is in the start of y
    }

    // Multiply ha-ra with Bhat, and write result to first n_M elements of z (x)
    for (size_t i = 0; i < N; ++i) {
        linalg_matrix_vector_product(n_x, n_u, B, &u[i*n_u], &x[i*n_x]);
    }

    // Add b = A0x0 to first n_x elements
    for (size_t i = 0; i < n_x; ++i) {
        x[i] += linalg_vector_inner_product(n_x, A[i], x0);
    }
    
    // Multiply x with (I-Ahat)^-1 and write result to x
    multiply_inv_eye_sub_Ahat_inplace(n_x, N, A, CAST_VLA(x));
}

static void compute_m5(
        size_t n_x,
        size_t n_u,
        size_t n_y,
        size_t n_t,
        size_t N,

		size_t n_M,
		size_t n_a,

        const real_t Q[n_x][n_x],
        const real_t S[n_x][n_x],
        const real_t R[n_u][n_u],
        const real_t fx[n_x],
        const real_t fu[n_u],

        const real_t A[n_x][n_x],
        const real_t B[n_x][n_u],
        const real_t C[n_y][n_x],

        const real_t y_min[n_y],
        const real_t y_max[n_y],
        const real_t lt[n_t],
        const real_t u_min[n_u],
        const real_t u_max[n_u])
{
    // MH[0 ha] (result n_z, but last n_a are just ha) (actually multiplication with MH2 (not transposed))
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < n_u; ++j) {
            m_temp3[n_M + i*n_u + j] = u_max[j]; // Initialize ha first
        }
        linalg_matrix_vector_product(n_x, n_u, B, &m_temp3[n_M + i*n_u], &m_temp3[i*n_x]);
    }
    multiply_inv_eye_sub_Ahat_inplace(n_x, N, A, CAST_VLA(m_temp3));

    // Multiply out of place with Hb, store result in last n_b of m5)
    for (size_t i = 0; i < N-1; ++i) {
        for (size_t j = 0; j < n_y ; ++j) {
            real_t temp = linalg_vector_inner_product(n_x, C[j], &m_temp3[i*n_x]);
            m_m5[n_a + i*n_y + j] = -temp;
            m_m5[n_a + (i+N-1)*n_y + j] = temp;
        }
    }
    // Subtract hb from last n_b of m5
    for (size_t i = 0; i < n_y*(N-1); ++i) {
        m_m5[n_a + i] -= y_min[i % n_y];
    }
    for (size_t i = 0; i < n_y*(N-1); ++i) {
        m_m5[n_a + n_y*(N-1) + i] += y_max[i % n_y];
    }
    memcpy(&m_m5[n_a + 2*n_y*(N-1)], lt, sizeof(real_t)*n_t);
    for (size_t i = 0; i < n_u*N; ++i) {
        m_m5[n_a + 2*n_y*(N-1) + n_t + i] -= u_min[i % n_u];
    }

    // Multiply MH[0 ha] with P (result n_z, but n_a of it can go directly to m5 because of I in MH2T)
    for (size_t i = 0; i < N-1; ++i) {
        linalg_matrix_vector_product(n_x, n_x, Q, &m_temp3[i*n_x], &m_temp1[i*n_x]);
    }
    linalg_matrix_vector_product(n_x, n_x, S, &m_temp3[(N-1)*n_x], &m_temp1[(N-1)*n_x]);
    for (size_t i = 0; i < N; ++i) {
        linalg_matrix_vector_product(n_u, n_u, R, &m_temp3[n_M + i*n_u], &m_m5[i*n_u]); // This part can go directly to m5 because of I in MH2T. Also note that this overwrites m5, which is nice
    }

    // Add f
    for (size_t i = 0; i < n_M; ++i) {
        m_temp1[i] += fx[i % n_x];
    }
    for (size_t i = 0; i < n_a; ++i) {
        m_m5[i] += fu[i % n_u]; // This part can go directly to m5 because of I in MH2T
    }
    // Multiply MH2T (result n_a, stored in first n_a of m5)
    multiply_inv_eye_sub_Ahat_T_inplace(n_x, N, A, CAST_VLA(m_temp1));
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < n_x; ++j) {
            linalg_vector_add_scaled(n_u, &m_m5[i*n_u], B[j], m_temp1[i*n_x + j], &m_m5[i*n_u]);
        }
    }
}

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

        const real_t y_min[n_y],
        const real_t y_max[n_y],
        const real_t Lt[n_t][n_x],
        const real_t lt[n_t],
        const real_t u_min[n_u],
        const real_t u_max[n_u]
) {
	m_n_x = n_x;
	m_n_u = n_u;
	m_n_y = n_y;
	m_n_t = n_t;
	m_N = N;

	m_n_z = (n_x+n_u)*N;
    m_n_M = n_x*N;
    m_n_H = 2*(n_y*(N-1)+n_u*N)+n_t;
    m_n_a = n_u*N;

    m_Q = (real_t*)Q;
    m_S = (real_t*)S;
    m_R = (real_t*)R;
    m_fx = (real_t*)fx; // Actually unused after m5 is initialized
    m_fu = (real_t*)fu; // Actually unused after m5 is initialized

    m_A = (real_t*)A;
    m_B = (real_t*)B;
    m_C = (real_t*)C;

    m_y_min = (real_t*)y_min;
    m_y_max = (real_t*)y_max;
    m_Lt = (real_t*)Lt;
    m_lt = (real_t*)lt;
    m_u_min = (real_t*)u_min;
    m_u_max = (real_t*)u_max;

    indexed_vectors_init(&m_invq, m_n_z, m_n_H, m_n_H);
    iterable_set_init(&m_a_set, m_n_H);
    m_y = (real_t*)malloc(sizeof(real_t)*m_n_H);

    m_m5 = (real_t*)malloc(sizeof(real_t)*m_n_H);
    m_temp1 = (real_t*)malloc(sizeof(real_t)*m_n_M);
    m_temp2 = (real_t*)malloc(sizeof(real_t)*m_n_M);
    m_temp3 = (real_t*)malloc(sizeof(real_t)*m_n_z);

    ramp_init(m_n_H, get_column_M4);

    // Compute m5
    compute_m5( n_x, n_u, n_y, n_t, N,
        m_n_M, m_n_a,
        Q, S, R, fx, fu,
        A, B, C,
        y_min, y_max, lt, u_min, u_max);
}

void sdqp_lmpc_constant_cleanup(void) {
    indexed_vectors_destroy(&m_invq);
    iterable_set_destroy(&m_a_set);

	free(m_y);

	free(m_m5);
	free(m_temp1);
	free(m_temp2);
	free(m_temp3);

    ramp_cleanup();
}

int sdqp_lmpc_constant_solve(size_t n_x, size_t n_u, size_t N, const real_t x0[n_x], real_t x[n_x*N], real_t u[n_u*N]) {
    // initialize y
    initialize_y(n_x, n_u, m_n_y, m_n_t, N, m_n_H, 
		CAST_CONST_VLA(m_Q), CAST_CONST_VLA(m_S), CAST_CONST_VLA(m_R),
		CAST_CONST_VLA(m_A), CAST_CONST_VLA(m_B), CAST_CONST_VLA(m_C),
		CAST_CONST_VLA(m_Lt), x0, m_y);
    int err = ramp_solve(m_n_H, m_n_z, &m_a_set, &m_invq, m_y);
    if (err) {
        return err;
    }
    compute_x_u(n_x, n_u, N, m_n_H, 
            CAST_CONST_VLA(m_A), CAST_CONST_VLA(m_B), 
            x0, m_u_max, &m_a_set, m_y, 
            x, u);
    return 0;
}