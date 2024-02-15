#include "matrix.h"
#include "mex.h"

#include <stddef.h>
#include <math.h>

#include "MexUtils.h"
#include "Utils.h"
#include "SdqpLmpcMmc.h"
#include "Types.h"

#define NRHS 17
#define NLHS 2

static int initialized = 0;

static void cleanup(void) {
    sdqp_lmpc_mmc_cleanup();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs!=NRHS) {
        mexErrMsgTxt("Incorrect number of inputs.");
    }
    if(nlhs!=NLHS) {
        mexErrMsgTxt("Incorrect number of outputs.");
    }

    for (int i = 0; i < NRHS; ++i) {
        if( !mxIsSingle(prhs[i]) || 
             mxIsComplex(prhs[i])) {
            mexErrMsgTxt("Input must be type single.");
        }
    }

    // Note that transposed matrices are needed because of MATLABs column-major memory layout.
    size_t n_x = (size_t)round(mxGetScalar(prhs[0]));
    size_t n_u = (size_t)round(mxGetScalar(prhs[1]));
    size_t N = (size_t)round(mxGetScalar(prhs[2]));
    real_t q1 = mxGetScalar(prhs[3]);
    real_t q2 = mxGetScalar(prhs[4]);

    if (mex_assert_1d_array_dimensions(prhs[5], n_x))
        mexErrMsgTxt("Input x_min must have n_x elements.");
    real_t *x_min = mxGetSingles(prhs[5]);
    if (mex_assert_1d_array_dimensions(prhs[6], n_x))
        mexErrMsgTxt("Input x_max must have n_x elements.");
    real_t *x_max = mxGetSingles(prhs[6]);

    real_t n_sm = mxGetScalar(prhs[7]);
    real_t insertion_index_deviation_allowance = mxGetScalar(prhs[8]);

    if (mex_assert_1d_array_dimensions(prhs[9], n_u))
        mexErrMsgTxt("Input u_min must have n_u elements.");
    real_t *u_min = mxGetSingles(prhs[9]);
    if (mex_assert_1d_array_dimensions(prhs[10], n_u))
        mexErrMsgTxt("Input u_max must have n_u elements.");
    real_t *u_max = mxGetSingles(prhs[10]);

    if (mex_assert_1d_array_dimensions(prhs[11], N))
        mexErrMsgTxt("Input x1_ref must have N elements.");
    real_t *x1_ref = mxGetSingles(prhs[11]);
    real_t x2_ref = mxGetScalar(prhs[12]);

    if (mex_assert_3d_array_dimensions(prhs[13], n_x, n_x, N))
        mexErrMsgTxt("Input A must be n_x by n_x by N.");
    real_t *A = mxGetSingles(prhs[13]);
    if (mex_assert_3d_array_dimensions(prhs[14], n_u, n_x, N))
        mexErrMsgTxt("Input B must be n_u by n_x by N.");
    real_t *B = mxGetSingles(prhs[14]);
    if (mex_assert_2d_array_dimensions(prhs[15], n_x, N))
        mexErrMsgTxt("Input d must be n_x by N.");
    real_t *d = mxGetSingles(prhs[15]);

    if (mex_assert_1d_array_dimensions(prhs[16], n_x))
        mexErrMsgTxt("Input x0 must have n_x elements.");
    real_t *x0 = mxGetSingles(prhs[16]);

    mwSize x_dims[] = {(mwSize)n_x, (mwSize)N};
    plhs[0] = mxCreateNumericArray(2, x_dims, mxSINGLE_CLASS, mxREAL);
    real_t *x = mxGetSingles(plhs[0]);

    mwSize u_dims[] = {(mwSize)n_u, (mwSize)N};
    plhs[1] = mxCreateNumericArray(2, u_dims, mxSINGLE_CLASS, mxREAL);
    real_t *u = mxGetSingles(plhs[1]);

    if (!initialized) {
        sdqp_lmpc_mmc_init(n_x, n_u, N, q1, q2, x_min, x_max, n_sm, insertion_index_deviation_allowance, u_min, u_max);
        // initialized = 1;
    }
    int err = sdqp_lmpc_mmc_solve(n_x, n_u, N, x1_ref, x2_ref, 
            CAST_CONST_3D_VLA(A, n_x, n_x), CAST_CONST_3D_VLA(B, n_x, n_u), CAST_CONST_2D_VLA(d, n_x), x0, 
            CAST_2D_VLA(x, n_x), CAST_2D_VLA(u, n_u));
    if (err) {
        mexErrMsgTxt("Error while solving.");
    }
    mexAtExit(cleanup);
}
 
