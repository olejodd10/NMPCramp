#include "matrix.h"
#include "mex.h"

#include <stddef.h>
#include <math.h>

#include "MexUtils.h"
#include "Utils.h"
#include "LtvMpc.h"
#include "Types.h"
#include "Ramp.h"

#define NRHS 21
#define NLHS 2

static int initialized = 0;

static void cleanup(void) {
    ltv_mpc_cleanup();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs!=NRHS) {
        mexErrMsgTxt("Incorrect number of inputs.");
    }
    if(nlhs!=NLHS) {
        mexErrMsgTxt("Incorrect number of outputs.");
    }

    for (int i = 0; i < NRHS; ++i) {
        if( !mxIsDouble(prhs[i]) || 
             mxIsComplex(prhs[i])) {
            mexErrMsgTxt("Input must be type double.");
        }
    }

    // Note that transposed matrices are needed because of MATLABs column-major memory layout.
    size_t n_x = (size_t)round(mxGetScalar(prhs[0]));
    size_t n_u = (size_t)round(mxGetScalar(prhs[1]));
    size_t n_y = (size_t)round(mxGetScalar(prhs[2]));
    size_t n_t = (size_t)round(mxGetScalar(prhs[3]));
    size_t N = (size_t)round(mxGetScalar(prhs[4]));

    if (mex_assert_2d_array_dimensions(prhs[5], n_x, n_x))
        mexErrMsgTxt("Input Q must be n_x by n_x.");
    real_t *Q = mxGetDoubles(prhs[5]);
    if (mex_assert_2d_array_dimensions(prhs[6], n_x, n_x))
        mexErrMsgTxt("Input S must be n_x by n_x.");
    real_t *S = mxGetDoubles(prhs[6]);
    if (mex_assert_2d_array_dimensions(prhs[7], n_u, n_u))
        mexErrMsgTxt("Input R must be n_u by n_u.");
    real_t *R = mxGetDoubles(prhs[7]);
    if (mex_assert_2d_array_dimensions(prhs[8], n_x, N))
        mexErrMsgTxt("Input fx must be n_x by N.");
    real_t *fx = mxGetDoubles(prhs[8]);
    if (mex_assert_2d_array_dimensions(prhs[9], n_u, N))
        mexErrMsgTxt("Input fu must be n_u by N.");
    real_t *fu = mxGetDoubles(prhs[9]);

    if (mex_assert_3d_array_dimensions(prhs[10], n_x, n_x, N))
        mexErrMsgTxt("Input A must be n_x by n_x by N.");
    real_t *A = mxGetDoubles(prhs[10]);
    if (mex_assert_3d_array_dimensions(prhs[11], n_u, n_x, N))
        mexErrMsgTxt("Input B must be n_u by n_x by N.");
    real_t *B = mxGetDoubles(prhs[11]);
    if (mex_assert_2d_array_dimensions(prhs[12], n_x, N))
        mexErrMsgTxt("Input d must be n_x by N.");
    real_t *d = mxGetDoubles(prhs[12]);
    if (mex_assert_2d_array_dimensions(prhs[13], n_x, n_y))
        mexErrMsgTxt("Input C must be n_x by n_y.");
    real_t *C = mxGetDoubles(prhs[13]);

    if (mex_assert_1d_array_dimensions(prhs[14], n_y))
        mexErrMsgTxt("Input y_min must have n_y elements.");
    real_t *y_min = mxGetDoubles(prhs[14]);
    if (mex_assert_1d_array_dimensions(prhs[15], n_y))
        mexErrMsgTxt("Input y_max must have n_y elements.");
    real_t *y_max = mxGetDoubles(prhs[15]);

    if (mex_assert_2d_array_dimensions(prhs[16], n_x, n_t))
        mexErrMsgTxt("Input Lt must be n_x by n_t.");
    real_t *Lt = mxGetDoubles(prhs[16]);
    if (mex_assert_1d_array_dimensions(prhs[17], n_t))
        mexErrMsgTxt("Input lt must have n_t elements.");
    real_t *lt = mxGetDoubles(prhs[17]);

    if (mex_assert_1d_array_dimensions(prhs[18], n_u))
        mexErrMsgTxt("Input u_min must have n_u elements.");
    real_t *u_min = mxGetDoubles(prhs[18]);
    if (mex_assert_1d_array_dimensions(prhs[19], n_u))
        mexErrMsgTxt("Input u_max must have n_u elements.");
    real_t *u_max = mxGetDoubles(prhs[19]);

    if (mex_assert_1d_array_dimensions(prhs[20], n_x))
        mexErrMsgTxt("Input x0 must have n_x elements.");
    real_t *x0 = mxGetDoubles(prhs[20]);

    plhs[0] = mxCreateDoubleMatrix((mwSize)n_x, (mwSize)N, mxREAL);
    real_t *x = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix((mwSize)n_u, (mwSize)N, mxREAL);
    real_t *u = mxGetDoubles(plhs[1]);

    if (!initialized) {
        ltv_mpc_init(n_x, n_u, n_y, n_t, N, 
                CAST_CONST_2D_VLA(Q, n_x), CAST_CONST_2D_VLA(S, n_x), CAST_CONST_2D_VLA(R, n_u), 
                CAST_CONST_2D_VLA(fx, n_x), CAST_CONST_2D_VLA(fu, n_u),
                CAST_CONST_2D_VLA(C, n_x), 
                y_min, y_max, CAST_CONST_2D_VLA(Lt, n_x), lt, u_min, u_max);
        initialized = 1;
    }
    int err = ltv_mpc_solve(n_x, n_u, N, 
            CAST_CONST_3D_VLA(A, n_x, n_x), CAST_CONST_3D_VLA(B, n_x, n_u), CAST_CONST_2D_VLA(d, n_x),
            x0, CAST_2D_VLA(x, n_x), CAST_2D_VLA(u, n_u));
    switch (err) {
        case RAMP_ERROR_INFEASIBLE:
            mexErrMsgTxt("Problem infeasible.");
            break;
        case RAMP_ERROR_RANK_2_UPDATE:
            mexErrMsgTxt("Unable to perform rank 2 update.");
            break;
        default:
            break;
    }
    mexAtExit(cleanup);
}
 
