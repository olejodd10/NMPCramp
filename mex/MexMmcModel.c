#include "matrix.h"
#include "mex.h"

#include <stddef.h>
#include <math.h>

#include "MexUtils.h"
#include "Utils.h"
#include "MmcModel.h"
#include "Types.h"

#define NRHS 15
#define NLHS 3

static int initialized = 0;

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

    real_t R = mxGetScalar(prhs[0]);
    real_t Rc = mxGetScalar(prhs[1]);
    real_t L = mxGetScalar(prhs[2]);
    real_t Lc = mxGetScalar(prhs[3]);
    real_t C = mxGetScalar(prhs[4]);
    real_t Ts = mxGetScalar(prhs[5]);
    real_t n_sm = mxGetScalar(prhs[6]);
    size_t N = (size_t)round(mxGetScalar(prhs[7]));

    if (mex_assert_2d_array_dimensions(prhs[8], N_X, N))
        mexErrMsgTxt("Input x must be n_x by N.");
    real_t *x = mxGetSingles(prhs[8]);
    if (mex_assert_2d_array_dimensions(prhs[9], N_U, N))
        mexErrMsgTxt("Input u must be n_u by N.");
    real_t *u = mxGetSingles(prhs[9]);

    if (mex_assert_1d_array_dimensions(prhs[10], N))
        mexErrMsgTxt("Input vf must have N elements.");
    real_t *vf = mxGetSingles(prhs[10]);
    real_t Vdc = mxGetScalar(prhs[11]);

    if (mex_assert_3d_array_dimensions(prhs[12], N_X, N_X, N))
        mexErrMsgTxt("Input A must be n_x by n_x by N.");
    plhs[0] = mxDuplicateArray(prhs[12]);
    real_t *A = mxGetSingles(plhs[0]);
    if (mex_assert_3d_array_dimensions(prhs[13], N_U, N_X, N))
        mexErrMsgTxt("Input B must be n_u by n_x by N.");
    plhs[1] = mxDuplicateArray(prhs[13]);
    real_t *B = mxGetSingles(plhs[1]);
    if (mex_assert_2d_array_dimensions(prhs[14], N_X, N))
        mexErrMsgTxt("Input d must be n_x by N.");
    plhs[2] = mxDuplicateArray(prhs[14]);
    real_t *d = mxGetSingles(plhs[2]);

    if (!initialized) {
        mmc_model_set_parameters(R, Rc, L, Lc, C, Ts, n_sm);
        mmc_model_get_init(N, CAST_3D_VLA(A, N_X, N_X), CAST_3D_VLA(B, N_X, N_U), CAST_2D_VLA(d, N_X));
        initialized = 1;
    }
    mmc_model_get(N, CAST_CONST_2D_VLA(x, N_X), CAST_CONST_2D_VLA(u, N_U), vf, Vdc, CAST_3D_VLA(A, N_X, N_X), CAST_3D_VLA(B, N_X, N_U), CAST_2D_VLA(d, N_X));
}
 
