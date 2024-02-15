#include "MexUtils.h"

#include "matrix.h"

int mex_assert_1d_array_dimensions(const mxArray *arr, mwSize d) {
    if ((mxGetM(arr) == 1 && mxGetN(arr) == d) || (mxGetM(arr) == d && mxGetN(arr) == 1)) {
        return 0;
    }
    return 1;
}

int mex_assert_2d_array_dimensions(const mxArray *arr, mwSize d1, mwSize d2) {
    if (mxGetNumberOfDimensions(arr) != 2) 
        return -1;
    const mwSize *dimensions = mxGetDimensions(arr);
    if (dimensions[0] != d1) {
        return 1;
    }
    if (dimensions[1] != d2) {
        return 2;
    }
    return 0;
}

int mex_assert_3d_array_dimensions(const mxArray *arr, mwSize d1, mwSize d2, mwSize d3) {
    if (mxGetNumberOfDimensions(arr) != 3) {
        if (mxGetNumberOfDimensions(arr) == 2 && d3 == 1) {
            return mex_assert_2d_array_dimensions(arr, d1, d2);
        } else {
            return -1;
        }
    }
    const mwSize *dimensions = mxGetDimensions(arr);
    if (dimensions[0] != d1) {
        return 1;
    } else if (dimensions[1] != d2) {
        return 2;
    } else if (dimensions[2] != d3) {
        return 3;
    }
    return 0;
}
