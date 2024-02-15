#pragma once

#include "mex.h"

int mex_assert_1d_array_dimensions(const mxArray *arr, mwSize d);

int mex_assert_2d_array_dimensions(const mxArray *arr, mwSize d1, mwSize d2);

int mex_assert_3d_array_dimensions(const mxArray *arr, mwSize d1, mwSize d2, mwSize d3);
