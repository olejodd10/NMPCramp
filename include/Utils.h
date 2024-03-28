#pragma once

#include "Types.h"

// Convenience macros to cast pointers to variable length array type, avoiding compiler warnings
#define CAST_2D_VLA(ptr, dim2) (real_t(*)[(dim2)])ptr
#define CAST_CONST_2D_VLA(ptr, dim2) (const real_t(*)[(dim2)])ptr

#define CAST_3D_VLA(ptr, dim2, dim3) (real_t(*)[(dim2)][(dim3)])ptr
#define CAST_CONST_3D_VLA(ptr, dim2, dim3) (const real_t(*)[(dim2)][(dim3)])ptr

#define MATRIX_ROW(M, n, i) (&M[n*i])
#define MATRIX_ELEMENT(M, n, i, j) (M[n*i+j])

#define PI 3.14159265358979323846
