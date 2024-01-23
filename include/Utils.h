#pragma once

#include "Types.h"

// Convenience macros to cast pointers to variable length array type, avoiding compiler warnings
#define CAST_VLA(ptr) (real_t(*)[])ptr
#define CAST_CONST_VLA(ptr) (const real_t(*)[])ptr

#define MATRIX_ROW(M, n, i) (&M[n*i])
#define MATRIX_ELEMENT(M, n, i, j) (M[n*i+j])
