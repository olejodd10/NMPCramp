#pragma once

#include <stddef.h>

#include "Types.h"

real_t linalg_vector_inner_product(size_t n, const real_t v1[n], const real_t v2[n]); 

void linalg_vector_negate(size_t n, const real_t vec[n], real_t res[n]); 

void linalg_vector_scale(size_t n, const real_t vec[n], real_t c, real_t res[n]); 

void linalg_vector_add(size_t n, const real_t v1[n], const real_t v2[n], real_t res[n]); 

void linalg_vector_add_scaled(size_t n, const real_t v1[n], const real_t v2[n], real_t c, real_t res[n]); 

void linalg_matrix_vector_product(size_t m, size_t n, const real_t mat[m][n], const real_t vec[n], real_t res[m]); 

void linalg_matrix_product(size_t m, size_t n, size_t p, const real_t m1[m][n], const real_t m2t[p][n], real_t res[m][p]); 

void linalg_matrix_negate(size_t m, size_t n, const real_t mat[m][n], real_t res[m][n]); 

void linalg_matrix_transpose(size_t m, size_t n, const real_t mat[m][n], real_t res[n][m]); 
