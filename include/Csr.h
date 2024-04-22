#pragma once

#include <stddef.h>

#include "Types.h"

typedef struct {
    real_t *values;
    size_t *col_indices;
    size_t *row_starts;
} csr_t;

void csr_init(size_t m, size_t n, const real_t mat[m][n], csr_t *csr);

void csr_destroy(csr_t *csr);

real_t csr_vector_multiply(size_t n, const csr_t *csr, size_t i, const real_t vec[n]);

void csr_matrix_multiply(size_t m, size_t n, const csr_t *csr, const real_t vec[n], real_t res[m]);
