#include "Csr.h"

#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include "Types.h"

#define CSR_EPS 1e-8

static size_t num_nonzero(size_t m, size_t n, const real_t mat[m][n]) {
    size_t result = 0;
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (fabs(mat[i][j]) > CSR_EPS) {
                ++result;
            }
        }
    }
    return result;
}

void csr_init(size_t m, size_t n, const real_t mat[m][n], csr_t *csr) {
    size_t nnz = num_nonzero(m, n, mat);

    csr->values = (real_t*)malloc(sizeof(real_t)*nnz);
    csr->col_indices = (size_t*)malloc(sizeof(size_t)*nnz);
    csr->row_starts = (size_t*)malloc(sizeof(size_t)*(m+1));

    size_t num_elements = 0;
    for (size_t i = 0; i < m; ++i) {
        csr->row_starts[i] = num_elements;
        for (size_t j = 0; j < n; ++j) {
            if (fabs(mat[i][j]) > CSR_EPS) {
                csr->values[num_elements] = mat[i][j];
                csr->col_indices[num_elements] = j; 
                ++num_elements;
            }
        }
    }
    csr->row_starts[m] = nnz;
}

void csr_destroy(csr_t *csr) {
    free(csr->values);
    free(csr->col_indices);
    free(csr->row_starts);
}

real_t csr_vector_multiply(size_t n, const csr_t *csr, size_t i, const real_t vec[n]) {
    real_t sum = 0.0;
    for (size_t k = csr->row_starts[i]; k < csr->row_starts[i+1]; ++k) {
        size_t j = csr->col_indices[k];
        sum += csr->values[k]*vec[j];
    }
    return sum;
}

void csr_matrix_multiply(size_t m, size_t n, const csr_t *csr, const real_t vec[n], real_t res[m]) {
    for (size_t i = 0; i < m; ++i) {
        res[i] = csr_vector_multiply(n, csr, i, vec);
    }
}
