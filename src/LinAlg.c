#include "LinAlg.h"

#include <stddef.h>

#include "Types.h"

real_t linalg_vector_inner_product(size_t n, const real_t v1[n], const real_t v2[n]) {
    real_t result = 0.0;
    for (size_t i = 0; i < n; ++i) {
        result += v1[i]*v2[i];
    }
    return result;
}

void linalg_vector_negate(size_t n, const real_t vec[n], real_t res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = -vec[i];
    }
}

void linalg_vector_scale(size_t n, const real_t vec[n], real_t c, real_t res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = c*vec[i];
    }
}

void linalg_vector_add(size_t n, const real_t v1[n], const real_t v2[n], real_t res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = v1[i] + v2[i];
    }
}

void linalg_vector_add_scaled(size_t n, const real_t v1[n], const real_t v2[n], real_t c, real_t res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = v1[i] + c*v2[i];
    }
}

void linalg_matrix_vector_product(size_t m, size_t n, const real_t mat[m][n], const real_t vec[n], real_t res[m]) {
    for (size_t i = 0; i < m; ++i) {
        res[i] = linalg_vector_inner_product(n, mat[i], vec);
    }
}

void linalg_matrix_product(size_t m, size_t n, size_t p, const real_t m1[m][n], const real_t m2t[p][n], real_t res[m][p]) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < p; ++j) {
            res[i][j] = linalg_vector_inner_product(n, m1[i], m2t[j]);
        }
    }
}

void linalg_matrix_negate(size_t m, size_t n, const real_t mat[m][n], real_t res[m][n]) {
    for (size_t i = 0; i < m; ++i) {
        linalg_vector_negate(n, mat[i], res[i]);
    }
}

void linalg_matrix_transpose(size_t m, size_t n, const real_t mat[m][n], real_t res[n][m]) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            res[j][i] = mat[i][j];
        }
    }
}
