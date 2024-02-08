#include "Simulate.h"

#include <stddef.h>

#include "Types.h"
#include "LinAlg.h"

void simulate_lti(size_t n, size_t p, const real_t a[n][n], const real_t x[n], const real_t b[n][p], const real_t u[p], real_t res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = linalg_vector_inner_product(n, a[i], x) + linalg_vector_inner_product(p, b[i], u);
    }
}

void simulate_affine(size_t n, size_t p, const real_t a[n][n], const real_t x[n], const real_t b[n][p], const real_t u[p], const real_t d[n], real_t res[n]) {
    for (size_t i = 0; i < n; ++i) {
        res[i] = linalg_vector_inner_product(n, a[i], x) + linalg_vector_inner_product(p, b[i], u) + d[i];
    }
}
