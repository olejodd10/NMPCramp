#include "LinComb.h"

#include <stddef.h>
#include <string.h>

#include "OrderedSet.h"
#include "IndexedVectors.h"
#include "LinAlg.h"
#include "Types.h"

real_t lincomb_element(size_t n_a, const ordered_set_t* a_set, const real_t coefficients[2*n_a], const indexed_vectors_t *cols, size_t i) {
    real_t result = 0.0;
    for (size_t j = 0; j < ordered_set_size(a_set); ++j) {
        result += coefficients[j]*indexed_vectors_get(cols, a_set->ordering[j])[i];
    }
    return result;
}

void lincomb_compute(size_t n_H, size_t n_a, const ordered_set_t* a_set, const real_t coefficients[2*n_a], const indexed_vectors_t *cols, real_t res[n_H]) {
    memset(res, 0, sizeof(real_t)*n_H);
    for (size_t i = 0; i < ordered_set_size(a_set); ++i) {
        linalg_vector_add_scaled(n_H, res, indexed_vectors_get(cols, a_set->ordering[i]), coefficients[i], res);
    }
}

void lincomb_add_scaled(size_t n_a, const ordered_set_t* a_set, const real_t coefs1[2*n_a], const real_t coefs2[2*n_a], real_t factor, real_t res[2*n_a]) {
    for (size_t i = 0; i < ordered_set_size(a_set); ++i) {
        res[i] = coefs1[i] + coefs2[i]*factor;
    }
}

void lincomb_scale(size_t n_a, const ordered_set_t* a_set, const real_t coefficients[2*n_a], real_t factor, real_t res[2*n_a]) {
    for (size_t i = 0; i < ordered_set_size(a_set); ++i) {
        res[i] = coefficients[i] * factor;
    }
}
