#include "LinComb.h"

#include <stddef.h>
#include <string.h>

#include "OrderedSet.h"
#include "IndexedVectors.h"
#include "LinAlg.h"
#include "Types.h"

real_t lincomb_element(size_t n_H, const ordered_set_t* a_set, const real_t coefficients[n_H], const indexed_vectors_t *cols, size_t i) {
    real_t result = 0.0;
    for (size_t n = 0, j = ordered_set_nth(a_set, n); j != ordered_set_end(a_set); j = ordered_set_nth(a_set, ++n)) {
        result += coefficients[j]*indexed_vectors_get(cols, j)[i];
    }
    return result;
}

void lincomb_compute(size_t n_H, const ordered_set_t* a_set, const real_t coefficients[n_H], const indexed_vectors_t *cols, real_t res[n_H]) {
    memset(res, 0, sizeof(real_t)*n_H);
    for (size_t n = 0, j = ordered_set_nth(a_set, n); j != ordered_set_end(a_set); j = ordered_set_nth(a_set, ++n)) {
        linalg_vector_add_scaled(n_H, res, indexed_vectors_get(cols, j), coefficients[j], res);
    }
}

void lincomb_add_scaled(size_t n_H, const ordered_set_t* a_set, const real_t coefs1[n_H], const real_t coefs2[n_H], real_t factor, real_t res[n_H]) {
    for (size_t n = 0, i = ordered_set_nth(a_set, n); i != ordered_set_end(a_set); i = ordered_set_nth(a_set, ++n)) {
        res[i] = coefs1[i] + coefs2[i]*factor;
    }
}

void lincomb_scale(size_t n_H, const ordered_set_t* a_set, const real_t coefficients[n_H], real_t factor, real_t res[n_H]) {
    for (size_t n = 0, i = ordered_set_nth(a_set, n); i != ordered_set_end(a_set); i = ordered_set_nth(a_set, ++n)) {
        res[i] = coefficients[i] * factor;
    }
}

