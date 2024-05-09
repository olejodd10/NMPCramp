#pragma once

#include <stddef.h>

#include "OrderedSet.h"
#include "IndexedVectors.h"
#include "Types.h"

real_t lincomb_element(size_t n_a, const ordered_set_t* a_set, const real_t coefficients[2*n_a], const indexed_vectors_t *cols, size_t i);

void lincomb_compute(size_t n_H, size_t n_a, const ordered_set_t* a_set, const real_t coefficients[2*n_a], const indexed_vectors_t *cols, real_t res[n_H]);

void lincomb_add_scaled(size_t n_a, const ordered_set_t* a_set, const real_t coefs1[2*n_a], const real_t coefs2[2*n_a], real_t factor, real_t res[2*n_a]);

void lincomb_scale(size_t n_a, const ordered_set_t* a_set, const real_t coefficients[2*n_a], real_t factor, real_t res[2*n_a]);
