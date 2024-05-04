#pragma once

#include <stddef.h>

#include "OrderedSet.h"
#include "IndexedVectors.h"
#include "Types.h"

real_t lincomb_element(size_t n_H, const ordered_set_t* a_set, const real_t coefficients[n_H], const indexed_vectors_t *cols, size_t i);

void lincomb_compute(size_t n_H, const ordered_set_t* a_set, const real_t coefficients[n_H], const indexed_vectors_t *cols, real_t res[n_H]);

void lincomb_add_scaled(size_t n_H, const ordered_set_t* a_set, const real_t coefs1[n_H], const real_t coefs2[n_H], real_t factor, real_t res[n_H]);

void lincomb_scale(size_t n_H, const ordered_set_t* a_set, const real_t coefficients[n_H], real_t factor, real_t res[n_H]);
