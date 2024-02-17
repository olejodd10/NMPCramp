#include "Ramp.h"

#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include "Types.h"
#include "IterableSet.h"
#include "IndexedVectors.h"
#include "LinAlg.h"

#define RAMP_EPS 1e-8

static uint8_t infeasiblity_error_enabled = 0;
static real_t infeasibility_error_min;
static real_t infeasibility_error_max;

static real_t *m_column_M4;
static real_t *m_v;

static void (*m_get_column_M4)(size_t, real_t*);

void ramp_init(size_t n_H, void (*get_column_M4)(size_t, real_t*)) {
    m_column_M4 = (real_t*)malloc(sizeof(real_t)*n_H);
    m_v = (real_t*)malloc(sizeof(real_t)*n_H);
    
    m_get_column_M4 = get_column_M4;
}

void ramp_cleanup(void) {
	free(m_column_M4);
	free(m_v);
}

void ramp_enable_infeasibility_error(real_t min, real_t max) {
    infeasiblity_error_enabled = 1;
    infeasibility_error_min = min;
    infeasibility_error_max = max;
}

void ramp_disable_infeasibility_error(void) {
    infeasiblity_error_enabled = 0;
}

// Returns n_H if none found
static inline size_t most_negative_index(size_t n_H, size_t n_a, const iterable_set_t* a_set, real_t y[n_H]) {
    real_t min = -RAMP_EPS;
    size_t index = n_H; // Invalid index, think of it as -1 but using an unsigned data type for efficiency
    for (size_t i = 0; i < n_a; ++i) {
        if (y[i] < min && !iterable_set_contains(a_set, i)) {
            min = y[i];
            index = i;
        }
    }
    for (size_t i = n_a; i < n_H; ++i) {
        if (y[i] < min && iterable_set_contains(a_set, i)) {
            min = y[i];
            index = i;
        }
    }
    return index;
}

// Returns n_H if none found
static inline size_t most_positive_index(size_t n_H, size_t n_a, const iterable_set_t* a_set, real_t y[n_H]) {
    real_t max = RAMP_EPS;
    size_t index = n_H; // Invalid index, think of it as -1 but using an unsigned data type for efficiency
    for (size_t i = 0; i < n_a; ++i) {
        if (y[i] > max && iterable_set_contains(a_set, i)) {
            max = y[i];
            index = i;
        }    
    }
    for (size_t i = n_a; i < n_H; ++i) {
        if (y[i] > max && !iterable_set_contains(a_set, i)) {
            max = y[i];
            index = i;
        }    
    }
    return index;
}

static int compute_v(size_t n_H, const iterable_set_t* a_set, const indexed_vectors_t *invq, size_t index, real_t q0, real_t v[n_H]) {
    // Compute matrix vector product
    // Sparse part
    m_get_column_M4(index, m_column_M4);
    for (size_t i = 0; i < n_H; ++i) {
        v[i] = iterable_set_contains(a_set, i) ? 0.0 : m_column_M4[i]; // 0.0 because the "dense part" computation takes care of the value
    }
    // Dense part
    for (size_t i = iterable_set_first(a_set); i != iterable_set_end(a_set); i = iterable_set_next(a_set, i)) {
        linalg_vector_add_scaled(n_H, v, indexed_vectors_get(invq, i), m_column_M4[i], v);
    }
    // At this point v is defined as in the paper
    real_t qdiv = q0+v[index];
    if (infeasiblity_error_enabled && 
        (fabs(qdiv) < infeasibility_error_min || 
         fabs(qdiv) > infeasibility_error_max)) {
        return RAMP_ERROR_INFEASIBLE;
    }
    linalg_vector_scale(n_H, v, -1.0/qdiv, v);
    return 0;
}

static size_t active_constraints(const iterable_set_t *a_set, size_t n_a) {
    size_t active_b_constraints = iterable_set_partition(a_set);
    size_t inactive_a_constraints = iterable_set_size(a_set) - active_b_constraints;
    size_t active_a_constraints = n_a - inactive_a_constraints;
    return active_a_constraints + active_b_constraints;
}

// Returns n_H if none found
static inline size_t rank_2_update_removal_index(size_t n_H, size_t n_a, const iterable_set_t* a_set, const indexed_vectors_t *invq, size_t i, const real_t y[n_H]) {
    real_t max = 0.0;
    size_t index = n_H;
    m_get_column_M4(i, m_column_M4);
    for (size_t j = 0; j < n_H; ++j) {
        if (j < n_a && iterable_set_contains(a_set, j) || j >= n_a && !iterable_set_contains(a_set, j)) {
            continue;
        }
        real_t divisor = 0.0;
        for (size_t k = iterable_set_first(a_set); k != iterable_set_end(a_set); k = iterable_set_next(a_set, k)) {
            // Note that the order of indices for neg_g_invh_gt doesn't matter since it's symmetric
            divisor += indexed_vectors_get(invq, k)[j] * m_column_M4[k];
        }
        if ((divisor < -RAMP_EPS) && (y[j]/divisor > max || index == n_H)) {
            max = y[j]/divisor;
            index = j;
        }
    }
    return index;
}

static inline void update_invq(size_t n_H, size_t index, const iterable_set_t* a_set, const real_t v[n_H], indexed_vectors_t *invq) {
    for (size_t i = iterable_set_first(a_set); i != iterable_set_end(a_set); i = iterable_set_next(a_set, i)) {
        real_t* invqi = indexed_vectors_get_mut(invq, i);
        linalg_vector_add_scaled(n_H, invqi, v, invqi[index], invqi);
    }
}

static inline void update_y(size_t n_H, size_t index, const real_t v[n_H], real_t y[n_H]) {
    linalg_vector_add_scaled(n_H, y, v, y[index], y);
}

static int active_set_remove(size_t n_H, size_t index, iterable_set_t *a_set, indexed_vectors_t *invq, real_t y[n_H]) {
    int err = compute_v(n_H, a_set, invq, index, 1.0, m_v);
    if (err) {
        return err;
    }
    iterable_set_remove(a_set, index);
    indexed_vectors_remove(invq, index);
    update_y(n_H, index, m_v, y);
    update_invq(n_H, index, a_set, m_v, invq);
    return 0;
}

static int active_set_insert(size_t n_H, size_t index, iterable_set_t *a_set, indexed_vectors_t *invq, real_t y[n_H]) {
    int err = compute_v(n_H, a_set, invq, index, -1.0, m_v);
    if (err) {
        return err;
    }
    update_y(n_H, index, m_v, y);
    update_invq(n_H, index, a_set, m_v, invq);
    iterable_set_insert(a_set, index);
    indexed_vectors_insert(invq, index, m_v);
    indexed_vectors_get_mut(invq, index)[index] += 1.0; // Pretend there was a unit vector in the column to start with
    return 0;
}

static int algorithm1(size_t n_H, size_t n_a, iterable_set_t *a_set, indexed_vectors_t *invq, real_t y[n_H]) {
    while (1) {
        size_t index = most_negative_index(n_H, n_a, a_set, y);
        if (index != n_H) {
            if (index < n_a) {
                int err = active_set_insert(n_H, index, a_set, invq, y);
                if (err) {
                    return err;
                }
            } else {
                int err = active_set_remove(n_H, index, a_set, invq, y);
                if (err) {
                    return err;
                }
            }
        } else {
            index = most_positive_index(n_H, n_a, a_set, y);
            if (index == n_H) {
                break;
            }

            if (active_constraints(a_set, n_a) == n_a) {
                size_t index2 = rank_2_update_removal_index(n_H, n_a, a_set, invq, index, y);
                if (index2 == n_H) {
                    return RAMP_ERROR_RANK_2_UPDATE;
                } else if (index2 < n_a) {
                    int err = active_set_insert(n_H, index2, a_set, invq, y);
                    if (err) {
                        return err;
                    }
                } else {
                    int err = active_set_remove(n_H, index2, a_set, invq, y);
                    if (err) {
                        return err;
                    }
                }
            }

            if (index < n_a) {
                int err = active_set_remove(n_H, index, a_set, invq, y);
                if (err) {
                    return err;
                }
            } else {
                int err = active_set_insert(n_H, index, a_set, invq, y);
                if (err) {
                    return err;
                }
            }
        }
    }
    return 0;
}

int ramp_solve(size_t n_H, size_t n_a, iterable_set_t *a_set, indexed_vectors_t *invq, real_t y[n_H]) {
    return algorithm1(n_H, n_a, a_set, invq, y);
}
