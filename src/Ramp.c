#include "Ramp.h"

#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Types.h"
#include "OrderedSet.h"
#include "IndexedVectors.h"
#include "LinAlg.h"

#define RAMP_EPS 1e-8

static uint8_t infeasiblity_error_enabled = 0;
static real_t infeasibility_error_min;
static real_t infeasibility_error_max;

static indexed_vectors_t m_invq;
static ordered_set_t m_a_set;

static indexed_vectors_t m_m4_cols;

static real_t *m_v;

static void (*m_get_column_M4)(size_t, real_t*);

static const real_t* get_column_M4(size_t index) {
    real_t* column_M4 = indexed_vectors_get_mut(&m_m4_cols, index);
    if (column_M4 == NULL) {
        indexed_vectors_insert(&m_m4_cols, index);
        column_M4 = indexed_vectors_get_mut(&m_m4_cols, index);
        m_get_column_M4(index, column_M4);
    } 
    return column_M4;
}

void ramp_init(size_t n_H, size_t n_a, void (*get_column_M4)(size_t, real_t*)) {
    indexed_vectors_init(&m_invq, 2*n_a, n_H, n_H);
    ordered_set_init(&m_a_set, n_H, n_a);

    indexed_vectors_init(&m_m4_cols, 2*n_a + 1, n_H, n_H);

    m_v = (real_t*)malloc(sizeof(real_t)*n_H);
    
    m_get_column_M4 = get_column_M4;
}

void ramp_cleanup(void) {
    indexed_vectors_destroy(&m_invq);
    ordered_set_destroy(&m_a_set);

    indexed_vectors_destroy(&m_m4_cols);

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
static inline size_t most_incorrect_active_constraint(size_t n_H, size_t n_a, const ordered_set_t* a_set, const real_t y[n_H]) {
    real_t min = -RAMP_EPS;
    size_t index = n_H; // Invalid index, think of it as -1 but using an unsigned data type for efficiency
    for (size_t i = 0; i < n_a; ++i) {
        if (-y[i] < min && !ordered_set_contains(a_set, i)) {
            min = -y[i];
            index = i;
        }
    }
    for (size_t i = n_a; i < n_H; ++i) {
        if (y[i] < min && ordered_set_contains(a_set, i)) {
            min = y[i];
            index = i;
        }
    }
    return index;
}

// Returns n_H if none found
static inline size_t most_incorrect_inactive_constraint(size_t n_H, size_t n_a, const ordered_set_t* a_set, const real_t y[n_H]) {
    real_t max = RAMP_EPS;
    size_t index = n_H; // Invalid index, think of it as -1 but using an unsigned data type for efficiency
    for (size_t i = 0; i < n_a; ++i) {
        if (-y[i] > max && ordered_set_contains(a_set, i)) {
            max = -y[i];
            index = i;
        }    
    }
    for (size_t i = n_a; i < n_H; ++i) {
        if (y[i] > max && !ordered_set_contains(a_set, i)) {
            max = y[i];
            index = i;
        }    
    }
    return index;
}

static int compute_v(size_t n_H, const ordered_set_t* a_set, const indexed_vectors_t *invq, size_t index, real_t q0, real_t v[n_H]) {
    // Compute matrix vector product
    // Sparse part
    const real_t *column_M4 = get_column_M4(index);
    for (size_t i = 0; i < n_H; ++i) {
        v[i] = ordered_set_contains(a_set, i) ? 0.0 : column_M4[i]; // 0.0 because the "dense part" computation takes care of the value
    }
    // Dense part
    for (size_t n = 0, i = ordered_set_nth(a_set, n); i != ordered_set_end(a_set); i = ordered_set_nth(a_set, ++n)) {
        linalg_vector_add_scaled(n_H, v, indexed_vectors_get(invq, i), column_M4[i], v);
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

static size_t active_constraints(const ordered_set_t *a_set, size_t n_a) {
    size_t active_b_constraints = ordered_set_partition(a_set);
    size_t inactive_a_constraints = ordered_set_size(a_set) - active_b_constraints;
    size_t active_a_constraints = n_a - inactive_a_constraints;
    return active_a_constraints + active_b_constraints;
}

// Returns n_H if none found
static inline size_t rank_2_update_removal_index(size_t n_H, size_t n_a, const ordered_set_t* a_set, const indexed_vectors_t *invq, size_t i, const real_t y[n_H]) {
    real_t min = -RAMP_EPS;
    size_t index = n_H;
    const real_t *column_M4 = get_column_M4(i);
    for (size_t j = 0; j < n_a; ++j) {
        if (ordered_set_contains(a_set, j)) {
            continue;
        }
        real_t numerator = column_M4[j];
        for (size_t n = 0, k = ordered_set_nth(a_set, n); k != ordered_set_end(a_set); k = ordered_set_nth(a_set, ++n)) {
            numerator += indexed_vectors_get(invq, k)[j] * column_M4[k];
        }
        real_t val = ordered_set_contains(a_set, i) ? -numerator/y[j] : numerator/y[j];
        if (val < min) {
            min = val;
            index = j;
        }
    }
    for (size_t j = n_a; j < n_H; ++j) {
        if (!ordered_set_contains(a_set, j)) {
            continue;
        }
        real_t numerator = 0.0;
        for (size_t n = 0, k = ordered_set_nth(a_set, n); k != ordered_set_end(a_set); k = ordered_set_nth(a_set, ++n)) {
            numerator += indexed_vectors_get(invq, k)[j] * column_M4[k];
        }
        real_t val = ordered_set_contains(a_set, i) ? -numerator/y[j] : numerator/y[j];
        if (val < min) {
            min = val;
            index = j;
        }
    }
    return index;
}

static inline void update_invq(size_t n_H, size_t index, const ordered_set_t* a_set, const real_t v[n_H], indexed_vectors_t *invq) {
    for (size_t n = 0, i = ordered_set_nth(a_set, n); i != ordered_set_end(a_set); i = ordered_set_nth(a_set, ++n)) {
        real_t* invqi = indexed_vectors_get_mut(invq, i);
        linalg_vector_add_scaled(n_H, invqi, v, invqi[index], invqi);
    }
}

static inline void update_y(size_t n_H, size_t index, const real_t v[n_H], real_t y[n_H]) {
    linalg_vector_add_scaled(n_H, y, v, y[index], y);
}

static int active_set_remove(size_t n_H, size_t index, ordered_set_t *a_set, indexed_vectors_t *invq, real_t y[n_H]) {
    int err = compute_v(n_H, a_set, invq, index, 1.0, m_v);
    if (err) {
        return err;
    }
    ordered_set_remove(a_set, index);
    indexed_vectors_remove(invq, index);
    update_y(n_H, index, m_v, y);
    update_invq(n_H, index, a_set, m_v, invq);
    indexed_vectors_remove(&m_m4_cols, index);
    return 0;
}

static int active_set_insert(size_t n_H, size_t index, ordered_set_t *a_set, indexed_vectors_t *invq, real_t y[n_H]) {
    int err = compute_v(n_H, a_set, invq, index, -1.0, m_v);
    if (err) {
        return err;
    }
    update_y(n_H, index, m_v, y);
    update_invq(n_H, index, a_set, m_v, invq);
    ordered_set_insert(a_set, index);
    indexed_vectors_insert(invq, index);
    memcpy(indexed_vectors_get_mut(invq, index), m_v, sizeof(real_t)*n_H);
    indexed_vectors_get_mut(invq, index)[index] += 1.0; // Pretend there was a unit vector in the column to start with
    return 0;
}

static int algorithm1(size_t n_H, size_t n_a, ordered_set_t *a_set, indexed_vectors_t *invq, real_t y[n_H]) {
    while (1) {
        size_t index = most_incorrect_active_constraint(n_H, n_a, a_set, y);
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
            index = most_incorrect_inactive_constraint(n_H, n_a, a_set, y);
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

int ramp_solve(size_t n_H, size_t n_a, int hotstart_variant, real_t y[n_H]) {
    switch (hotstart_variant) {
        case RAMP_HOTSTART_M4_UNCHANGED:
            memcpy(m_v, y, sizeof(real_t)*n_H); // Use m_v temporarily
            for (size_t n = 0, i = ordered_set_nth(&m_a_set, n); i != ordered_set_end(&m_a_set); i = ordered_set_nth(&m_a_set, ++n)) {
                linalg_vector_add_scaled(n_H, y, indexed_vectors_get(&m_invq, i), m_v[i], y);
                y[i] -= m_v[i];
            }
            break;
        default:
            ordered_set_clear(&m_a_set);
            indexed_vectors_clear(&m_invq);
            indexed_vectors_clear(&m_m4_cols);
            break;
    }
    return algorithm1(n_H, n_a, &m_a_set, &m_invq, y);
}
