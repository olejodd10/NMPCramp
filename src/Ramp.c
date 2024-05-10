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
#include "LinComb.h"

#define RAMP_EPS 1e-8

static uint8_t infeasiblity_error_enabled = 0;
static real_t infeasibility_error_min;
static real_t infeasibility_error_max;

static ordered_set_t m_a_set;

static indexed_vectors_t m_m4_cols;
static indexed_vectors_t m_invq_coefs;

static size_t m_n_H;

static real_t *m_v_coefs;
static real_t *m_y_coefs;
static char *m_temp;

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
    ordered_set_init(&m_a_set, n_H, n_a);

    indexed_vectors_init(&m_m4_cols, 2*n_a + 1, n_H, n_H);
    indexed_vectors_init(&m_invq_coefs, 2*n_a, 2*n_a, n_H);

    m_n_H = n_H;

    m_v_coefs = (real_t*)malloc(sizeof(real_t)*2*n_a);
    m_y_coefs = (real_t*)malloc(sizeof(real_t)*2*n_a);
    m_temp = (char*)malloc(sizeof(char)*n_H);
    
    m_get_column_M4 = get_column_M4;
}

void ramp_cleanup(void) {
    ordered_set_destroy(&m_a_set);

    indexed_vectors_destroy(&m_m4_cols);
    indexed_vectors_destroy(&m_invq_coefs);
    free(m_v_coefs);
    free(m_y_coefs);
    free(m_temp);
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

static int compute_v(size_t n_a, const ordered_set_t* a_set, const indexed_vectors_t *invq_coefs, size_t index, real_t q0, real_t v_coefs[2*n_a]) {
    // Compute matrix vector product
    // Sparse part
    const real_t *column_M4 = get_column_M4(index);
    memset(v_coefs, 0, sizeof(real_t)*2*n_a);
    size_t lim = ordered_set_size(a_set);
    for (size_t n = 0; n < lim; ++n) {
        size_t i = ordered_set_nth(a_set, n);
        lincomb_add_scaled(n_a, a_set, v_coefs, indexed_vectors_get(invq_coefs, i), column_M4[i], v_coefs);
    }
    // On inserts the new index will not be inserted yet, so we need to compensate manually
    real_t vi = lincomb_element(n_a, a_set, v_coefs, &m_m4_cols, index) + column_M4[index];
    real_t qdiv = q0 + vi;
    if (infeasiblity_error_enabled && 
        (fabs(qdiv) < infeasibility_error_min || 
         fabs(qdiv) > infeasibility_error_max)) {
        return RAMP_ERROR_INFEASIBLE;
    }
    lincomb_scale(n_a, a_set, v_coefs, -1.0/qdiv, v_coefs);
    if (ordered_set_contains(a_set, index)) {
        v_coefs[ordered_set_whereis(a_set, index)] += -1.0/qdiv;
    } else {
        v_coefs[ordered_set_size(a_set)] += -1.0/qdiv;
    }
    return 0;
}

static size_t active_constraints(const ordered_set_t *a_set, size_t n_a) {
    size_t active_b_constraints = ordered_set_partition(a_set);
    size_t inactive_a_constraints = ordered_set_size(a_set) - active_b_constraints;
    size_t active_a_constraints = n_a - inactive_a_constraints;
    return active_a_constraints + active_b_constraints;
}

// Returns n_H if none found
static inline size_t rank_2_update_removal_index(size_t n_H, size_t n_a, const ordered_set_t* a_set, const indexed_vectors_t *invq_coefs, size_t i, const real_t y[n_H]) {
    real_t min = -RAMP_EPS;
    size_t index = n_H;
    const real_t *column_M4 = get_column_M4(i);
    for (size_t j = 0; j < n_a; ++j) {
        if (ordered_set_contains(a_set, j)) {
            continue;
        }
        real_t numerator = column_M4[j];
        for (size_t n = 0, k = ordered_set_nth(a_set, n); k != ordered_set_end(a_set); k = ordered_set_nth(a_set, ++n)) {
            const real_t *invq_coefs_k = indexed_vectors_get(invq_coefs, k);
            real_t invq_kj = lincomb_element(n_a, a_set, invq_coefs_k, &m_m4_cols, j) + (k == j ? 1.0 : 0.0);
            numerator += invq_kj * column_M4[k];
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
            const real_t *invq_coefs_k = indexed_vectors_get(invq_coefs, k);
            real_t invq_kj = lincomb_element(n_a, a_set, invq_coefs_k, &m_m4_cols, j) + (k == j ? 1.0 : 0.0);
            numerator += invq_kj * column_M4[k];
        }
        real_t val = ordered_set_contains(a_set, i) ? -numerator/y[j] : numerator/y[j];
        if (val < min) {
            min = val;
            index = j;
        }
    }
    return index;
}

static inline void update_invq(size_t n_a, size_t index, const ordered_set_t* a_set, const real_t v_coefs[2*n_a], indexed_vectors_t *invq_coefs) {
    size_t lim = ordered_set_size(a_set);
    for (size_t n = 0; n < lim; ++n) {
        size_t i = ordered_set_nth(a_set, n);
        if (i == index) {
            continue;
        }
        real_t *invq_coefs_i = indexed_vectors_get_mut(invq_coefs, i);
        real_t factor = lincomb_element(n_a, a_set, invq_coefs_i, &m_m4_cols, index); // + (i == index ? 1.0 : 0.0)
        if (!ordered_set_contains(a_set, index)) {
            invq_coefs_i[ordered_set_size(a_set)] = factor*v_coefs[ordered_set_size(a_set)]; // This is irrelevant if index is leaving, and important if index is joining
        }
        lincomb_add_scaled(n_a, a_set, invq_coefs_i, v_coefs, factor, invq_coefs_i);
    }
}

static void compute_y(size_t n_H, size_t n_a, const ordered_set_t* a_set, const real_t y_coefs[2*n_a], const real_t m[n_H], real_t y[n_H]) {
    lincomb_compute(n_H, n_a, a_set, y_coefs, &m_m4_cols, y);
    linalg_vector_add(n_H, y, m, y);
}

static inline void update_y_coefs(size_t n_H, size_t n_a, size_t index, const ordered_set_t* a_set, const real_t v_coefs[2*n_a], real_t y_coefs[2*n_a], const real_t m[n_H], real_t y[n_H]) {
    real_t yi = y[index]; 
    // real_t yi = lincomb_element(n_H, a_set, y_coefs, &m_m4_cols, index) + m[index];
    lincomb_add_scaled(n_a, a_set, y_coefs, v_coefs, yi, y_coefs);
    if (!ordered_set_contains(a_set, index)) {
        y_coefs[ordered_set_size(a_set)] = yi*v_coefs[ordered_set_size(a_set)]; // Irrelevant if removal, important if insertion
    }
}

static int active_set_remove(size_t n_H, size_t n_a, size_t index, ordered_set_t *a_set, indexed_vectors_t *invq_coefs, const real_t m[n_H], real_t y[n_H]) {
    int err = compute_v(n_a, a_set, invq_coefs, index, 1.0, m_v_coefs); // OBS! Krever gammel invq_coefs/invq, mens v_i krever ny active set og kolonne
    if (err) {
        return err;
    }
    update_y_coefs(n_H, n_a, index, a_set, m_v_coefs, m_y_coefs, m, y);
    update_invq(n_a, index, a_set, m_v_coefs, invq_coefs);
    indexed_vectors_remove(invq_coefs, index);
    indexed_vectors_remove(&m_m4_cols, index);

    // v coefs no longer needed, but invq coefs and y coefs need to be swapped
    m_y_coefs[ordered_set_whereis(a_set, index)] = m_y_coefs[ordered_set_size(a_set)-1];
    size_t lim = ordered_set_size(&m_a_set);
    for (size_t n = 0; n < lim; ++n) {
        size_t i = ordered_set_nth(&m_a_set, n);
        if (i == index) {
            continue;
        }
        real_t *invq_coefs_i = indexed_vectors_get_mut(invq_coefs, i);
        invq_coefs_i[ordered_set_whereis(a_set, index)] = invq_coefs_i[ordered_set_size(a_set)-1];
    }

    ordered_set_remove(a_set, index);
    compute_y(n_H, n_a, a_set, m_y_coefs, m, y);
    return 0;
}

static int active_set_insert(size_t n_H, size_t n_a, size_t index, ordered_set_t *a_set, indexed_vectors_t *invq_coefs, const real_t m[n_H], real_t y[n_H]) {
    int err = compute_v(n_a, a_set, invq_coefs, index, -1.0, m_v_coefs);
    if (err) {
        return err;
    }
    update_invq(n_a, index, a_set, m_v_coefs, invq_coefs);
    indexed_vectors_insert(invq_coefs, index); // Create slot so memory can be manipulated
    memcpy(indexed_vectors_get_mut(invq_coefs, index), m_v_coefs, sizeof(real_t)*2*n_a);
    update_y_coefs(n_H, n_a, index, a_set, m_v_coefs, m_y_coefs, m, y);
    ordered_set_insert(a_set, index); // Needed because y has element here, lincomb needs that
    compute_y(n_H, n_a, a_set, m_y_coefs, m, y);
    return 0;
}

static int algorithm1(size_t n_H, size_t n_a, ordered_set_t *a_set, indexed_vectors_t *invq_coefs, const real_t m[n_H], real_t y[n_H]) {
    while (1) {
        size_t index = most_incorrect_active_constraint(n_H, n_a, a_set, y);
        if (index != n_H) {
            if (index < n_a) {
                int err = active_set_insert(n_H, n_a, index, a_set, invq_coefs, m, y);
                if (err) {
                    return err;
                }
            } else {
                int err = active_set_remove(n_H, n_a, index, a_set, invq_coefs, m, y);
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
                size_t index2 = rank_2_update_removal_index(n_H, n_a, a_set, invq_coefs, index, y);
                if (index2 == n_H) {
                    return RAMP_ERROR_RANK_2_UPDATE;
                } else if (index2 < n_a) {
                    int err = active_set_insert(n_H, n_a, index2, a_set, invq_coefs, m, y);
                    if (err) {
                        return err;
                    }
                } else {
                    int err = active_set_remove(n_H, n_a, index2, a_set, invq_coefs, m, y);
                    if (err) {
                        return err;
                    }
                }
            }

            if (index < n_a) {
                int err = active_set_remove(n_H, n_a, index, a_set, invq_coefs, m, y);
                if (err) {
                    return err;
                }
            } else {
                int err = active_set_insert(n_H, n_a, index, a_set, invq_coefs, m, y);
                if (err) {
                    return err;
                }
            }
        }
    }
    return 0;
}

int ramp_solve(size_t n_H, size_t n_a, int hotstart_variant, const real_t m[n_H], real_t y[n_H]) {
    switch (hotstart_variant) {
        case RAMP_HOTSTART_M4_UNCHANGED:
            memset(m_y_coefs, 0, sizeof(real_t)*2*n_a);
            size_t lim = ordered_set_size(&m_a_set);
            for (size_t n = 0; n < lim; ++n) {
                size_t i = ordered_set_nth(&m_a_set, n);
                lincomb_add_scaled(n_a, &m_a_set, m_y_coefs, indexed_vectors_get(&m_invq_coefs, i), m[i], m_y_coefs);
            }
            break;
        case RAMP_HOTSTART_M4_CHANGED:
            memset(m_temp, 0, sizeof(char)*n_H);
            lim = ordered_set_size(&m_a_set);
            for (size_t n = 0; n < lim; ++n) {
                size_t i = ordered_set_nth(&m_a_set, n);
                m_temp[i] = 1;
            }
            indexed_vectors_clear(&m_m4_cols);
            indexed_vectors_clear(&m_invq_coefs);
            ordered_set_clear(&m_a_set);
            for (size_t i = 0; i < n_H; ++i) {
                if (m_temp[i]) {
                    int err = compute_v(n_a, &m_a_set, &m_invq_coefs, i, -1.0, m_v_coefs);
                    if (err) {
                        return err;
                    }
                    update_invq(n_a, i, &m_a_set, m_v_coefs, &m_invq_coefs);
                    indexed_vectors_insert(&m_invq_coefs, i); 
                    memcpy(indexed_vectors_get_mut(&m_invq_coefs, i), m_v_coefs, sizeof(real_t)*2*n_a);
                    real_t yi = lincomb_element(n_a, &m_a_set, m_y_coefs, &m_m4_cols, i) + m[i];
                    lincomb_add_scaled(n_a, &m_a_set, m_y_coefs, m_v_coefs, yi, m_y_coefs);
                    m_y_coefs[ordered_set_size(&m_a_set)] = yi*m_v_coefs[ordered_set_size(&m_a_set)]; // Irrelevant if removal, important if insertion
                    ordered_set_insert(&m_a_set, i); 
                }
            }
            break;
        default:
            indexed_vectors_clear(&m_m4_cols);
            indexed_vectors_clear(&m_invq_coefs);
            ordered_set_clear(&m_a_set);
            break;
    }
    compute_y(n_H, n_a, &m_a_set, m_y_coefs, m, y);
    return algorithm1(n_H, n_a, &m_a_set, &m_invq_coefs, m, y);
}
