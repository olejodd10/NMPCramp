#pragma once

#include <stddef.h>

#include "Types.h"

#define RAMP_ERROR_INFEASIBLE -1 // Infeasibility limits breached
#define RAMP_ERROR_RANK_2_UPDATE -2 // Rank 2 update fail

#define RAMP_HOTSTART_NONE 0 // Reset active set on every iteration
#define RAMP_HOTSTART_M4_UNCHANGED 1 // Keep active set between iterations

/**
 * @brief Initialize ramp solver
 *
 * @param[in] n_H Height of matrix H
 * @param[in] n_a Total number of constraints in a-partition
 * @param[in] get_column_M4 Function pointer to procedure for getting column of M4
 */
void ramp_init(size_t n_H, size_t n_a, void (*get_column_M4)(size_t, real_t*));

/**
 * @brief Cleanup ramp solver
 */
void ramp_cleanup(void);

/**
 * @brief Enable infeasibility errors
 *
 * @param[in] min Lower bound
 * @param[in] max Upper bound
 */
void ramp_enable_infeasibility_error(real_t min, real_t max);

/**
 * @brief Disable infeasibility errors
 */
void ramp_disable_infeasibility_error(void);

/**
 * @brief Solve ramp equation
 *
 * @param[in] n_H Height of matrix H
 * @param[in] n_a Total number of constraints in a-partition
 * @param[in] hotstart_variant Hotstart variant, see RAMP_HOTSTART_...-defines in Ramp.h
 * @param[inout] y m on input, converged y on output
 * @return 0 on success, <0 on error
 */
int ramp_solve(size_t n_H, size_t n_a, int hotstart_variant, real_t y[n_H]);
