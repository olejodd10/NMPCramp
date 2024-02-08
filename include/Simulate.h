#pragma once

#include <stddef.h>

#include "Types.h"

/**
 * @brief Calculate the next state of an LTI system.
 *
 * @param[in] n Number of system states
 * @param[in] p Number of system inputs
 * @param[in] a Flattened n by n array representing system matrix
 * @param[in] x Length n array representing current system state
 * @param[in] b Flattened n by p array representing input matrix
 * @param[in] u Length p array representing system inputs
 * @param[out] res Length n array for next system state
 * @warning x and res can not be the same arrays
 */
void simulate_lti(size_t n, size_t p, const real_t a[n][n], const real_t x[n], const real_t b[n][p], const real_t u[p], real_t res[n]);

void simulate_affine(size_t n, size_t p, const real_t a[n][n], const real_t x[n], const real_t b[n][p], const real_t u[p], const real_t d[n], real_t res[n]);
