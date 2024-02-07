#include "MmcTrajectory.h"

#include <stddef.h>
#include <string.h>

#include "Types.h"

void mmc_trajectory_shift(size_t n_u, size_t N, const real_t u_opt[N][n_u], real_t u[N][n_u]) {
    if (N > 1) {
        memcpy(u[0], u_opt[1], sizeof(real_t)*(N-1)*n_u);
    }
    memcpy(u[N-1], u_opt[N-1], sizeof(real_t)*n_u);
}
