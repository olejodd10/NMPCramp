#include "Extrapolate.h"

#include <stddef.h>
#include <math.h>

#include "Types.h"

void extrapolate_sine(size_t N, real_t amp, real_t freq, real_t Ts, real_t phase, real_t offset, real_t sine[N]) {
    real_t omega = 2.0*M_PI*freq*Ts;
    for (size_t i = 0; i < N; ++i) {
        sine[i] = amp * sin(omega * (real_t)i + phase) + offset;
    }
}
