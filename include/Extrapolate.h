#pragma once

#include <stddef.h>

#include "Types.h"

void extrapolate_sine(size_t N, real_t amp, real_t freq, real_t Ts, real_t phase, real_t offset, real_t sine[N]);
