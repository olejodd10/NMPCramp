#include <stdlib.h>
#include <stddef.h>

#include "LmpcConstant.h"

int main(int argc, char *argv[]) {
    if (argc < 5) {
        return 1;
    } else {
        const char *input_dir = argv[1];
        const char *output_dir = argv[2];
        size_t N = (size_t)atoi(argv[3]);
        size_t simulation_timesteps = (size_t)atoi(argv[4]);
        return lmpc_constant_simulate(input_dir, output_dir, N, simulation_timesteps);
    }
    return 1;
}
