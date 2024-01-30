#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Simulate.h"
#include "Types.h"
#include "Utils.h"
#include "Csv.h"
#include "SdqpLmpcVarying.h"

static int lmpc_varying_simulate(const char* input_dir, const char* output_dir, size_t N, size_t simulation_timesteps) {
    char path[128];

    // Parse system dimensions
    sprintf(path, "%s/A.csv", input_dir);
    size_t n_x = csv_parse_matrix_width(path);
    sprintf(path, "%s/B.csv", input_dir);
    size_t n_u = csv_parse_matrix_width(path);
    sprintf(path, "%s/C.csv", input_dir);
    size_t n_y = csv_parse_matrix_height(path);
    sprintf(path, "%s/Lt.csv", input_dir);
    size_t n_t = csv_parse_matrix_height(path);

    // Allocate memory for system matrices and vectors
	real_t *Q = (real_t*)malloc(sizeof(real_t)*n_x*n_x);
	real_t *S = (real_t*)malloc(sizeof(real_t)*n_x*n_x);
	real_t *R = (real_t*)malloc(sizeof(real_t)*n_u*n_u);
	real_t *fx = (real_t*)malloc(sizeof(real_t)*n_x);
	real_t *fu = (real_t*)malloc(sizeof(real_t)*n_u);

	real_t *A = (real_t*)malloc(sizeof(real_t)*N*n_x*n_x);
	real_t *B = (real_t*)malloc(sizeof(real_t)*N*n_x*n_u);
	real_t *C = (real_t*)malloc(sizeof(real_t)*n_y*n_x);

    real_t *y_min = (real_t*)malloc(sizeof(real_t)*n_y);
    real_t *y_max = (real_t*)malloc(sizeof(real_t)*n_y);
	real_t *Lt = (real_t*)malloc(sizeof(real_t)*n_t*n_x);
	real_t *lt = (real_t*)malloc(sizeof(real_t)*n_t);
    real_t *u_min = (real_t*)malloc(sizeof(real_t)*n_u);
    real_t *u_max = (real_t*)malloc(sizeof(real_t)*n_u);

    real_t *x = (real_t*)malloc(sizeof(real_t)*N*n_x);
    real_t *u = (real_t*)malloc(sizeof(real_t)*N*n_u);

    real_t *xout = (real_t*)malloc(sizeof(real_t)*(simulation_timesteps+1)*n_x);
    real_t *uout = (real_t*)malloc(sizeof(real_t)*simulation_timesteps*n_u);

    // Parse system matrices and vectors
    sprintf(path, "%s/Q.csv", input_dir);
    if (csv_parse_matrix(path, n_x, n_x, CAST_2D_VLA(Q, n_x))) { 
        printf("Error while parsing input matrix Q.\n");
        return 1;
    }
    sprintf(path, "%s/S.csv", input_dir);
    if (csv_parse_matrix(path, n_x, n_x, CAST_2D_VLA(S, n_x))) { 
        printf("Error while parsing input matrix S.\n");
        return 1;
    }
    sprintf(path, "%s/R.csv", input_dir);
    if (csv_parse_matrix(path, n_u, n_u, CAST_2D_VLA(R, n_u))) { 
        printf("Error while parsing input matrix R.\n");
        return 1;
    }
    sprintf(path, "%s/fx.csv", input_dir);
    if (csv_parse_vector(path, n_x, fx)) { 
        printf("Error while parsing input vector fx.\n");
        return 1;
    }
    sprintf(path, "%s/fu.csv", input_dir);
    if (csv_parse_vector(path, n_u, fu)) { 
        printf("Error while parsing input vector fu.\n");
        return 1;
    }

    for (size_t i = 0; i < N; ++i) {
        sprintf(path, "%s/A.csv", input_dir);
        if (csv_parse_matrix(path, n_x, n_x, CAST_2D_VLA(&A[i*n_x*n_x], n_x))) { 
            printf("Error while parsing input matrix A.\n");
            return 1;
        }
        sprintf(path, "%s/B.csv", input_dir);
        if (csv_parse_matrix(path, n_x, n_u, CAST_2D_VLA(&B[i*n_x*n_u], n_u))) { 
            printf("Error while parsing input matrix B.\n");
            return 1;
        }
    }
    sprintf(path, "%s/C.csv", input_dir);
    if (csv_parse_matrix(path, n_y, n_x, CAST_2D_VLA(C, n_x))) { 
        printf("Error while parsing input matrix C.\n");
        return 1;
    }

    sprintf(path, "%s/y_min.csv", input_dir);
    if (csv_parse_vector(path, n_y, y_min)) {
        printf("Error while parsing input vector y_min.\n");
        return 1;
    }
    sprintf(path, "%s/y_max.csv", input_dir);
    if (csv_parse_vector(path, n_y, y_max)) {
        printf("Error while parsing input vector y_max.\n");
        return 1;
    }
    sprintf(path, "%s/Lt.csv", input_dir);
    if (csv_parse_matrix(path, n_t, n_x, CAST_2D_VLA(Lt, n_x))) { 
        printf("Error while parsing input matrix Lt.\n");
        return 1;
    }
    sprintf(path, "%s/lt.csv", input_dir);
    if (csv_parse_vector(path, n_t, lt)) { 
        printf("Error while parsing input vector lt.\n");
        return 1;
    }
    sprintf(path, "%s/u_min.csv", input_dir);
    if (csv_parse_vector(path, n_u, u_min)) {
        printf("Error while parsing input vector u_min.\n");
        return 1;
    }
    sprintf(path, "%s/u_max.csv", input_dir);
    if (csv_parse_vector(path, n_u, u_max)) {
        printf("Error while parsing input vector u_max.\n");
        return 1;
    }

    sprintf(path, "%s/x0.csv", input_dir);
    if (csv_parse_vector(path, n_x, xout)) { 
        printf("Error while parsing input vector x0.\n");
        return 1;
    }

    // Initialize solver
    sdqp_lmpc_varying_init(n_x, n_u, n_y, n_t, N, 
            CAST_CONST_2D_VLA(Q, n_x), CAST_CONST_2D_VLA(S, n_x), CAST_CONST_2D_VLA(R, n_u), fx, fu, 
            CAST_CONST_2D_VLA(C, n_x), 
            y_min, y_max, CAST_CONST_2D_VLA(Lt, n_x), lt, u_min, u_max);

    // Simulate
    for (size_t i = 0; i < simulation_timesteps; ++i) {
        int err = sdqp_lmpc_varying_solve(n_x, n_u, N, 
                CAST_CONST_3D_VLA(A, n_x, n_x), CAST_CONST_3D_VLA(B, n_x, n_u), &xout[i*n_x], 
                x, u);
        if (err) {
            printf("Error while solving: %d\n", err);
            return 1;
        }
        simulate_lti(n_x, n_u, CAST_CONST_2D_VLA(&A[0*n_x*n_x], n_x), &xout[i*n_x], CAST_CONST_2D_VLA(&B[0*n_x*n_u], n_u), u, &xout[(i+1)*n_x]);
        memcpy(&uout[i*n_u], u, sizeof(real_t)*n_u);
    }

    // Save output
    sprintf(path, "%s/xout.csv", output_dir);
    if (csv_save_matrix(path, simulation_timesteps+1, n_x, CAST_CONST_2D_VLA(xout, n_x)) < 0) {
        printf("Error while saving xout.\n");
        return 1;
    }
    sprintf(path, "%s/uout.csv", output_dir);
    if (csv_save_matrix(path, simulation_timesteps, n_u, CAST_CONST_2D_VLA(uout, n_u)) < 0) {
        printf("Error while saving uout.\n");
        return 1;
    }

    // Cleanup
    sdqp_lmpc_varying_cleanup();

	free(Q);
	free(S);
	free(R);
	free(fx);
	free(fu);

	free(A);
	free(B);
	free(C);

    free(y_min);
    free(y_max);
	free(Lt);
	free(lt);
    free(u_min);
    free(u_max);

	free(x);
	free(u);

    free(xout);
    free(uout);

    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        return 1;
    } else {
        const char *input_dir = argv[1];
        const char *output_dir = argv[2];
        size_t N = (size_t)atoi(argv[3]);
        size_t simulation_timesteps = (size_t)atoi(argv[4]);
        return lmpc_varying_simulate(input_dir, output_dir, N, simulation_timesteps);
    }
    return 1;
}
