#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Simulate.h"
#include "Types.h"
#include "Utils.h"
#include "Csv.h"
#include "Timer.h"
#include "LtvMpc.h"

static int ltv_mpc_simulate(const char* input_dir, const char* output_dir, size_t N, size_t simulation_timesteps) {
    char path[128];

    // Parse system dimensions
    sprintf(path, "%s/A.csv", input_dir);
    size_t n_x = csv_parse_matrix_width(path);
    sprintf(path, "%s/B.csv", input_dir);
    size_t n_u = csv_parse_matrix_width(path);
    sprintf(path, "%s/C.csv", input_dir);
    size_t n_y = csv_parse_matrix_height(path);
    sprintf(path, "%s/_Lt.csv", input_dir);
    size_t n_t = csv_parse_matrix_height(path);

    // Allocate memory for system matrices and vectors
	real_t *Q = (real_t*)malloc(sizeof(real_t)*n_x*n_x);
	real_t *S = (real_t*)malloc(sizeof(real_t)*n_x*n_x);
	real_t *R = (real_t*)malloc(sizeof(real_t)*n_u*n_u);
	real_t *fx = (real_t*)malloc(sizeof(real_t)*N*n_x);
	real_t *fu = (real_t*)malloc(sizeof(real_t)*N*n_u);

	real_t *A = (real_t*)malloc(sizeof(real_t)*N*n_x*n_x);
	real_t *B = (real_t*)malloc(sizeof(real_t)*N*n_x*n_u);
	real_t *d = (real_t*)malloc(sizeof(real_t)*N*n_x);
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
    real_t *tout = (real_t*)malloc(sizeof(real_t)*simulation_timesteps);

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
    if (csv_parse_vector(path, n_x, &fx[0*n_x])) { 
        printf("Error while parsing input vector fx.\n");
        return 1;
    }
    sprintf(path, "%s/fu.csv", input_dir);
    if (csv_parse_vector(path, n_u, &fu[0*n_u])) { 
        printf("Error while parsing input vector fu.\n");
        return 1;
    }
    for (size_t i = 1; i < N; ++i) {
        memcpy(&fx[i*n_x], &fx[0*n_x], sizeof(real_t)*n_x);
        memcpy(&fu[i*n_u], &fu[0*n_u], sizeof(real_t)*n_u);
    }

    sprintf(path, "%s/A.csv", input_dir);
    if (csv_parse_matrix(path, n_x, n_x, CAST_2D_VLA(&A[0*n_x*n_x], n_x))) { 
        printf("Error while parsing input matrix A.\n");
        return 1;
    }
    sprintf(path, "%s/B.csv", input_dir);
    if (csv_parse_matrix(path, n_x, n_u, CAST_2D_VLA(&B[0*n_x*n_u], n_u))) { 
        printf("Error while parsing input matrix B.\n");
        return 1;
    }
    for (size_t i = 1; i < N; ++i) {
        memcpy(&A[i*n_x*n_x], &A[0*n_x*n_x], sizeof(real_t)*n_x*n_x);
        memcpy(&B[i*n_x*n_u], &B[0*n_x*n_u], sizeof(real_t)*n_x*n_u);
    }
    memset(d, 0, sizeof(real_t)*N*n_x);
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
    sprintf(path, "%s/_Lt.csv", input_dir);
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
    ltv_mpc_init(n_x, n_u, n_y, n_t, N, 
            CAST_CONST_2D_VLA(Q, n_x), CAST_CONST_2D_VLA(S, n_x), CAST_CONST_2D_VLA(R, n_u), 
            CAST_CONST_2D_VLA(fx, n_x), CAST_CONST_2D_VLA(fu, n_u),
            CAST_CONST_2D_VLA(C, n_x), 
            y_min, y_max, CAST_CONST_2D_VLA(Lt, n_x), lt, u_min, u_max);

    // Simulate
    long long simulation_time_ns = 0;
    for (size_t i = 0; i < simulation_timesteps; ++i) {
        timer_reset();
        int err = ltv_mpc_solve(n_x, n_u, N, 
                CAST_CONST_3D_VLA(A, n_x, n_x), CAST_CONST_3D_VLA(B, n_x, n_u), CAST_CONST_2D_VLA(d, n_x), &xout[i*n_x], 
                CAST_2D_VLA(x, n_x), CAST_2D_VLA(u, n_u));
        if (err) {
            printf("Error while solving: %d\n", err);
            return 1;
        }
        long timestep_elapsed_ns = timer_elapsed_ns();
        tout[i] = ((real_t)timestep_elapsed_ns)/1000.0;
        simulation_time_ns += (long long)timestep_elapsed_ns;
        simulate_lti(n_x, n_u, CAST_CONST_2D_VLA(&A[0*n_x*n_x], n_x), &xout[i*n_x], CAST_CONST_2D_VLA(&B[0*n_x*n_u], n_u), u, &xout[(i+1)*n_x]);
        memcpy(&uout[i*n_u], u, sizeof(real_t)*n_u);
    }

    // Timer output
    printf("%ld timesteps with horizon %ld finished in %lld ms\n", simulation_timesteps, N, simulation_time_ns/1000000);

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
    sprintf(path, "%s/toutN%ld.csv", output_dir, N);
    if (csv_save_vector(path, simulation_timesteps, tout) < 0) {
        printf("Error while saving tout.\n");
        return 1;
    }

    // Cleanup
    ltv_mpc_cleanup();

	free(Q);
	free(S);
	free(R);
	free(fx);
	free(fu);

	free(A);
	free(B);
    free(d);
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
    free(tout);

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
        return ltv_mpc_simulate(input_dir, output_dir, N, simulation_timesteps);
    }
    return 1;
}
